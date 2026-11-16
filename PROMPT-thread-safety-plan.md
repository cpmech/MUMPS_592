# Prompt — design a thread-safety plan for MUMPS 5.9.2

You are a senior numerical-software engineer with Fortran, OpenMP, MPI and sparse-direct-solver
experience. Your task is to **design a plan**, not to implement it: make MUMPS 5.9.2 safe to call
concurrently from several threads of one process, well enough that a client library can stop
serializing every MUMPS call onto a single dedicated thread. Read the whole brief before starting,
then produce the deliverable in §9.

## 1. Why this work exists

The stack is `pmsim` (Rust FEM) → `russell_sparse` (Rust wrapper) → a small C shim
(`russell_sparse/c_code/interface_mumps.c`) → `dmumps_c` / `zmumps_c` from a locally built MUMPS
5.9.2 (`/usr/local/lib/mumps`, MKL BLAS, sequential `libmpiseq`, `PLAT=_cpmech`, static and shared,
s/d/c/z).

- Typical runs: 192 dof up to ~10^5, one process, MUMPS's own OpenMP threads enabled
  (`ICNTL(16) = 0`, `ICNTL(48) = 1`). Test suites and a GUI call solvers from several threads.
- Today every MUMPS FFI call (`JOB = -1/1/2`, `JOB = 3`, `JOB = -2`) is dispatched to one
  process-wide `mumps-worker` thread in `russell_sparse/src/mumps_lock.rs`; that serializes all
  MUMPS work and keeps one persistent OpenMP master.
- Measured cost of that workaround (do not re-derive; use as the baseline to beat): ≈7 µs per MUMPS
  FFI call; +2 % on a 192-dof benchmark at `OMP_NUM_THREADS=1`, +17–22 % on the same benchmark at
  `OMP=8`, below noise at 2187 dof. It also forbids two instances from solving at the same time.
- **Goal of the plan:** the library itself provides the guarantee, so the worker thread and all
  client-side locking can be deleted from `russell_sparse`.

## 2. The defect, as measured (confirm, cite, extend — do not re-attribute)

Two failure classes, observed on this workstation (32 logical CPUs, Intel `libiomp5`): class 1 on
MUMPS 5.9.0, class 2 on 5.9.2.

1. **Concurrent entry corrupts state.** Two threads inside `dmumps_c` at once produce `SIGSEGV`
   inside `libdmumps_cpmech.so` / `libiomp5.so` (`__kmp_reap_worker`) and MUMPS's own
   `PB allocation in MUMPS_LOAD_INIT`.
2. **A client mutex does not fix it.** With one MUMPS call at a time, the crash persists once
   `N_caller_threads × OMP_NUM_THREADS > 128`: `OMP=10,12` pass 6/6; `OMP=13` fails 1/6; `OMP=14`
   fails 3/6; `OMP=16` fails 6/6; `OMP=32` fails 4/6; `KMP_ALL_THREADPRIVATE=512` still fails 0/10 at
   `OMP=16`. `gdb` + `objdump` place the fault in `DMUMPS_SOL_L0OMP_R` (`src/dsol_omp_m.F`), whose
   module `DMUMPS_SOL_L0OMP_M` uses `!$OMP THREADPRIVATE(THREAD_ID, IS_ROOT_OF_L0_OMP)`
   (`src/mumps_l0_omp_m.F:20`): `__kmpc_threadprivate_cached(gtid=145)` runs past libiomp5's
   `__kmp_default_tp_capacity = 128` (`4 × 32` CPUs), together with caller-thread OpenMP team reaping
   in `__kmp_reap_worker`.

So the core problem is **ownership of the OpenMP master thread and its teams**, not just mutual
exclusion: a per-process lock inside MUMPS is *not* a solution by itself. The single persistent
worker is the current workaround (20/20 clean at `OMP=1`, `16` and `32`).

## 3. What "thread-safe" must mean (acceptance criteria)

Define the guarantee precisely, and say which of these the design reaches.

- N threads (≥8), each owning its own instance, any interleaving of `JOB=-1 → 1/2 → 3 → -2`, with no
  client lock, no crash, no corruption, and each solution equal to the serial run (bit-identical at
  `OMP=1`; within solver round-off otherwise).
- MUMPS's own OpenMP parallelism stays enabled (`ICNTL(16)=0`, `ICNTL(48)=1`) and repeated
  create/solve/destroy churn does not exhaust the OpenMP runtime (`__kmp_default_tp_capacity`).
- At least s/d/c/z and `SYM=0/2`; the four arithmetics must not diverge in behavior.
- The sequential `libmpiseq` build is the target; MPI builds must not regress.
- A documented boundary: exactly which features are **not** covered (OOC, save/restore, BLR,
  `ICNTL(17)` MPI-to-OpenMP, METIS/Scotch orderings), with the failure mode named — an explicit
  boundary is acceptable, a silent one is not.
- Ideally independent instances run **concurrently**; if the design only makes concurrent calls
  *safe* (still serialized), say so and quantify what concurrency costs.

## 4. Investigate — starting map (verify every line yourself; this is not a conclusion)

There is no named `COMMON` block in `src/*.F`; mutable state lives in module data, in explicit `SAVE`
statements inside routines (e.g. `src/dsol_aux.F`), and in the C interface. A first scan found
candidates in ~83 files.
Start with these, classify each as per-instance / process-global-read-only / process-global-mutable,
and name the callers and the failure if two instances share it:

- OpenMP and thread state: `src/mumps_l0_omp_m.F` (`THREADPRIVATE`, `NB_CORE_PER_THREAD*`),
  `src/omp_tps_common_m.F`, `src/sol_omp_common_m.F`, `src/fac_future_niv2_mod.F` (`FUTURE_NIV2`).
- Communication buffers: `src/d{m,s,c,z}mumps_comm_buffer.F` (`BUF_LMAX_ARRAY`, `BUF_MAX_ARRAY`),
  `src/mumps_comm_buffer_common.F` (`SIZEofINT`, `SIZEofREAL`, `BUF_MYID`, `SIZE_RBUF_BYTES`).
- Load balancing and memory: `src/mumps_load.F` (`LOAD_FLOPS`, `BUF_LOAD_RECV`, …),
  `src/mumps_memory_mod.F` (`ISIZE`, `I8SIZE`, …), `src/mumps_pivnul_mod.F` (`PIVNUL_LIST`).
- Out-of-core: `src/d…ooc.F`, `src/d…ooc_buffer.F` (`PANEL_FLAG`, buffer counters, `SAVE`),
  `src/mumps_ooc_common.F` (`TYPEF_L`, `TYPEF_U`, `TYPEF_CB`, sequences).
- Scratch with explicit `SAVE` in routines: `src/dsol_aux.F` (`SAVE ITER, J, JLAST, JUMP`, …; same
  in `ssol_aux.F`, `csol_aux.F`, `zsol_aux.F`), `src/dstatic_ptr_m.F` (`DMUMPS_TMP_PTR`),
  `src/lr_common.F` (precision settings), `src/mumps_static_mapping.F` (`cv_*`), `src/ana_blk_m.F`,
  `src/ana_*par.F`.
- The C interface: `src/mumps_c.c` — six file-scope static pointers (lines ~244–249:
  `*_COLSCA_STATIC`, `*_ROWSCA_STATIC`, `*_COLSCA_LOC_STATIC`, `*_ROWSCA_LOC_STATIC`,
  `*_RHSINTR_STATIC`, `*_SINGULAR_VALUES_STATIC`) plus `MUMPS_SET_TMP_PTR` → the Fortran static in
  `src/dstatic_ptr_m.F`.
- Libraries below MUMPS: `libseq` (verify static state in `mpic.c`/`mpi.f`), `PORD/lib` (appears
  argument-driven — confirm), METIS/Scotch thread safety, MKL BLAS/MKL threading.
- Existing OpenMP locks to reuse or extend: `OMP_INIT_LOCK`/`OMP_SET_LOCK` in `src/dmumps_ooc.F`
  (~lines 151, 3018, 3076) and `LOCK_FOR_SCATTER` in `src/dsol_fwd_aux.F:900,915`.

Also read what the library promises today: the user guide (extracted text section around “Several
instances of MUMPS can be handled simultaneously”) and any statement about thread safety; the plan
must correct the documentation if the guarantee changes.

## 5. Design options to evaluate (at least these)

a. Process-wide lock inside MUMPS — the control, not a solution (fails class 2);
b. MUMPS-owned single master thread, dispatched inside `dmumps_c` — moves today's workaround into the
   library: zero client setup, still serialized;
c. Per-instance master thread (one worker per instance) — allows concurrent instances; analyze team
   count, `tp_capacity` and teardown;
d. True reentrancy — every mutable global moves into the instance (or is passed down); no dedicated
   threads at all;
e. Hybrid — reentrant data state for a supported subset (sequential, in-core, assembled matrix, no
   BLR/save-restore/OOC) plus MUMPS-owned OpenMP team(s);
f. Scope-restricted documented mode — `ICNTL(16)=1` per instance plus full reentrancy (cheapest
   correct option, worst performance for large problems).

For each: what changes and where (`file:line`), ABI impact, effort in person-weeks (not days),
performance model (per-call cost, solve throughput, memory, thread count), risk, and the experiment
that would falsify it. Recommend one with a fallback, and split it into phases, each with a gate and a
stop condition.

## 6. Constraints

- **ABI.** Do not change `DMUMPS_STRUC_C` layout or existing entry-point signatures; the header in
  `/usr/local/include/mumps` must keep working for the existing client shim. A new capability must be
  detectable from a client (e.g., an extra exported symbol queried at build time), so a stock library
  and the patched one are distinguishable.
- **Build.** `~/01-Code/rust/russell/zscripts/{arch,debian,macos,windows}-compile-mumps.bash` download
  the official tarball (`VERSION="5.9.2"`), copy `zscripts/makefiles-mumps/MakefileMKL.inc` (or
  `Makefile.inc`) over `Makefile.inc`, build `d`, `dshared`, `z`, `zshared`, and install to
  `/usr/local/{lib,include}/mumps`. The plan must say where the patch lives (a maintained fork vs. a
  patch series applied after `tar xzf`), how it is versioned (`MUMPS_SCRIPT_VERSION`), and how a stock
  tarball still builds unchanged.
- **License.** CeCILL-C: modifications may be redistributed if the license and attribution are kept;
  state the obligation and where the modified files are marked.
- **No regressions for stock users.** Removing the `mumps-worker` from `russell_sparse` must be gated
  behind the capability check; stock MUMPS keeps the worker.
- **American English**, and **every claim about MUMPS internals carries `file:line`**. No invented
  function or module names.
- **Design only.** Do not edit this tree; build probes from the `/tmp` copy that `build.bash` makes
  (see §8). Do not touch `pmsim` or `russell` — a client-side change is not the fix.

## 7. Probes worth running first (cheap, and they decide the design)

Every probe can be built and run against the installed library today; §8 covers building this tree
once the source changes.

1. A standalone reproducer, no Rust: `examples/dorival_thread_safe.F` already does it (its header
   documents the flags; without a lock it fails 12/12 runs of 2 and 8 threads with a `SIGSEGV`, an
   allocatable-array abort or `INFOG(1) = -13`, and with a lock around the whole cycle it passes
   6/6). Confirm both failure classes against the installed stock library, record exit codes,
   `dmesg`/journal lines and the faulting offset.
2. Does `ICNTL(48)=0` (no L0 tree parallelism) remove the `OMP ≥ 13` failure with many caller threads?
   Does `ICNTL(16)=1`? Does `ICNTL(16)=k` fixed per instance change the threshold?
3. Is the threshold driven by instances, caller threads, or `N_threads × OMP`? Does a caller thread
   that never enters MUMPS consume team capacity?
4. What static state does `libmpiseq` keep, and is `dmumps_c` reentrant when only one instance is
   active per thread?
5. Which of the six C statics in `src/mumps_c.c` are touched per `dmumps_c` call, and are the
   `set_tmp_ptr`/`assign_colsca` paths exercised by an assembled-matrix run?
6. Scope of the existing OpenMP locks (`dmumps_ooc.F`, `dsol_fwd_aux.F`): per instance or per process?

## 8. Building this tree

**A compiled 5.9.2 is already installed**, in `/usr/local/lib/mumps` and `/usr/local/include/mumps`
(built from this upstream release with the same build files: `PLAT = _cpmech`, s/d/c/z, static and
shared). The reproducer of §7 can be written, compiled and run today against that library, without
building anything first.

**`build.bash` in this directory compiles the checkout** for when the source changes:

```bash
./build.bash                         # Intel oneAPI (ifx/icx) + MKL, the workstation default
./build.bash 0                       # gfortran/Clang + system BLAS and LAPACK
PREFIX=/tmp/mumps-probe BUILD_DIR=/tmp/mumps-build LDCONFIG=0 ./build.bash
```

It copies the working tree to `${BUILD_DIR:-/tmp}/MUMPS_<version>`, so uncommitted edits are what
gets built and the checkout stays clean; selects `zscripts/MakefileMKL.inc` or
`zscripts/Makefile.inc` as `Makefile.inc`; builds `d`, `dshared`, `z`, `zshared` with `make clean`
between them; and installs the libraries and the headers under `${PREFIX:-/usr/local}` (registering
the directory with `ldconfig` only for `/usr/local`). `TARGETS` selects a subset of the four,
`JOBS` sets the make parallelism (default 8), and `KEEP_BUILD=1` leaves the `/tmp` copy for
inspection after a successful install.

Points that matter for a probe:

- **Do not install a probe build over `/usr/local`.** Use a scratch prefix, link against
  `-L/tmp/mumps-probe/lib/mumps`, and run with `LD_LIBRARY_PATH=/tmp/mumps-probe/lib/mumps`. The
  stock install is the reference the client stack (`russell_sparse`, `pmsim`) links against, and
  leaving it untouched means a crash can be reproduced against the unmodified library at any time.
- **Headers must come from a finished build:** every client header includes the generated
  `mumps_int_def.h`, and `make clean` deletes it. `build.bash` copies the headers before its last
  clean; a manual `make` has to do the same.
- **Rebuilds are cheap and always clean.** The script copies the tree anew, so there are no stale
  objects and no incremental-make risk (MUMPS's Fortran modules make a partial rebuild unsafe). With
  `make -j8` the four variants take about 75 s here, and `TARGETS="d dshared" ./build.bash` — the loop
  for a probe that calls `dmumps_c` — takes about 40 s; `JOBS` raises the parallelism.
- **A C probe links the shared libraries**: the static archives are Intel Fortran objects, and a
  plain `gcc` link fails on the Fortran runtime (`for_write_seq_lis`, …). Build `dshared` and link
  the `.so`, or build the probe with `ifx` if it is to use the archives.

## 9. Deliverable

One Markdown file, `PLAN-thread-safety.md`, in this directory (`~/01-Code/fortran/MUMPS_5.9.2/`),
with:

1. `Status:` line (state, date, what is local).
2. Goal and scope, including the feature boundary.
3. What exists (library build, client usage, the two failure classes with the numbers above).
4. Root cause and fault model (`file:line`, OpenMP runtime behavior).
5. State-inventory table: symbol, file:line, instance/global, callers, failure if shared.
6. Design options as numbered decisions, each with a recommendation and its cost; the owner answers by
   number.
7. Recommended design in phases, each with a gate, a stop condition, and the test that proves it.
8. ABI and capability detection; what `russell_sparse` must do to drop the worker.
9. Test plan with exact commands: the standalone C reproducer, the churn matrix
   (threads × `OMP_NUM_THREADS` × instances × arithmetics), the `russell_sparse` stress test
   (`russell_sparse/tests/test_mumps.rs::test_mumps_concurrent_stress`, 10 threads × 150 iterations)
   as the integration gate, and 20/20 clean full-suite runs.
10. Performance protocol and the numbers in §1 as the baseline to beat (hyperfine, `taskset`,
    `OMP=1` for bit identity and CSV diff, `OMP=8` for the threaded penalty).
11. Build and packaging (patch location, version string, how CI and crates.io users are unaffected).
12. Risks, open questions, and a Kanban (`In progress` / `Backlog` / `Done`) of the work the plan
    implies.

Then give a ≤300-word summary in chat: the recommended design, its effort, the first three things to
verify, and any question that is genuinely the owner's to answer (numbered).

## 10. Inputs

- MUMPS source: `~/01-Code/fortran/MUMPS_5.9.2` — git repository on `main`, the upstream 5.9.2
  tarball as the first commit plus `build.bash` and `zscripts/`; 25 MB, 492 files in `src/`;
  `doc/userguide_5.9.2.pdf`; `libseq/`; `PORD/`.
- Installed library and headers: `/usr/local/lib/mumps`, `/usr/local/include/mumps` — the compiled
  stock 5.9.2 the client stack links against (§8).
- Client: `~/01-Code/rust/russell` — `russell_sparse/src/mumps_lock.rs` (the worker, 208 lines),
  `russell_sparse/c_code/interface_mumps.c` and `interface_complex_mumps.c`,
  `russell_sparse/tests/test_mumps.rs`, `zscripts/*-compile-mumps.bash`,
  `zscripts/makefiles-mumps/MakefileMKL.inc`.
- Consumer: `~/01-Code/rust/pmsim` — `src/fem/simulator.rs` (the concurrency tests).
- Project record (read-only, useful context and the measurement details):
  `~/03-Writings/ObsidianVault/Note-Mumps-Thread-Safety.md` and
  `~/03-Writings/ObsidianVault/Plan-Mumps-Thread-Safety.md`.

## 11. Do not

- Do not present "wrap it in a mutex" as the solution — a client mutex already exists and does not
  fix failure class 2.
- Do not propose client-side (pmsim/russell) changes as the fix; the point is to remove them.
- Do not claim thread safety for features you scoped out.
- Do not write code now; a plan that survives review is the deliverable.
