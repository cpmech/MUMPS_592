# MUMPS 5.9.2 Thread-Safety Plan

## 1. Status

`Status:` **Proposed — root-cause analysis and runtime probes complete (2026-11-17); Phase 0 in progress.** Phase 0 dispatched to Kleo on branch `ts-phase0` from `main` `d7b6b1f` (brief at `/home/dorival/.hermes/profiles/ada/cache/scratch/brief-kleo-mumps-phase0.md`); Ada reviews. Local workspace state: [PLAN-thread-safety.md](file:///home/dorival/01-Code/fortran/MUMPS_5.9.2/PLAN-thread-safety.md) added on branch `main` over commit `2da3511fce243579b55b85218e52f92f9184f499`. No MUMPS source files, `pmsim` files, or `russell` files have been modified.

---

## 2. Goal, Acceptance Criteria, and Feature Boundary

### 2.1 Goal

Make MUMPS 5.9.2 (`PLAT=_cpmech`, sequential MPI stub `libmpiseq_cpmech.{a,so}`, real and complex arithmetics `s/d/c/z`, static `.a` and shared `.so`) thread-safe under concurrent calls from arbitrary OS threads of a single process, eliminating the client-side `mumps-worker` serialization thread and mutex from `russell_sparse` (`russell_sparse/src/mumps_lock.rs`) while preserving MUMPS's internal OpenMP parallelism (`ICNTL(16) = 0`, `ICNTL(48) = 1`) and ABI compatibility.

### 2.2 Acceptance Criteria

| Criterion | Exact Guarantee Required | Phase 1 (Serialized Persistent Master) | Phase 2 (Concurrent $K$-Worker Pool) |
| :--- | :--- | :--- | :--- |
| **AC-1: Concurrent entry & churn** | $N \ge 8$ OS caller threads (`pthread`, Rust `std::thread`, or OpenMP threads), each owning its own instance and running arbitrary interleaved `JOB = -1 -> 1/2 -> 3 -> -2` cycles with **no client lock**, complete with zero `SIGSEGV`, zero heap corruption, zero spurious `INFOG(1) < 0` errors, and zero team exhaustion under churn. | **Met** (serialized safely inside `libmumps_common`) | **Met** (parallel across $K$ persistent masters) |
| **AC-2: Numerical accuracy** | At `OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`, solutions are bit-identical (`max |x_ts - x_ref| == 0.0`) to stock MUMPS 5.9.2; at `OMP_NUM_THREADS > 1`, within IEEE-754 floating-point reduction round-off ($\le 10^{-13}$ relative in double precision). | **Met** | **Met** |
| **AC-3: Internal OpenMP preserved** | `ICNTL(16) = 0` (runtime thread count) and `ICNTL(48) = 1` (L0 tree-level OpenMP) remain active and functional without exhausting OpenMP runtime threadprivate tables (`__kmp_default_tp_capacity`) or crashing on caller-thread exit (`__kmp_reap_worker`). | **Met** | **Met** |
| **AC-4: Arithmetic & symmetry matrix** | Uniform behavior across `s`, `d`, `c`, `z` and all symmetry modes (`SYM = 0` unsymmetric, `SYM = 1` SPD, `SYM = 2` general symmetric). | **Met** | **Met** |
| **AC-5: Sequential & MPI scope** | Targets sequential `libmpiseq` builds (`-DMUMPS_SEQ`); standard MPI builds (`-DMUMPS_MPI`) must compile and run with zero regressions. | **Met** | **Met** |
| **AC-6: True concurrency** | Independent instances factorize and solve concurrently on distinct core subsets to achieve multi-instance wall-clock speedup. | **Not met** (serialized in library) | **Met** (up to $K$ concurrent solves) |

### 2.3 Explicit Feature Boundary

Every MUMPS 5.9.2 feature is classified below. Features outside the Phase 2 concurrent boundary are explicitly serialized or guarded—never left silently racy.

| MUMPS Feature | Control Parameters | Phase 1 Status | Phase 2 Status & Rationale / Failure Mode if Unprotected |
| :--- | :--- | :--- | :--- |
| **In-core assembled sparse solve** | `ICNTL(5)=0`, `ICNTL(18)=0`, `ICNTL(22)=0` | **Thread-safe** (serialized) | **Concurrent** (primary target for `russell_sparse` and `pmsim`). |
| **Dense / sparse / multi-RHS solves & transpose** | `ICNTL(20)=0..3`, `ICNTL(21)=0`, `ICNTL(9)=0,1` | **Thread-safe** | **Concurrent** once C/F77 statics (`src/mumps_c.c:244-249`, `src/mumps_common.c:20-26`) and `DMUMPS_TMP_PTR` (`src/[sdcz]static_ptr_m.F:16`) are removed. |
| **Scaling (automatic & user-provided)** | `ICNTL(6)=0..7`, `ICNTL(8)=-2..8` | **Thread-safe** | **Concurrent** (requires removal of `[SDCZ]MUMPS_COLSCA_STATIC` and `ROWSCA_STATIC` in `src/mumps_c.c:244-245`, since dummy length-1 scaling arrays are allocated even when `ICNTL(8)=0` at `src/dmumps_driver.F:1926-1934`). |
| **Iterative refinement & error analysis** | `ICNTL(10) != 0`, `ICNTL(11)=1,2` | **Thread-safe** | **Concurrent** once explicit `SAVE` reverse-communication variables in `[SDCZ]MUMPS_SOL_B` (`src/[sdcz]sol_aux.F:319`), `OMEGA` (`line 883`), and `LCOND` (`line 958`) move to local state structs. |
| **Ordering: AMD, AMF, QAMD** | `ICNTL(7)=0,2,6` | **Thread-safe** | **Concurrent** (pure Fortran `src/ana_orderings.F`; no module/static globals). |
| **Ordering: PORD** | `ICNTL(7)=4` | **Thread-safe** | **Concurrent** (`nm -A /usr/local/lib/mumps/libpord_cpmech.a` confirms zero `.bss`/`.data`/`COMMON` mutable globals). |
| **Ordering: METIS** | `ICNTL(7)=5` | **Thread-safe** | **Concurrent**: `nm -D /usr/lib/libmetis.so` shows zero writable globals in `libmetis.so`; `nm -D /usr/lib/libGKlib.so` shows TLS PRNG state and signal `gk_jbuf`/`gk_jbufs` (unused in library mode). |
| **Ordering: SCOTCH** | `ICNTL(7)=3` | **Out of scope** | **Out of scope**: disabled in `zscripts/MakefileMKL.inc:25` (`ORDERINGSF = -Dmetis -Dpord`). |
| **Schur complement & null pivot detection** | `ICNTL(19)>0`, `ICNTL(24)=1` | **Thread-safe** | **Concurrent** after eliminating `DMUMPS_TMP_PTR` at `src/dmumps_driver.F:1662-1664, 2022-2025`. (`MUMPS_PIVNUL_MOD` data lives on `dfac_driver.F:185` stack). |
| **Block Low-Rank (`BLR`)** | `ICNTL(35)>0` | **Thread-safe** (serialized) | **Concurrent** if `BLR_ARRAY` (`src/[sdcz]mumps_lr_data_m.F:208-209`), `MUMPS_LR_COMMON` (`src/lr_common.F:20-42`), `MUMPS_LR_STATS` (`src/lr_stats.F:30-128`), and `mumps_sol_omp_memory_manager.c:32-34` move into `idintr` / TLS; otherwise **serialized** with a BLR lock (`arena` double-free / `BLR_ARRAY` overwrite). |
| **Out-of-Core (`OOC`)** | `ICNTL(22)=1` | **Serialized (single active OOC instance)** | **Serialized via OOC mutex on Worker 0** (out of scope for concurrent execution): global file handles/buffers in `src/mumps_ooc_common.F:18-34`, `src/[sdcz]mumps_ooc.F:19-55`, `src/[sdcz]mumps_ooc_buffer.F:20-35`, `src/[sdcz]mumps_sol_es.F:27`, and `src/mumps_io*.c` persist across `JOB = 2` and `JOB = 3`; concurrent OOC corrupts temp factor files and segfaults `io_queue`. |
| **Save / Restore** | `JOB=7,8,-3` | **Serialized** | **Serialized via OOC mutex** (shares OOC state checks in `src/dmumps_save_restore.F:23`). |
| **MPI-to-OpenMP & Parallel Analysis** | `ICNTL(17)>0`, `ICNTL(28)=2` | **Out of scope** (MPI-only) | **Out of scope**: `libmpiseq` runs `NPROCS=1` (`ICNTL(28)=2` falls back to sequential at `src/dana_driver.F:847-865`); `[SDCZ]MUMPS_SEND_BUF` `SAVE` state at `src/[sdcz]ana_aux_par.F:3103-3106` is unused in sequential mode. |
| **Configuration file** | `DMUMPS_CONFIG_FROM_FILE` | **Fixed in Phase 0** | Fixed by removing implicit `SAVE` from initialized locals at `src/dmumps_config_file.F:26-28`. |
| **GPU offload** | `-DUSE_GPU`, `ICNTL(51)>0` | **Out of scope** | **Out of scope**: not compiled (`src/[sdcz]mumps_gpu.c:32-38` global cuBLAS handle). |

### 2.4 Required User Guide Documentation Corrections

Review of `doc/userguide_5.9.2.pdf` reveals two documentation gaps that the patched release notes / manual errata must correct:
1. **Section 1 ("Several instances of MUMPS can be handled simultaneously", `doc/userguide_5.9.2.pdf`, `pdftotext` line 402):** Stock 5.9.2 supports multiple *allocated* instances in memory only when calls from host code are serialized on one OS thread (or under a lock keeping `N_callers * OMP_NUM_THREADS <= 128`). Under this plan, MUMPS guarantees simultaneous handling from **arbitrary OS threads**.
2. **Section 3.14 ("Multithreaded parallelism", `pdftotext` line 1192) & Section 5.22 (`ICNTL(16)`, `pdftotext` line 4760) & Section 5.26 (`ICNTL(48)`, `pdftotext` line 265):** The manual currently contains zero discussion of thread safety or `libiomp5` team ownership. It must document: (a) in Phase 1, concurrent C/Fortran calls are internally dispatched to MUMPS's persistent OpenMP master worker; (b) in Phase 2, up to `MUMPS_NUM_WORKERS` instances run concurrently while OOC (`ICNTL(22)=1`) remains serialized; and (c) `ICNTL(16)` is restored upon return from `dmumps` (`src/dmumps_driver.F:770, 2348-2350`).

---

## 3. What Exists (Build, Client Stack, Baselines, and Observed Failures)

### 3.1 Current Workstation Build & Installation

- **Repository:** `/home/dorival/01-Code/fortran/MUMPS_5.9.2` (upstream MUMPS 5.9.2 from `VERSION:1-2` [`Fri Nov 13 10:00:55 UTC 2026`] plus [build.bash](file:///home/dorival/01-Code/fortran/MUMPS_5.9.2/build.bash) and `zscripts/`).
- **Default toolchain (`zscripts/MakefileMKL.inc:7-70`, selected by `build.bash:38,123`):**
  - `PLAT = _cpmech` (`line 7`)
  - Compilers: `CC = icx` (`line 45`), `FC = ifx` (`line 46`), `FL = ifx` (`line 47`)
  - Flags: `OPTF = -O -nofor-main -qopenmp -DGEMMT_AVAILABLE` (`line 63`), `OPTL = -O -nofor-main -qopenmp` (`line 64`), `OPTC = -O -qopenmp` (`line 65`)
  - Sequential MPI & BLAS/LAPACK: `LIBSEQ = $(LAPACK) -L$(topdir)/libseq -lmpiseq$(PLAT)` (`line 54`), `LAPACK = -L$(MKLROOT) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core` (`line 51`), `LIBOTHERS = -lpthread` (`line 57`), orderings `-Dmetis -Dpord` (`line 25`).
- **Build & install (`build.bash:50, 148-178`):** Copies working tree to `${BUILD_DIR:-/tmp}/MUMPS_5.9.2`, builds `TARGETS="d dshared z zshared"` (which build all four arithmetics `s, d, c, z` per `src/Makefile:16-19`), and installs static `.a` and shared `.so` libraries into `/usr/local/lib/mumps` and headers (`include/*.h`, including generated `mumps_int_def.h`) into `/usr/local/include/mumps`.

### 3.2 Client Stack & Baseline Cost of `mumps-worker`

- **Call chain:** `pmsim` (Rayon / `cargo test` multi-threaded harness) $\rightarrow$ `russell_sparse` $\rightarrow$ `interface_mumps.c` / `interface_complex_mumps.c` $\rightarrow$ `dmumps_c` / `zmumps_c` (`JOB = -1, 1, 2, 3, -2`).
- **Current workaround:** `russell_sparse/src/mumps_lock.rs` serializes every MUMPS FFI call onto a single dedicated `mumps-worker` OS thread over `std::sync::mpsc::sync_channel(0)`.
- **Measured baseline performance penalty (from prompt §1 on this 32-CPU workstation):**
  - Round-trip dispatch latency per MUMPS FFI call: **≈ 7 µs**.
  - `cook_elastic_2d_fem` (192 dof, 20 steps, `OMP_NUM_THREADS=1`): **176.5 ± 2.2 ms $\rightarrow$ 180.0 ± 1.0 ms (+2.0 %)**.
  - `cook_elastic_2d_fem` (192 dof, 20 steps, `OMP_NUM_THREADS=8 MKL_NUM_THREADS=8`): **450.5 ± 22.8 ms $\rightarrow$ 552.1 ± 15.6 ms (+16.9 % to +22.5 %)** due to cross-thread wakeup and cache migration away from the Rayon caller thread.
  - `cook_elastic_3d_fem` (2187 dof, 20 steps, `OMP_NUM_THREADS=1`): **3974 ± 14 ms $\rightarrow$ 3972 ± 14 ms (0.0 %)**.
  - Concurrent multi-instance throughput: **0 % overlap** (strictly serialized).

### 3.3 The Two Observed Failure Classes

1. **Failure Class 1 — Concurrent entry into stock MUMPS corrupts process-global state:**
   - Observed in `pmsim` `test_cryer`: null-pointer dereference at `0x0` in `/usr/local/lib/mumps/libdmumps_cpmech.so` at offset `0x274520` (`dmumps_f77_+0x8b0`).
   - Observed in [dorival_thread_safe.F](file:///home/dorival/01-Code/fortran/MUMPS_5.9.2/examples/dorival_thread_safe.F) (`8` threads, `locked=0`, `n=100`): fails **6/6** runs (and `2` threads, `locked=0` fails **6/6**) with `SIGSEGV`, allocatable-array aborts, `PB allocation in MUMPS_LOAD_INIT` (`INFOG(1) = -13` at `JOB = 1/2`), or `MPI_ABORT`.
   - With `locked=2` (or `locked=1` when reusing 8 fixed OpenMP outer threads with `ICNTL(16)=1`): passes **6/6**.
2. **Failure Class 2 — Serializing calls with a mutex across dynamic OS caller threads fails once `N_caller_threads × OMP_NUM_THREADS > 128`:**
   - Observed in `russell_sparse` full suite guarded by a process-wide `Mutex<()>`:
     - `OMP_NUM_THREADS=10, 12`: **10/10 pass** (`~10 callers × 12 = 120 <= 128`).
     - `OMP_NUM_THREADS=13`: **9/10 pass (1/6 / 1/10 failure onset)** (`~10 callers × 13 = 130 > 128`).
     - `OMP_NUM_THREADS=14`: **3/6 fail**.
     - `OMP_NUM_THREADS=16`: **0/6 pass**.
     - `OMP_NUM_THREADS=32`: **2/6 pass, 4/6 fail**.
     - `KMP_ALL_THREADPRIVATE=512` at `OMP_NUM_THREADS=16`: **0/10 pass** (shifts crash from `__kmpc_threadprivate_cached` to `__kmp_reap_worker`).
   - Observed faulting sites (`gdb` + `objdump`):
     1. `DMUMPS_SOL_L0OMP_R` in `libdmumps_cpmech.so` (`src/dsol_omp_m.F:151`) calling `__kmpc_threadprivate_cached(gtid=145)` past `libiomp5`'s `__kmp_default_tp_capacity = 128`.
     2. `libiomp5.so` `__kmp_reap_worker` reading `0x18(%rdi)` at corrupted pointer `%rdi = 0x56` when short-lived OS caller threads tear down their OpenMP teams on exit.
   - Control experiment: single persistent worker (`mumps-worker`): **20/20 pass** at `OMP_NUM_THREADS=16` and `32`.

---

## 4. Root Cause and Fault Model (`file:line` and OpenMP Runtime Behavior)

### 4.1 Fault 1a — `dmumps_f77_` Pointer-Table Reallocation Race (`0x274520`)

In `src/dmumps_f77.F:102-107` (and identically in `src/[scz]mumps_f77.F`), `DMUMPS_F77` stores all active `DMUMPS_STRUC` Fortran pointers in three `SAVE` procedure-global variables:
```fortran
      TYPE (DMUMPS_STRUC_PTR), DIMENSION (:), POINTER, SAVE ::
     &  mumps_par_array
      INTEGER, SAVE :: DMUMPS_STRUC_ARRAY_SIZE = 0
      INTEGER, SAVE :: N_INSTANCES = 0
```
- On `JOB = -1` (`src/dmumps_f77.F:136-163`), if all slots are full (`line 145`), `DMUMPS_F77` allocates a larger `mumps_par_array_bis`, copies pointers, executes `DEALLOCATE(mumps_par_array)` (`line 149`), and rebinds `mumps_par_array => mumps_par_array_bis` (`line 150`).
- On `JOB = -2` (`src/dmumps_f77.F:449-457`), once `N_INSTANCES` reaches `0`, it executes `DEALLOCATE(mumps_par_array)` (`line 455`).
- Disassembly of `/usr/local/lib/mumps/libdmumps_cpmech.so` at `dmumps_f77_+0x89f` (`0x27450f..0x274528`, corresponding to `src/dmumps_f77.F:182` `IF (.NOT. associated(mumps_par_array(INSTANCE_NUMBER)%PTR))`):
  ```text
  27450f: 48 8b 05 ea f1 0a 00  mov  0xaf1ea(%rip),%rax   # dmumps_f77_$MUMPS_PAR_ARRAY (0x323700)
  274516: 49 63 ce              movslq %r14d,%rcx         # INSTANCE_NUMBER
  274519: 48 2b 0d 20 f2 0a 00  sub  0xaf220(%rip),%rcx   # lower bound
  274520: 48 0f af 0d 10 f2 0a  imul 0xaf210(%rip),%rcx   # stride
  274528: 48 8b 1c 08           mov  (%rax,%rcx,1),%rbx   # FAULT: %rax == 0 when deallocated concurrently
  ```
  When Thread B executes `DEALLOCATE(mumps_par_array)` (`line 149` or `line 455`), the descriptor base at `0x323700` becomes `0x0`, and Thread A's concurrent read at `0x274528` (attributed to bundle `0x274520`) segfaults at address `0x0`.

### 4.2 Fault 1b — C/Fortran Interface Static Pointer Laundering

- `src/mumps_c.c:244-249` declares 6 file-scope `static` pointers per arithmetic (`[SDCZ]MUMPS_COLSCA_STATIC`, `ROWSCA_STATIC`, `COLSCA_LOC_STATIC`, `ROWSCA_LOC_STATIC`, `RHSINTR_STATIC`, `SINGULAR_VALUES_STATIC`).
- `src/mumps_common.c:20-26` declares 6 process-global `static` pointers shared across **all four arithmetics** (`MUMPS_MAPPING`, `MUMPS_PIVNUL_LIST`, `MUMPS_SYM_PERM`, `MUMPS_UNS_PERM`, `MUMPS_GLOB2LOC_RHS`, `MUMPS_GLOB2LOC_SOL`).
- On every `[sdcz]mumps_c` call, `src/mumps_c.c:432-436` zeroes 5 of the 6 arithmetic statics (`RHSINTR_STATIC` is omitted), calls `MUMPS_F77` (`lines 600-639`), which sets all 12 static pointers via `MUMPS_ASSIGN_*` / `MUMPS_NULLIFY_C_*` (`src/dmumps_f77.F:380-440`), and then copies them into `mumps_par` (`src/mumps_c.c:651-688`).
- Concurrent calls—even between `dmumps_c` and `zmumps_c`—overwrite `MUMPS_SYM_PERM`, `MUMPS_COLSCA_STATIC`, etc., attaching Instance B's internal arrays to Instance A's struct.

### 4.3 Fault 1c — Module-Global State Across `JOB = 1, 2, 3`

- During `JOB = 1`, `MUMPS_DISTRIBUTE` (`src/mumps_static_mapping.F:93`, called from `src/dana_aux.F:4070`) and `MUMPS_INIT_ARCH_PARAMETERS` (`src/mumps_static_mapping.F:4317`, called from `src/dana_driver.F:572`) allocate, mutate, and deallocate 65 module globals (`src/mumps_static_mapping.F:24-91`), causing `SIGSEGV` or `Attempt to DEALLOCATE unallocated 'mem_distribtmp'` at `src/mumps_static_mapping.F:4392`.
- During `JOB = 2`, `MUMPS_LOAD_INIT` (`src/dfac_driver.F:2175` $\rightarrow$ `src/mumps_load.F:717-815`) and `MUMPS_LOAD_END` (`src/dfac_driver.F:5032` $\rightarrow$ `src/mumps_load.F:1142-1190`) allocate/deallocate `FUTURE_NIV2` (`src/fac_future_niv2_mod.F:15`) and 75 module variables in `MUMPS_LOAD` (`src/mumps_load.F:34-126`). Concurrent `JOB = 2` calls cause `PB allocation in MUMPS_LOAD_INIT` (`INFOG(1) = -13` at `src/mumps_load.F:726`), `MPI_ABORT` (`line 1062`), or out-of-bounds reads in `MUMPS_UPPER_PREDICT` (`line 4948`).
- During `JOB = 2` and `JOB = 3`, `DMUMPS_TMP_PTR` (`src/dstatic_ptr_m.F:16`) is written and read via `DMUMPS_DM_SET_PTR` (`src/dfac_mem_dynamic.F:483-486`, guarded by `CRITICAL(STATIC_PTR_ACCESS)`), **unprotected** in `src/dsol_fwd.F:156-157`, **unprotected** in `src/dsol_bwd.F:166-167`, and **unprotected** for Schur setup at `src/dmumps_driver.F:1662-1664, 2022-2025`.

### 4.4 Fault 2a — `libiomp5` `THREADPRIVATE` Cache Overflow (`DMUMPS_SOL_L0OMP_R`)

- `src/mumps_l0_omp_m.F:16-20` holds the **only** `!$OMP THREADPRIVATE` directive in the entire MUMPS 5.9.2 codebase:
  ```fortran
  LOGICAL, DIMENSION(:), POINTER :: NB_CORE_PER_THREAD_CHANGED
  INTEGER, DIMENSION(:), POINTER :: NB_CORE_PER_THREAD
  INTEGER :: THREAD_ID = 1
  LOGICAL :: IS_ROOT_OF_L0_OMP
  !$OMP   THREADPRIVATE ( THREAD_ID , IS_ROOT_OF_L0_OMP )
  ```
  `rg` across `src/` confirms `NB_CORE_PER_THREAD_CHANGED`, `NB_CORE_PER_THREAD`, and `IS_ROOT_OF_L0_OMP` are **completely unreferenced** dead code; only `THREAD_ID` is used, solely in `src/[sdcz]sol_omp_m.F:36,110-111,151,328,409-410,450` and `src/[sdcz]sol_lr.F:19,132,883,1991,2342`.
- Disassembly of `/usr/local/lib/mumps/libdmumps_cpmech.so` (`dmumps_sol_l0omp_m_mp_dmumps_sol_l0omp_r_` parallel outlined functions at `0x19e240..0x19e30f`):
  ```text
  # First !$OMP PARALLEL (src/dsol_omp_m.F:109-113)
  19e257: lea  0x1876e2(%rip),%r8    # &__tpv_ptr_mumps_l0_omp_m_mp_thread_id_ (0x325940)
  19e263: call __kmpc_threadprivate_cached@plt
  19e268: mov  %rax,%r14
  19e26b: movl $0x1,(%rax)           # THREAD_ID = 1
  19e271: call omp_get_thread_num@plt
  19e276: inc  %eax
  19e278: mov  %eax,(%r14)           # THREAD_ID = omp_get_thread_num() + 1

  # Second !$OMP PARALLEL (src/dsol_omp_m.F:139-151)
  19e2eb: lea  0x18764e(%rip),%r8    # &__tpv_ptr_mumps_l0_omp_m_mp_thread_id_
  19e2f7: call __kmpc_threadprivate_cached@plt
  19e2fc: mov  %rax,0x48(%rbx)       # later read at THREAD_ID_P = THREAD_ID (line 151) -> FAULT if %rax == 0
  ```
- Disassembly of `libiomp5.so` `<__kmpc_threadprivate_cached>` (`0xf4fc0..0xf5050`):
  - Every distinct OS thread entering an OpenMP parallel region becomes a separate OpenMP root thread owning a hot team of `OMP_NUM_THREADS` worker threads, consuming `OMP_NUM_THREADS` global thread IDs (`gtid`).
  - On a 32-CPU workstation, `libiomp5` initializes `__kmp_default_tp_capacity = max(128, 4 * 32) = 128` (`__kmp_tp_capacity`).
  - On the first call to `__kmpc_threadprivate_cached`, `libiomp5` allocates the 128-slot cache table at `__tpv_ptr_mumps_l0_omp_m_mp_thread_id_` (`0xf4fe6..0xf503c`) and sets `__kmp_tp_cached = 1` (`0xf5036`). Once `__kmp_tp_cached == 1`, `libiomp5` does **not** reallocate existing cache tables when `__kmp_threads` expands beyond 128.
  - Thus, once cumulative `N_caller_threads * OMP_NUM_THREADS > 128` (e.g., 10 `cargo test` threads at `OMP=13`, assigning `gtid = 145`), `__kmpc_threadprivate_cached` indexes past the 128-entry cache or returns `NULL`, segfaulting at `src/dsol_omp_m.F:151`.

### 4.5 Fault 2b — Caller-Thread Hot-Team Teardown Race (`__kmp_reap_worker`)

When transient OS caller threads (e.g. `cargo test` worker threads or short-lived `std::thread::spawn` tasks) exit after calling MUMPS, `libiomp5`'s POSIX TLS destructor runs `__kmp_internal_end_thread` $\rightarrow$ `__kmp_reap_worker` to join and destroy that caller thread's `OMP_NUM_THREADS` hot team (`0x18(%rdi)` read on thread descriptor). Under concurrent/sequential churn from many short-lived caller threads at `OMP_NUM_THREADS=16`, `__kmp_reap_worker` faults even when `KMP_ALL_THREADPRIVATE=512` expands the threadprivate cache (`0/10` pass at `OMP=16`). **Consequence:** Removing `!$OMP THREADPRIVATE` alone is necessary but **not sufficient**; MUMPS must execute OpenMP regions only from persistent MUMPS-owned OpenMP master thread(s) that never exit mid-process.

### 4.6 Probe Findings (§7 Verification Summary)

| Probe | Result & Verified Evidence |
| :--- | :--- |
| **Probe 1 (`dorival_thread_safe.F`)** | Confirmed against stock MUMPS: `locked=0` fails 6/6 (`SIGSEGV` / `INFOG(1)=-13` / `MPI_ABORT`), while `locked=1` and `locked=2` pass 6/6 (`ICNTL(16)=1`, fixed 8 OpenMP outer threads, zero `gtid > 128` overflow). |
| **Probe 2 (`ICNTL(48)=0`, `ICNTL(16)=1`, `ICNTL(16)=k`)** | `ICNTL(48)=0` sets `KEEP(400)=0` (`src/dana_driver.F:498-501`) and disables `DO_L0OMP_FWD`/`BWD` (`src/dsol_c.F:253-258`), bypassing `DMUMPS_SOL_L0OMP_R` and `__kmpc_threadprivate_cached` in non-BLR solves; however node-level OpenMP (`src/dfac_front_LU_type1.F:672,745,1047`) and MKL still allocate teams per OS caller thread, leaving `__kmp_reap_worker` exposed. `ICNTL(16)=k` shifts the overflow threshold to `N_callers * k > 128`. Note: `src/dmumps_driver.F:770, 2348-2350` already restores `PREVIOUS_OMP_THREADS_NUM` on exit (`GOTO 499/500` fall through to `line 2346`). |
| **Probe 3 (`128` threshold driver)** | Driven strictly by `N_OS_threads_entering_OMP * OMP_NUM_THREADS > 128` (`__kmp_tp_capacity`), **not** instance count. An OS thread holding multiple `DMUMPS_STRUC` instances uses one team (`OMP_NUM_THREADS` `gtid`s); an OS thread that never enters MUMPS consumes 0 `gtid`s. |
| **Probe 4 (`libmpiseq` & single-instance reentrancy)** | `nm -A /usr/local/lib/mumps/libmpiseq_cpmech.a` shows only `mpif_libseq_` (`libseq/mpif.h:108`), which is strictly a read-only sentinel address checked via `MUMPS_CHECKADDREQUAL` (`libseq/mpi.f:795`, `libseq/mpic.c:40-50`). `libmpiseq` is 100% thread-safe. However, `dmumps_c` is **not** reentrant across threads even with 1 instance per thread due to shared Fortran/C modules. Crucially, for in-core (`ICNTL(22)=0`), MUMPS **is** phase-reentrant on a single thread (`A JOB=1` $\rightarrow$ `B JOB=1` $\rightarrow$ `A JOB=2` $\rightarrow$ `B JOB=2` $\rightarrow$ `A JOB=3` $\rightarrow$ `B JOB=3`) because inter-phase handles are encoded into `id%FDM_F_ENCODING` (`src/front_data_mgt_m.F:143-233`) and `id%BLRARRAY_ENCODING` (`src/dmumps_lr_data_m.F:328-406`). |
| **Probe 5 (`mumps_c.c` statics & `DMUMPS_TMP_PTR` in assembled runs)** | Even with `ICNTL(6)=0, ICNTL(8)=0`, `src/dmumps_driver.F:1926-1934` allocates dummy length-1 `id%COLSCA(1)` and `id%ROWSCA(1)`, so `DMUMPS_COLSCA_STATIC`, `DMUMPS_ROWSCA_STATIC`, and `MUMPS_SYM_PERM` (`src/dmumps_f77.F:396-425`) are actively written on every assembled run. `DMUMPS_TMP_PTR` is exercised on every factorization (`src/dfac_mem_dynamic.F:483-486`) and **unprotected** on every solve (`src/dsol_fwd.F:156-157`, `src/dsol_bwd.F:166-167`). |
| **Probe 6 (OpenMP lock scopes)** | `LOCK_FOR_L0OMP` (`src/[sdcz]mumps_ooc.F:19`) is per-arithmetic module-global (`OMP_INIT_LOCK` at `line 151`, `OMP_SET_LOCK`/`OMP_TEST_LOCK`/`OMP_UNSET_LOCK` at `lines 3018, 3022, 3076`, never destroyed). `LOCK_FOR_SCATTER` (`src/sol_omp_common_m.F:24`) is process-global across all arithmetics (`MUMPS_SOL_L0OMP_LI` at `lines 50-55` called at `src/dsol_driver.F:3198`; `OMP_SET_LOCK`/`UNSET_LOCK` at `src/dsol_fwd_aux.F:900,915`; `MUMPS_SOL_L0OMP_LD` at `lines 78-85` called at `src/dsol_driver.F:6468`). |

---

## 5. Complete State-Inventory Table

Audit of `nm -A` across `/usr/local/lib/mumps/lib{s,d,c,z}mumps_cpmech.a`, `libmumps_common_cpmech.a`, `libpord_cpmech.a`, and `libmpiseq_cpmech.a` plus source inspection reveals 23 distinct state groups:

| # | Module / File | Symbols (Verified via `nm` & Source) | Exact `file:line` | Classification | Callers (`JOB` Phase) | Failure Mode if Shared & Action |
| :- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | `MUMPS_L0_OMP_M` | `THREAD_ID` (`!$OMP THREADPRIVATE`), dead `IS_ROOT_OF_L0_OMP`, `NB_CORE_PER_THREAD*` | `src/mumps_l0_omp_m.F:16-20` | `THREADPRIVATE` + dead globals | `[sdcz]sol_omp_m.F:36,109-151,328,407-450`; `[sdcz]sol_lr.F:19,132,883,1991,2342` (`JOB=3`) | Overflows `libiomp5` 128-entry TP cache (`gtid >= 128`). **Delete directive**; use `OMP_GET_THREAD_NUM()+1` and `OMP_GET_ANCESTOR_THREAD_NUM(1)`. |
| 2 | `MUMPS_SOL_L0OMP_M` | `LOCK_FOR_SCATTER(:)` (`OMP_LOCK_KIND`) | `src/sol_omp_common_m.F:24` | `process-global-mutable` | `src/[sdcz]sol_driver.F:3198,6468`; `src/[sdcz]sol_fwd_aux.F:900,915` (`JOB=3`) | Double `ALLOCATE` abort or lock destroyed mid-solve. Move array into `idintr` (`src/[sdcz]mumps_intr_types.F:56-65`). |
| 3 | `MUMPS_TPS_M` / `[SDCZ]MUMPS_TPS_M` | Derived types `MUMPS_TPS_T`, `[SDCZ]MUMPS_TPS_T` only | `src/omp_tps_common_m.F:14-28`; `src/domp_tps_m.F` | `per-instance` (type defs only) | `src/dfac_omp_m.F:36,45-46` (`JOB=2`) | Instances passed on call stack (`MUMPS_TPS_ARR`). **Safe as-is.** |
| 4 | `[SDCZ]MUMPS_F77` | `mumps_par_array(:)`, `[SDCZ]MUMPS_STRUC_ARRAY_SIZE`, `N_INSTANCES` | `src/[sdcz]mumps_f77.F:102-107` | `process-global-mutable` | Every `[sdcz]mumps_c` call (`JOB=-1..6,-2`) | `DEALLOCATE(mumps_par_array)` (`lines 149, 455`) causes `SIGSEGV` at `0x274520` (`line 182`). Protect slot table with mutex and never free live slot array. |
| 5 | `mumps_c.c` (per arith.) | `[SDCZ]MUMPS_COLSCA_STATIC`, `ROWSCA_STATIC`, `COLSCA_LOC_STATIC`, `ROWSCA_LOC_STATIC`, `RHSINTR_STATIC`, `SINGULAR_VALUES_STATIC` | `src/mumps_c.c:244-249` | `process-global-mutable` | Every `[sdcz]mumps_c` call (`lines 432-436, 657-688`) | Cross-instance pointer aliasing (`colsca`, `rowsca`, etc.). Pass output pointer struct directly via `MUMPS_F77`. |
| 6 | `mumps_common.c` (all arith.) | `MUMPS_MAPPING`, `MUMPS_PIVNUL_LIST`, `MUMPS_SYM_PERM`, `MUMPS_UNS_PERM`, `MUMPS_GLOB2LOC_RHS`, `MUMPS_GLOB2LOC_SOL` | `src/mumps_common.c:20-26` | `process-global-mutable` | Every `[sdcz]mumps_c` call via `src/[sdcz]mumps_f77.F:380-440` and `src/mumps_c.c:651-661` | Cross-instance and cross-arithmetic (`d` vs `z`) pointer overwrite. Pass via `MUMPS_F77` output struct. |
| 7 | `[SDCZ]MUMPS_STATIC_PTR_M` | `[SDCZ]MUMPS_TMP_PTR(:)` | `src/[sdcz]static_ptr_m.F:16` | `process-global-mutable` | `src/[sdcz]fac_mem_dynamic.F:483-486`; `src/[sdcz]sol_fwd.F:156-157`; `src/[sdcz]sol_bwd.F:166-167`; `src/[sdcz]mumps_driver.F:1662-1664, 2022-2025` | Unprotected in solve (`sol_fwd`/`sol_bwd`); overwrites `A_PTR` across instances. Mark dummy args `TARGET` (`A_PTR => A`) or use `C_F_POINTER`. |
| 8 | `MUMPS_STATIC_MAPPING` | 65 globals (`cv_*`, `mem_distribtmp`, `table_of_process`, `mem_distribmpi`, `allowed_nodes`, `score`) + dead `init1, init2` | `src/mumps_static_mapping.F:24-91, 3957-3958` | `process-global-mutable` | `MUMPS_INIT_ARCH_PARAMETERS` (`src/dana_driver.F:572`), `MUMPS_DISTRIBUTE` (`src/dana_aux.F:4070`), `MUMPS_END_ARCH_CV` (`src/dana_driver.F:2462`) (`JOB=1`) | `Attempt to DEALLOCATE unallocated 'mem_distribtmp'` at `line 4392` or `SIGSEGV`. Move into derived type passed on stack/`idintr`. |
| 9 | `MUMPS_LOAD` & `MUMPS_FUTURE_NIV2` | 75 `SAVE` variables/arrays/pointers (`LOAD_FLOPS`, `BUF_LOAD_RECV`, `WLOAD`, `NB_SON`, `DM_MEM`, etc.) + `FUTURE_NIV2(:)` | `src/mumps_load.F:34-126`; `src/fac_future_niv2_mod.F:15` | `process-global-mutable` | `MUMPS_LOAD_INIT` (`src/dfac_driver.F:2175`), `MUMPS_LOAD_END` (`src/dfac_driver.F:5032`) (`JOB=2`) | `PB allocation in MUMPS_LOAD_INIT` (`line 726`), `MPI_ABORT` (`line 1062`), `SIGSEGV` in `MUMPS_UPPER_PREDICT` (`line 4948`). Move into `idintr`. |
| 10 | `MUMPS_BUF_COMMON` & `[SDCZ]MUMPS_BUF` | `SIZEofINT`, `SIZEofREAL`, `BUF_MYID`, `SIZE_RBUF_BYTES`, `BUF_CB`, `BUF_SMALL`, `BUF_LOAD`, `BUF_LMAX_ARRAY`, `BUF_MAX_ARRAY(:)` | `src/mumps_comm_buffer_common.F:38-46`; `src/[sdcz]mumps_comm_buffer.F:19-21` | `process-global-mutable` | `MUMPS_BUF_INI_MYID` (`src/dfac_driver.F:2163`), all `JOB=2` comm/load calls | Concurrent buffer deallocation/reallocation crash. Move into `idintr`. |
| 11 | `MUMPS_FRONT_DATA_MGT_M` | `FDM_A`, `FDM_F` (encoded into `id%FDM_F_ENCODING` between `JOB=2` and `JOB=3`) | `src/front_data_mgt_m.F:70, 143-233` | `process-global-mutable` (during active `JOB=2,3`) | `src/dfac_driver.F:4942`, `src/dsol_driver.F:1961` (`JOB=2,3`) | Overlapping `JOB=2/3` corrupts `FDM_F` stack/counts. Store `FDM_A`/`FDM_F` on `idintr` / active worker state. |
| 12 | `MUMPS_FAC_DESCBAND_DATA_M` & `MUMPS_FAC_MAPROW_DATA_M` | `INODE_WAITED_FOR`, `FDBD_ARRAY`, `FMRD_ARRAY` | `src/fac_descband_data_m.F:17,28`; `src/fac_maprow_data_m.F:45` | `process-global-mutable` | `JOB=2` factorization | Overlapping `JOB=2` overwrites front descriptor tables. Move into `idintr`. |
| 13 | `[SDCZ]MUMPS_LR_DATA_M`, `MUMPS_LR_COMMON`, `MUMPS_LR_STATS` | `BLR_ARRAY` (encoded to `id%BLRARRAY_ENCODING` at `lines 328-406`), 17 precision `SAVE` vars, 67 stats vars | `src/[sdcz]mumps_lr_data_m.F:208-209`; `src/lr_common.F:20-42`; `src/lr_stats.F:30-128` | `process-global-mutable` | `src/dfac_driver.F:4903-4911, 5314-5315` (`JOB=1,2,3` when BLR active; `MUMPS_LR_COMMON` written on every `JOB=2,3`) | Concurrent BLR corrupts mixed-precision tables and `BLR_ARRAY`. Move into `idintr` (or serialize BLR in Phase 2). |
| 14 | `mumps_sol_omp_memory_manager.c` | `static uint8_t **arena`, `nb_slots`, `surface` | `src/mumps_sol_omp_memory_manager.c:32-34` | `process-global-mutable` | `src/dsol_driver.F:1737, 6785`; `src/dsol_lr.F:146` (`JOB=3` BLR) | Double-alloc returns error 1 (`line 42`) or frees live workspace. Make `_Thread_local` or store handle in `idintr`. |
| 15 | `[sdcz]sol_aux.F` | `[SDCZ]MUMPS_SOL_B`: `SAVE ITER, J, JLAST, JUMP`; `OMEGA`: `SAVE OM1, OLDOMG`; `LCOND`: `SAVE LCOND1, LCOND2, JUMP, DXIMAX, DXMAX` | `src/[sdcz]sol_aux.F:319, 883, 958` | `process-global-mutable` | `src/dsol_driver.F:5846, 6080, 6117` (`JOB=3` with `ICNTL(10)!=0` or `ICNTL(11)>0`) | Concurrent iterative refinement / error analysis corrupts `JUMP` state machine. Pass state derived type from stack. |
| 16 | `[sdcz]fac_process_message.F` | `INTEGER, SAVE :: RECURS = 0` in `[SDCZ]MUMPS_TRY_RECVTREAT` | `src/[sdcz]fac_process_message.F:650` | `process-global-mutable` | `JOB=2` message processing | Concurrent `RECURS` corruption disables early exit at `line 795`. Move to instance/worker state. |
| 17 | `[SDCZ]MUMPS_CONFIG_FILE` | `PARAM=' '`, `INDICE=' '`, `VAL=' '` (implicit `SAVE` from declaration initialization) | `src/[sdcz]mumps_config_file.F:26-28` | `process-global-mutable` | `src/dmumps_driver.F:608` (every `JOB != -1,-4,-3,7,8`) | Data race if config file enabled. Move `=' '` init to executable body. |
| 18 | `MUMPS_MEMORY_MOD` | `ISIZE, I8SIZE, SSIZE, DSIZE, CSIZE, ZSIZE` | `src/mumps_memory_mod.F:23` | `process-global-read-only` (after idempotent init) | `MUMPS_MEMORY_SET_DATA_SIZES` (`src/dmumps_driver.F:1157`) | Benign race setting constants `4,8,4,8,8,16`. Replace with `PARAMETER` constants or `c_sizeof`. |
| 19 | `MUMPS_PIVNUL_MOD` | Derived type `PIVNUL_LIST_STRUCT_T`; `CRITICAL(critical_pivnul)` | `src/mumps_pivnul_mod.F:15-18, 34` | `per-instance` | `src/dfac_driver.F:185, 1988-1994, 4674-4676, 5070-5072` | Instance lives on `DFAC_DRIVER` stack (`id%PIVNUL_LIST` output). **Safe.** |
| 20 | `[SDCZ]MUMPS_PARALLEL_ANALYSIS` | `MP, MPG, LP, ...` and `DMUMPS_SEND_BUF` `SAVE`: `INIT, SPACE, PENDING, REQ, CPNT` | `src/[sdcz]ana_aux_par.F:51-54, 3103-3106` | `process-global-mutable` | Only when `ICNTL(28)=2` & `NPROCS > 1` (`src/dana_driver.F:847-865`) | Unused under `libmpiseq` (`NPROCS=1`). Out of scope for concurrent multi-instance. |
| 21 | Out-of-Core (`OOC`) subsystem | `MUMPS_OOC_COMMON`, `[SDCZ]MUMPS_OOC` (`LOCK_FOR_L0OMP`), `[SDCZ]MUMPS_OOC_BUFFER`, `[SDCZ]MUMPS_SOL_ES`, `mumps_io*.c` (`mumps_files`, `io_queue`, `err_mutex`, `MUMPS_OOC_STORE_*`) | `src/mumps_ooc_common.F:18-34`; `src/[sdcz]mumps_ooc.F:19-55`; `src/[sdcz]mumps_ooc_buffer.F:20-35`; `src/[sdcz]mumps_sol_es.F:27`; `src/mumps_io.h:31-34`; `src/mumps_io.c:24-26`; `src/mumps_io_basic.c:23-34`; `src/mumps_io_err.c:23-28`; `src/mumps_io_thread.c:22-38` | `process-global-mutable` (persists `JOB=2` $\rightarrow$ `JOB=3`) | `ICNTL(22)=1` and `DMUMPS_SOL_ES_INIT` (`KEEP(201)>0`) | Concurrent OOC destroys factor files and async I/O queue. Serialize via OOC mutex on Worker 0 in Phase 2. |
| 22 | `libmpiseq` | `COMMON/mpif_libseq/MPI_IN_PLACE` | `libseq/mpif.h:107-108`; `libseq/mpi.f:795`; `libseq/mpic.c:40-50` | `process-global-read-only` | Address comparison sentinel only | Confirmed via `nm`: zero mutable state. **Safe.** |
| 23 | Orderings (`PORD`, `METIS`, `AMD/AMF/QAMD`) & `GPU` | Zero globals in `libpord` & `libmetis.so`; `libGKlib.so` (`isaac_*` TLS + `gk_jbuf`); GPU `dmumps_handle` | `src/ana_orderings.F`; `src/ana_orderings_wrappers_m.F:57-344`; `src/mumps_metis.c:109-196`; `src/dmumps_gpu.c:32-45`; `src/mumps_register_thread.c:59-60` | PORD/METIS/AMD: **reentrant**; GPU: **uncompiled/out of scope** | `JOB=1` (`ICNTL(7)`), `ICNTL(51)` | PORD/AMD/METIS safe concurrently (`libGKlib` uses `.tbss` per-thread RNG); GPU uncompiled under `MakefileMKL.inc`. |

---

## 6. Design Options as Numbered Decisions

Owners can approve or adjust the design by answering Decisions 1 through 4 by number.

### Decision 1 — Core Concurrency Architecture (Options 1a–1f)

| Option | Design Summary | Exact Changes (`file:line`) | ABI Impact | Effort (Person-Weeks) | Performance Model | Risk & Falsification Experiment |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **1a (Control — Rejected)** | Process-wide `pthread_mutex_t` inside `[sdcz]mumps_c`. | `src/mumps_c.c:365` | None | **0.5 pw** | +0.1 µs/call mutex; 0 % concurrency; `N_callers * OMP` teams in `libiomp5`. | **Falsified today:** Fails Failure Class 2 (`0/6` at `OMP=16` in `russell_sparse`) due to `__kmpc_threadprivate_cached` overflow (`src/mumps_l0_omp_m.F:20`) and `__kmp_reap_worker`. |
| **1b** | MUMPS-owned single persistent worker thread dispatched inside `libmumps_common` + remove `THREADPRIVATE`. | New `src/mumps_thread_safe.c`; `src/mumps_c.c:244-249, 365`; `src/mumps_l0_omp_m.F:16-20`; `src/[sdcz]sol_omp_m.F:109-151, 407-450`; `src/[sdcz]sol_lr.F:132,883,1991,2342` | Additive capability symbol only | **1.5 pw** | **≈4–6 µs/call** C mutex/condvar hop (saves Rust wrapper layers vs `mumps-worker`); 1 OpenMP hot team (`<= OMP_NUM_THREADS` `gtid`s); 0 % concurrent solve overlap; memory same as 1 master + $N$ instances. | Low risk; reproduces `mumps-worker` `20/20` reliability inside MUMPS itself. Falsified only if 1-instance throughput is insufficient for multi-instance `pmsim` workloads. |
| **1c (Rejected)** | One worker OS thread spawned per `DMUMPS_STRUC` instance (`JOB=-1` $\rightarrow$ `pthread_create`, `JOB=-2` $\rightarrow$ `pthread_join`) + data reentrancy. | Same data changes as Phase 2 + per-instance queue in `src/mumps_c.c` | Same | **4.5 pw** | Unbounded threads: $M$ instances $\times$ `OMP_NUM_THREADS` (`16 * 16 = 256` on 32 cores $\rightarrow$ 8× oversubscription + `>128` `gtid`s); every `JOB=-2` triggers `__kmp_reap_worker`. | **Falsified** by churning 20 instances $\times$ `OMP=16`: triggers Failure Class 2 (`__kmp_reap_worker`) on `JOB=-2` teardown. |
| **1d (Rejected)** | 100% full reentrancy across all of MUMPS (including OOC, BLR, MPI-to-OMP, `mumps_io_thread.c`, `mumps_c.c`) executing directly on caller OS threads (no worker pool). | ~83 files, >250 globals (Table rows 1–21) | Same | **16.0–24.0 pw** | 0 µs dispatch overhead, but short-lived caller OS threads spawn & reap OpenMP hot teams repeatedly and oversubscribe CPUs under `Rayon`/`cargo test`. | **Falsified** on `libiomp5`: even with zero Fortran globals, running OpenMP on >128 churning OS caller threads triggers `__kmp_reap_worker` faults. Also fails cost/benefit due to OOC complexity unused by clients. |
| **1e (Recommended)** | **Hybrid Phased:** **Phase 0+1:** Remove `THREADPRIVATE`, `TMP_PTR`, & C/F77 glue statics + MUMPS-owned persistent worker pool ($K=1$). **Phase 2:** Move in-core module state (Table rows 2, 7–16) into `idintr` / TLS and scale persistent pool to **$K$ immortal OpenMP masters** (`instance_number % K`). | Phase 0+1: Option 1b files + `src/[sdcz]static_ptr_m.F`, `src/mumps_common.c:20-26`, `src/[sdcz]mumps_f77.F:102-182`. Phase 2: `src/[sdcz]mumps_intr_types.F:56-65` + Table rows 2, 7–16, OOC mutex. | Additive capability symbol only; `DMUMPS_STRUC_C` untouched | **1.5 pw (Phase 0+1)**; **+4.5 pw (Phase 2)**; **6.0 pw total** | **Phase 1:** ≈4–6 µs/call, 1 persistent OpenMP team. **Phase 2:** Bounded $K$ teams ($K \times \text{OMP} \le N_{\text{cores}}$), zero team teardown (`__kmp_reap_worker` never fires on churn), **up to $K\times$ concurrent solve speedup**, zero client locks. | Phase 1 de-risks immediate deployment in 1.5 pw; Phase 2 delivers true concurrency without `libiomp5` churn crashes. Falsification: Phase 2 must deliver $\ge 4\times$ speedup on 8 concurrent 3D instances (`OMP=4`) and 20/20 clean stress runs. |
| **1f (Rejected)** | Scope-restricted mode: force `ICNTL(16)=1`, `ICNTL(48)=0` + full data reentrancy on caller threads. | Force `ICNTL(16)=1`, `ICNTL(48)=0` in `src/dmumps_driver.F:766` + Table rows 1–16 | Violates AC-3 (`ICNTL(16)=0`, `ICNTL(48)=1`) | **4.5 pw** | Loses intra-instance OpenMP speedup (`4–16×` slower on single large 3D matrices $\sim 10^5$ dof). | **Falsified** by large 3D single-solve slowdown (violates §3 AC-3) and still triggers `__kmp_reap_worker` if >128 OS threads enter single-thread OpenMP regions and exit. |

**Recommendation for Decision 1:** Adopt **Option 1e** (implementing Phase 0 + Phase 1 first at **1.5 person-weeks**, then Phase 2 at **+4.5 person-weeks**). Fallback if `pmsim` rarely solves independent instances simultaneously: stop cleanly after **Phase 1** (1.5 pw), which already deletes `russell_sparse/src/mumps_lock.rs`.

### Decision 2 — ABI Capability Detection Symbol
- Expose `int mumps_is_thread_safe(void)` in `libmumps_common_cpmech.{a,so}` (returning `1` in Phase 1, `2` in Phase 2) and `#define MUMPS_THREAD_SAFE_VERSION` in `include/mumps_c_types.h`.
- **Recommendation:** Approve (zero change to `DMUMPS_STRUC_C`).

### Decision 3 — Fortran Direct Entry (`dorival_thread_safe.F`) vs. C-Entry Only (`[sdcz]mumps_c`)
- Option A: Intercept only `[sdcz]mumps_c` (protects all C/Rust clients; direct Fortran callers `CALL DMUMPS(ID)` run unlocked unless reentrant in Phase 2).
- Option B: Route both `[sdcz]mumps_c` AND Fortran `[SDCZ]MUMPS` (`src/[sdcz]mumps_driver.F` renamed internally to `[SDCZ]MUMPS_IMPL`, with `[SDCZ]MUMPS` calling the C worker dispatcher via `ISO_C_BINDING`) through the persistent pool.
- **Recommendation:** Approve **Option B** (+0.2 pw included in Phase 1) so `examples/dorival_thread_safe.F` (`locked=0`) passes out of the box alongside C/Rust callers.

### Decision 4 — Phase 2 Concurrency Pool Width ($K$)
- Set default persistent worker count $K = \max(1, \min(8, N_{\text{cpus}} / \text{OMP\_NUM\_THREADS}))$ overridable via `MUMPS_NUM_WORKERS`, and automatically clamp per-worker `omp_set_num_threads` or keep user `OMP_NUM_THREADS` if `MUMPS_NUM_WORKERS` is explicitly tuned so total `libiomp5` `gtid`s stay far below 128 (`<= 32..64`). Note that calls already running on a Rayon thread whose `OMP=1` (e.g. `pmsim` tiny 192-dof tests) could alternatively fast-path when $K=N$ reentrant, except `libiomp5` teardown safety still favors the persistent $K$-worker pool.
- **Recommendation:** Approve default $K = 4$ (when `OMP_NUM_THREADS <= 8`) / auto-scaled, env-overridable via `MUMPS_NUM_WORKERS`.

---

## 7. Recommended Design in Phases (Gates, Stop Conditions, and Proving Tests)

### Phase 0 — Eliminate `!$OMP THREADPRIVATE`, `DMUMPS_TMP_PTR`, and C/F77 Interface Statics (**0.5 person-weeks**)

1. **Delete `!$OMP THREADPRIVATE` and dead globals (`src/mumps_l0_omp_m.F:14-21`):**
   - Remove `NB_CORE_PER_THREAD_CHANGED`, `NB_CORE_PER_THREAD`, `THREAD_ID`, `IS_ROOT_OF_L0_OMP`, and line 20's `!$OMP THREADPRIVATE ( THREAD_ID , IS_ROOT_OF_L0_OMP )`.
   - In `src/[sdcz]sol_omp_m.F:36, 109-113, 139-151, 328, 407-412, 439-450`, remove `USE MUMPS_L0_OMP_M, ONLY: THREAD_ID` and initialize local `PRIVATE(THREAD_ID_P)` directly via `THREAD_ID_P = OMP_GET_THREAD_NUM() + 1` inside the `!$OMP PARALLEL` region at lines `139` and `439` (the exact pattern already used in `src/dfac_omp_m.F:255-267`!). Delete the now-redundant dummy `!$OMP PARALLEL` blocks at `lines 109-113` and `lines 407-412`.
   - In `src/[sdcz]sol_lr.F:19, 132, 883, 1991, 2342`, replace `thread_id - 1` with `OMP_GET_ANCESTOR_THREAD_NUM(1)` (checking `OMP_GET_ACTIVE_LEVEL() > 1` or passing `THREAD_ID_P - 1` down from `sol_omp_m.F`, zero TLS required).
2. **Eliminate C/Fortran static pointer laundering (`src/mumps_c.c`, `src/mumps_common.c`, `src/[sdcz]mumps_f77.F`):**
   - Remove the 6 per-arithmetic statics in `src/mumps_c.c:244-249` and 6 shared statics in `src/mumps_common.c:20-26`. Pass a stack-allocated C output-pointer array `void *out_ptrs[12]` as an explicit argument in the internal C$\leftrightarrow$F77 call (`src/mumps_c.c:600-639` $\rightarrow$ `src/[sdcz]mumps_f77.F:380-440`).
   - Protect `mumps_par_array` (`src/[sdcz]mumps_f77.F:102-182, 449-457`) with a small C mutex exported from `libmumps_common` (`MUMPS_F77_LOCK`/`UNLOCK`), store each instance in an immutable bucket/slot or keep a C-side hash map/free-list of `DMUMPS_STRUC` containers, and **never** call `DEALLOCATE(mumps_par_array)` (`lines 149, 455`) while other instances exist (avoiding the `0x274520` null descriptor fault forever).
3. **Eliminate `DMUMPS_TMP_PTR` (`src/[sdcz]static_ptr_m.F:14-32`):**
   - In `src/[sdcz]sol_fwd.F:156-157` and `src/[sdcz]sol_bwd.F:166-167`, give dummy argument `A` the `TARGET` attribute and use `A_PTR => A`.
   - In `src/[sdcz]fac_mem_dynamic.F:483-486` and `src/[sdcz]mumps_driver.F:1662-1664, 2022-2025`, replace `DMUMPS_SET_TMP_PTR` / `DMUMPS_GET_TMP_PTR` with standard Fortran 2003 `ISO_C_BINDING` `C_LOC` / `C_F_POINTER`. Delete `CRITICAL(STATIC_PTR_ACCESS)` (`src/dfac_mem_dynamic.F:483,486`).
4. **Fix minor `SAVE` hazards:**
   - Separate variable declarations from `' '` assignment in `src/[sdcz]mumps_config_file.F:26-28`.
   - Make `ISIZE, I8SIZE, SSIZE, DSIZE, CSIZE, ZSIZE` `INTEGER(8), PARAMETER` in `src/mumps_memory_mod.F:23`.
- **Phase 0 Gate:** `nm -A lib{s,d,c,z}mumps_cpmech.a libmumps_common_cpmech.a | grep -E '__tpv_ptr|STATIC_PTR|_STATIC\b|MUMPS_SYM_PERM'` returns **0 symbols**.
- **Phase 0 Proving Test:** Existing MUMPS `examples/*simpletest` and `dsimpletest_save_restore` pass cleanly on `/tmp/mumps-probe`.

### Phase 1 — MUMPS-Owned Persistent OpenMP Master Dispatcher (**1.0 person-week; 1.5 person-weeks cumulative**)

1. Add `src/mumps_thread_safe.c` (added to `OBJS_C_COMMON` in `src/Makefile:98-120` and linked automatically into `libmumps_common_cpmech.{a,so}`).
2. Implement `pthread_once` lazy initialization of **1** persistent POSIX worker thread (`mumps-omp-master-0`) with a lock-free / `pthread_mutex_t` + `pthread_cond_t` work queue (`max(4, N)` zero-allocation stack-resident request struct per waiter).
3. Route `[sdcz]mumps_c` (`src/mumps_c.c:365`) and optional Fortran entry `[SDCZ]MUMPS` (`src/[sdcz]mumps_driver.F`) into `mumps_ts_dispatch(arith, struc_ptr)` (executing inline without queueing if `pthread_self()` is already a MUMPS worker thread, preventing nested deadlocks).
4. Expose `int mumps_is_thread_safe(void)` returning `1` (see Section 8).
- **Phase 1 Gate:** Zero kernel/OpenMP faults and 100% pass rate on all churn tests (Section 9) with `ICNTL(16)=0, ICNTL(48)=1` enabled.
- **Phase 1 Stop Condition:** **Safe production milestone.** At 1.5 person-weeks total, `russell_sparse` can delete `mumps-worker` (`mumps_lock.rs`) immediately. Proceed to Phase 2 only when concurrent multi-instance solve speedup is required (Decision 1).
- **Phase 1 Proving Test:** `examples/dorival_thread_safe 16 25 200 0` (`locked=0`), C pthread churn (`32 threads × 20 cycles × OMP=16`), and 20/20 `russell_sparse` + `pmsim` test runs at `OMP_NUM_THREADS=16, 32`.

### Phase 2 — In-Core Module Reentrancy + Bounded $K$-Worker Concurrent Pool (**4.5 person-weeks; 6.0 person-weeks cumulative**)

1. **Encapsulate per-phase module state into `idintr` (`src/[sdcz]mumps_intr_types.F:56-65`) and stack structs:**
   - **Mapping (`JOB=1`):** Pack the 65 globals in `src/mumps_static_mapping.F:24-91` into a derived type `MUMPS_STATIC_MAPPING_T` allocated cleanly on the stack of `MUMPS_DISTRIBUTE` / `MUMPS_INIT_ARCH_PARAMETERS` (`src/mumps_static_mapping.F:93,4317`) and passed explicitly; delete `SAVE init1, init2` (`lines 3957-3958`).
   - **Load & comm buffers (`JOB=2,3`):** Encapsulate `MUMPS_LOAD` (`src/mumps_load.F:34-126`), `MUMPS_FUTURE_NIV2` (`src/fac_future_niv2_mod.F:15`), `MUMPS_BUF_COMMON` (`src/mumps_comm_buffer_common.F:38-46`), and `[SDCZ]MUMPS_BUF` (`src/[sdcz]mumps_comm_buffer.F:19-21`) into `idintr` (or a per-worker master state slot passed down to `[SDCZ]MUMPS_FAC_DRIVER`, never using `!$OMP THREADPRIVATE`).
   - **Front/BLR/Solve modules:** Move `FDM_A, FDM_F` (`src/front_data_mgt_m.F:70`), `FDBD_ARRAY`, `INODE_WAITED_FOR` (`src/fac_descband_data_m.F:17,28`), `FMRD_ARRAY` (`src/fac_maprow_data_m.F:45`), `LOCK_FOR_SCATTER` (`src/sol_omp_common_m.F:24`), `RECURS` (`src/[sdcz]fac_process_message.F:650`), and BLR globals (`src/[sdcz]mumps_lr_data_m.F:208-209`, `src/lr_common.F:20-42`, `src/lr_stats.F:30-128`) into `idintr` / per-worker master state; make `arena, nb_slots, surface` in `src/mumps_sol_omp_memory_manager.c:32-34` `_Thread_local`; refactor `[sdcz]sol_aux.F:319,883,958` `SAVE` variables into explicit `INOUT` state structs held in `[SDCZ]MUMPS_SOL_DRIVER`.
   - **Guard OOC / Save-Restore:** Add a mutex/routing check in `src/mumps_thread_safe.c` so any call with `ICNTL(22) > 0` (`OOC`) or `JOB = 7, 8, -3` runs serialized on Worker 0 while in-core calls run freely in parallel across Workers `0..K-1`.
2. **Widen worker pool (`src/mumps_thread_safe.c`) to $K$ persistent OpenMP masters** (`MUMPS_NUM_WORKERS`, default $K=4$ capped so $K \times \text{OMP} \le 64$) and bump `mumps_is_thread_safe()` return value to `2`. Note: Since inter-phase state stays inside `id` / `idintr`, any idle worker can service any `JOB`, or instances pin to `instance_number % K` for L2/NUMA cache locality.
- **Phase 2 Gate:** Zero global mutable symbols remain in non-OOC/non-MPI-par-ana `.o` files (`nm` check) and ThreadSanitizer / Helgrind reports 0 races on concurrent in-core assembled solves.
- **Phase 2 Proving Test:** 8 caller threads concurrently factoring & solving 8 independent 2,187-dof 3D FEM matrices (`OMP_NUM_THREADS=4`, `MUMPS_NUM_WORKERS=8`) achieve **$\ge 4.0\times$ wall-clock speedup** over Phase 1 ($K=1$) and bit-identical answers at `OMP=1`.

---

## 8. ABI Preservation and Build-Time / Run-Time Capability Detection

1. **Zero Breaking Changes:** `DMUMPS_STRUC_C` and `[SDCZ]MUMPS_STRUC_C` (`include/[sdcz]mumps_c.h`), `MUMPS_VERSION` (`"5.9.2"`), and `[sdcz]mumps_c` function signatures remain **100% unchanged**. Existing client binaries link seamlessly against `lib*mumps_cpmech.so`.
2. **New Capability Symbols (in `include/mumps_c_types.h` and `libmumps_common_cpmech.{a,so}`):**
   ```c
   /* Thread-safety capability query (absent in upstream stock MUMPS 5.9.2):
    *   0 = not thread-safe (never returned by this symbol)
    *   1 = Phase 1 (persistent worker-serialized inside MUMPS; client locks unnecessary)
    *   2 = Phase 2 (concurrent K-worker execution for in-core assembled matrices)
    */
   #define MUMPS_THREAD_SAFE_VERSION 1
   int MUMPS_CALL mumps_is_thread_safe(void);
   ```
3. **What `russell_sparse` Must Do to Drop `mumps-worker`:**
   - In `russell_sparse/build.rs`, perform a compile-and-link probe (or `nm -g /usr/local/lib/mumps/libmumps_common_cpmech.a | grep mumps_is_thread_safe`): if `mumps_is_thread_safe` resolves, emit `cargo:rustc-cfg=mumps_thread_safe`.
   - Under `#[cfg(mumps_thread_safe)]`, `russell_sparse` calls `solver_mumps_*` directly on the caller's thread (no `mumps-worker`, no channel, no `Mutex<()>`).
   - Under `#[cfg(not(mumps_thread_safe))]` (e.g., stock Debian `libmumps-seq-dev` or stock upstream tarball in CI / crates.io users), `russell_sparse/src/mumps_lock.rs` remains compiled and active as the safe fallback.

---

## 9. Test Plan (Exact Commands)

### 9.1 Binary Audit Gate
```bash
nm -A /usr/local/lib/mumps/lib{s,d,c,z}mumps_cpmech.a /usr/local/lib/mumps/libmumps_common_cpmech.a \
  | awk '$2 ~ /^[bBdDgGsScCvV]$/' \
  | grep -E '__tpv_ptr_|STATIC_PTR|COLSCA_STATIC|ROWSCA_STATIC|SINGULAR_VALUES_STATIC|RHSINTR_STATIC|MUMPS_SYM_PERM|MUMPS_MAPPING'
# Expected Phase 0+: 0 matches
```

### 9.2 Standalone Fortran Reproducer (`examples/dorival_thread_safe.F`)
```bash
cd /home/dorival/01-Code/fortran/MUMPS_5.9.2/examples
make dorival_thread_safe
for omp in 1 4 8 13 16 32; do
  OMP_NUM_THREADS=$omp timeout 120 ./dorival_thread_safe 16 25 200 0 || exit 1
done
```

### 9.3 Standalone C `pthread` Churn Matrix (`examples/c_thread_safe_test.c` to be added during implementation)
Build against `/tmp/mumps-probe/lib/mumps` shared libraries (or installed `/usr/local/lib/mumps`) and run the cross-product matrix with ephemeral OS threads (spawned + joined per iteration so `libiomp5` team-reap stress is exercised):
- Caller OS threads: `T ∈ {1, 2, 8, 16, 32}` (both persistent pool callers and newly spawned/joined per-cycle threads)
- OpenMP runtime: `OMP_NUM_THREADS ∈ {1, 4, 8, 13, 16, 32}`
- Active instances per thread: `I ∈ {1, 8, 32}` (interleaving `JOB = -1, 1, 2, 3, -2`)
- Arithmetics & symmetries: `{smumps_c, dmumps_c, cmumps_c, zmumps_c}` $\times$ `SYM ∈ {0, 1, 2}` + iterative refinement (`ICNTL(10)=3`) and condition number (`ICNTL(11)=1`):
```bash
LD_LIBRARY_PATH=/tmp/mumps-probe/lib/mumps ./c_thread_safe_test \
  --threads 1,2,8,16,32 --omp 1,4,8,13,16,32 --instances 1,8,32 --ephemeral-churn
```

### 9.4 Integration Stress Gate (`russell_sparse` and `pmsim`)
```bash
# 1. Explicit stress test (10 threads x 150 iterations):
OMP_NUM_THREADS=16 cargo test --manifest-path ~/01-Code/rust/russell/Cargo.toml \
  -p russell_sparse --test test_mumps test_mumps_concurrent_stress -- --nocapture

# 2. 20/20 full-suite clean runs across critical OMP counts:
for omp in 1 13 16 32; do
  for run in $(seq 1 20); do
    OMP_NUM_THREADS=$omp cargo test --manifest-path ~/01-Code/rust/russell/Cargo.toml -p russell_sparse || exit 1
  done
done

for run in $(seq 1 20); do
  OMP_NUM_THREADS=16 cargo test --manifest-path ~/01-Code/rust/pmsim/Cargo.toml || exit 1
done

# 3. Kernel fault audit:
journalctl -k --since "30 minutes ago" | grep -i -E "segfault|traps|libdmumps|libiomp5"
# Expected: empty output
```

---

## 10. Performance Protocol and Baseline to Beat

### 10.1 Measurement Protocol
- Pin workers/benchmarks with `taskset -c 0-7` on an idle machine and run `hyperfine --warmup 3 --runs 20`.
- Baselines from Section 3.2 to compare against:
  - Per-call overhead: **≈ 7 µs** (`mumps-worker` baseline).
  - `cook_elastic_2d_fem` (192 dof, 20 steps, `OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`): **180.0 ms** (`mumps-worker`) / **176.5 ms** (stock unlocked single-thread).
  - `cook_elastic_2d_fem` (192 dof, 20 steps, `OMP_NUM_THREADS=8 MKL_NUM_THREADS=8`): **552.1 ms** (`mumps-worker`) / **450.5 ms** (stock unlocked).
  - `cook_elastic_3d_fem` (2187 dof, 20 steps, `OMP_NUM_THREADS=1`): **3972 ms**.

### 10.2 Commands & Bit-Identity Check
```bash
# 1. Bit-identity & CSV diff against stock MUMPS 5.9.2 at OMP=1 MKL=1:
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 cargo run --release --manifest-path ~/01-Code/rust/pmsim/Cargo.toml --bin ...
diff -u baseline_stock_omp1.csv patched_omp1.csv   # Required: 0 byte diff

# 2. 192-dof & 2187-dof latency benchmarks (OMP=1 and OMP=8):
taskset -c 0-7 hyperfine --warmup 3 --runs 20 \
  'OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ./target/release/bench_cook_2d' \
  'OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 ./target/release/bench_cook_2d' \
  'OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ./target/release/bench_cook_3d' \
  'OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 ./target/release/bench_cook_3d'

# 3. Phase 2 concurrent throughput benchmark (8 independent 2187-dof solves, OMP=4):
# Target: >= 4.0x wall-clock speedup over Phase 1 (K=1)
OMP_NUM_THREADS=4 MUMPS_NUM_WORKERS=8 ./c_concurrent_3d_bench
```

---

## 11. Build, Licensing, and Packaging

1. **License Compliance (CeCILL-C V1 — `LICENSE:1-13`, `doc/CeCILL-C_V1-en.txt`):**
   - MUMPS is licensed under **CeCILL-C V1**, a file-level weak-copyleft license compatible with downstream MIT/Apache-2.0 crates (`russell_sparse`, `pmsim`).
   - **Obligation:** Retain existing copyright notices, keep `LICENSE` and `doc/CeCILL-C_V1-*.txt`, and insert a prominent modification header in every touched/added file inside `src/` and `include/` stating the author, date, and nature of the thread-safety changes (Article 5.2). Downstream applications dynamically or statically linking `lib*mumps_cpmech.{a,so}` remain under their own license.
2. **Fork, Version String, and Tarball Packaging:**
   - **Canonical development fork:** `/home/dorival/01-Code/fortran/MUMPS_5.9.2` (`main` branch, tagged `v5.9.2-ts1` / `v5.9.2-ts2`). Keep `VERSION:1` as `MUMPS 5.9.2` (`build.bash:77`) so `MUMPS_VERSION` matches client headers while bumping `MUMPS_SCRIPT_VERSION` in `russell/zscripts/*-compile-mumps.bash` (e.g. `"5.9.2-ts1"`).
   - **Portable patch series:** Export `git diff 2da3511..v5.9.2-ts1 -- src/ include/ > zscripts/patches/mumps-5.9.2-thread-safe.patch`. The `russell/zscripts/{arch,debian,macos,windows}-compile-mumps.bash` installers download upstream `MUMPS_5.9.2.tar.gz`, apply the patch cleanly after `tar xzf`, copy `MakefileMKL.inc` / `Makefile.inc`, and build identically.
   - **Zero regression for stock / CI / crates.io users:** Upstream `MUMPS_5.9.2.tar.gz` continues to build unchanged when built directly; `russell_sparse` builds against stock system MUMPS (`libmumps-seq-dev`, Homebrew, etc.) automatically fall back to `mumps_lock.rs` whenever `mumps_is_thread_safe` is absent.

---

## 12. Risks, Open Questions for the Owner, and Kanban

### 12.1 Risks and Mitigations

| Risk | Severity | Mitigation |
| :--- | :--- | :--- |
| **`libiomp5` vs `libgomp` divergence** (`zscripts/MakefileMKL.inc` uses `ifx`/`icx` with `libiomp5`; `zscripts/Makefile.inc` uses `gfortran`/`clang` with `libgomp`). | High | Eliminate `!$OMP THREADPRIVATE` completely (Phase 0) and run OpenMP regions only inside persistent C `pthread` masters (`src/mumps_thread_safe.c`), testing both `./build.bash` and `./build.bash 0`. |
| **Nested caller OpenMP regions** (e.g., `examples/dorival_thread_safe.F`, where callers originate inside `!$OMP PARALLEL` with `OMP_MAX_ACTIVE_LEVELS=1`). | High | Dispatching from `libmumps_common` onto dedicated POSIX `pthread` masters detaches MUMPS from the caller's OpenMP team hierarchy, avoiding nested-region serialization (`OMP_GET_MAX_THREADS() = 1`) cleanly. |
| **MUMPS Fortran runtime state / I/O unit collisions** (if two Phase 2 workers log to `ICNTL(1..4)` Fortran units simultaneously). | Medium | Keep default quiet units or protect unit allocation (`MUMPS_FIND_UNIT` in `src/tools_common.F`) with a C mutex when error/diagnostic printing is enabled. |
| **Internal `DMUMPS_INTR_STRUC` layout growth** in Phase 2 (`src/[sdcz]mumps_intr_types.F:56-65`). | Low | `DMUMPS_INTR_STRUC` is private Fortran storage serialized via `id%INTR_ENCODING`; public C struct `DMUMPS_STRUC_C` (`include/[sdcz]mumps_c.h`) never changes layout. |

### 12.2 Numbered Open Questions for Owner Decision

1. **Stopping Gate (Phase 1 vs. Phase 2):** Approve **Option 1e** executing **Phase 0 + Phase 1** first (**1.5 person-weeks**) to eliminate client-side locking in `russell_sparse`, and decide whether **Phase 2** (**+4.5 person-weeks**, **6.0 person-weeks total**) should follow immediately for concurrent multi-instance speedup in `pmsim`.
2. **Fortran Entry Point Routing:** In Phase 1, should direct Fortran calls (`CALL DMUMPS(id)` as in `examples/dorival_thread_safe.F`) also route through `libmumps_common`'s persistent worker dispatcher (Option B, recommended), or only C/Rust calls (`[sdcz]mumps_c`)?
3. **Phase 2 OOC Serialization Policy:** Given that `russell_sparse` and `pmsim` run strictly in-core (`ICNTL(22)=0`), confirm that serializing OOC (`ICNTL(22)=1`) and Save/Restore (`JOB=7,8,-3`) on Worker 0 under a mutex in Phase 2 is acceptable.
4. **Target `MUMPS_NUM_WORKERS` Default for Phase 2:** Confirm default $K = 4$ (or $\min(8, N_{\text{cpus}} / \text{OMP\_NUM\_THREADS})$) overridable via `MUMPS_NUM_WORKERS` environment variable and `mumps_set_num_workers(int)`.

### 12.3 Kanban

| Done | In Progress | Backlog |
| :--- | :--- | :--- |
| - Complete static/`SAVE`/`COMMON`/`THREADPRIVATE` inventory (`nm -A` + source audit across `src/`, `libseq/`, `PORD/`, `libmetis.so`, `libGKlib.so`).<br>- Disassembly & root-cause verification of Fault 1a (`dmumps_f77_` at `0x274520`), Fault 1b/c (`MUMPS_STATIC_MAPPING`, `MUMPS_LOAD`, `mumps_c.c`, `mumps_common.c`, `DMUMPS_TMP_PTR`), Fault 2a (`DMUMPS_SOL_L0OMP_R` `__kmpc_threadprivate_cached`), and Fault 2b (`__kmp_reap_worker`).<br>- All 6 verification probes completed.<br>- Authored `PLAN-thread-safety.md`. | - **Phase 0 (0.5 pw), Kleo:** branch `ts-phase0` from `main` `d7b6b1f`, Ada reviewing; items 1–4 of §7; gate = zero symbols in the §9.1 audit, the five example runs and `test_mumps` green. | - Owner answers to Questions 1–4 (Decisions 1–4, §6).<br>- **Phase 1 (1.0 pw, blocked on `GO-PHASE-1`):** Implement `src/mumps_thread_safe.c` persistent OpenMP master worker + `mumps_is_thread_safe()` symbol; pass 20/20 `dorival_thread_safe`, C pthread churn, `russell_sparse`, and `pmsim` gates; patch `russell_sparse/build.rs` behind `#[cfg(mumps_thread_safe)]`.<br>- **Phase 2 (4.5 pw):** Encapsulate `MUMPS_STATIC_MAPPING`, `MUMPS_LOAD`, buffers, `FDM_F`, `LOCK_FOR_SCATTER`, `sol_aux.F`, & BLR state into `idintr`/TLS; enable $K$-worker pool + OOC mutex; verify $\ge 4\times$ concurrent speedup. |
