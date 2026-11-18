C
C  This file is part of MUMPS 5.9.2, released
C  on Fri Nov 13 10:00:55 UTC 2026
C
C  THREAD_SAFE_TEST: several MUMPS instances running from several threads
C
C  This is the reproducer for the MUMPS thread-safety defect.  Every thread of an
C  OpenMP team owns its own DMUMPS instance and churns the whole lifecycle
C  (JOB = -1, 1, 2, 3, -2) concurrently with the others, then verifies the
C  solution against the exact one.  A stock MUMPS keeps process-global state in
C  its modules and its OpenMP runtime, so the runs are expected to crash, hang or
C  return wrong results; a patched MUMPS must pass with every thread count.
C
C  Usage:
C
C     ./thread_safe_test [nthreads] [niter] [n] [locked] [nested]
C                           [mumps_threads]
C
C     nthreads       caller threads (default 8)
C     niter          solver cycles per thread (default 20), each cycle is one
C                    full JOB = -1, 1, 2, 3, -2 sequence
C     n              order of the tridiagonal test matrix (default 100)
C     locked         0  no synchronization at all (default)
C                    1  a named critical section around every CALL DMUMPS,
C                       i.e. the calls are serialized but the phases of
C                       different instances still interleave
C                    2  a named critical section around the whole cycle,
C                       i.e. one instance at a time end to end
C     nested         1 enables OpenMP nesting (OMP_SET_MAX_ACTIVE_LEVELS(2)),
C                    so every caller is the master of its own MUMPS team;
C                    0 keeps MUMPS's parallel regions at one thread inside
C                    our parallel region (default 0)
C     mumps_threads  ICNTL(16) forwarded to every instance.  When omitted it
C                    is 1 with nested = 0 (MUMPS's own OpenMP team cannot be
C                    entered from inside our parallel region, so its internal
C                    parallelism has to be off; 8 nested-in-1 deadlocks in
C                    dmumps_fac_omp_m) and 0, the OpenMP default, with
C                    nested = 1
C
C  Build: this example needs an OpenMP-enabled compiler and is NOT part of the
C  `make d` target: a broken example must not stop the library build, and a
C  platform makefile that compiles MUMPS without -fopenmp (macOS, with
C  Accelerate) cannot link omp_get_thread_num and friends.  On a machine whose
C  toolchain has OpenMP, build it on demand once the libraries exist:
C
C     make dshared && make thread_safe_test
C
C  The test matrix is A = tridiag(-1, 4, -1) and the right-hand side is
C  b = A * (1, ..., 1), so the exact solution is x = (1, ..., 1).  Each cycle
C  checks max |x - 1| <= 1e-10; a wrong solution counts as a failure even when
C  the process does not crash, and a MUMPS error return (INFOG(1) < 0) is
C  reported with its JOB.
C
C  With locked = 0 or 1 the process may deadlock instead of failing (MUMPS or
C  its OpenMP runtime waiting on a lock a corrupted instance never releases),
C  so wrap a run in `timeout`: the exit code 124 is itself a reproduction.
C
C  Observed with a stock MUMPS 5.9.2 build (MKL, OMP_NUM_THREADS = 8,
C  10 cycles, n = 100):
C
C     8 threads, locked = 0   6/6 runs fail: SIGSEGV, an allocatable-array
C                             abort, or a MUMPS error at JOB = 1
C     2 threads, locked = 0   6/6 runs fail: SIGSEGV, the same abort, or
C                             INFOG(1) = -13 at JOB = 1
C     8 threads, locked = 2   6/6 runs pass
C
C  The OpenMP-team failure (threadprivate caches past libiomp5's
C  __kmp_default_tp_capacity, __kmp_reap_worker) needs plain OS caller threads,
C  which an OpenMP team does not provide: nested = 1 with 10 callers at
C  OMP_NUM_THREADS = 16 still passes here.  This example reproduces the
C  concurrent-entry corruption.
C
C  Build against an in-tree or installed library (from this directory):
C
C     ifx -O -qopenmp -I../include -I../libseq
C         thread_safe_test.F -o thread_safe_test
C         -L../lib -ldmumps -lsmumps
C         -lmumps_common -lpord -lmpiseq -lmetis
C         -L$MKLROOT/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread
C         -lmkl_core -lpthread
C
C  or `make thread_safe_test` after building the library tree.
C
C  Examples:
C
C     timeout 60 ./thread_safe_test 8 20 100       # the raw defect
C     timeout 60 ./thread_safe_test 8 20 100 1     # per-call mutex
C     timeout 60 ./thread_safe_test 8 20 100 2     # per-cycle mutex
C     OMP_NUM_THREADS=16 timeout 60 ./thread_safe_test 10 10 100 2 1
C                                            # nested teams, tp_capacity
C
      PROGRAM THREAD_SAFE_TEST
      USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER NTH, NITER, N, LOCKED, NESTED, MTH
      INTEGER I, IERR, NF
      DOUBLE PRECISION, ALLOCATABLE :: EMAX(:)

      NTH = 8
      NITER = 20
      N = 100
      LOCKED = 0
      NESTED = 0
      MTH = -1
      CALL GET_INT_ARG(1, NTH)
      CALL GET_INT_ARG(2, NITER)
      CALL GET_INT_ARG(3, N)
      CALL GET_INT_ARG(4, LOCKED)
      CALL GET_INT_ARG(5, NESTED)
      CALL GET_INT_ARG(6, MTH)
      IF (NTH < 1 .OR. NITER < 1 .OR. N < 2) THEN
         WRITE(6,'(A)') ' invalid arguments'
         CALL FLUSH(6)
         CALL EXIT(2)
      END IF
      IF (LOCKED < 0 .OR. LOCKED > 2) THEN
         WRITE(6,'(A)') ' locked must be 0, 1 or 2'
         CALL FLUSH(6)
         CALL EXIT(2)
      END IF
      IF (MTH < 0) THEN
         IF (NESTED == 0) THEN
            MTH = 1
         ELSE
            MTH = 0
         END IF
      END IF

      WRITE(6,'(A,I3,A,I4,A,I5,A,I1,A,I1,A,I4)')
     &    ' threads = ', NTH, '  cycles = ', NITER, '  n = ', N,
     &    '  locked = ', LOCKED, '  nested = ', NESTED,
     &    '  ICNTL(16) = ', MTH

      IF (NESTED == 1) THEN
         CALL OMP_SET_NESTED(.TRUE.)
         CALL OMP_SET_MAX_ACTIVE_LEVELS(2)
      END IF

      CALL MPI_INIT(IERR)
      ALLOCATE(EMAX(NTH))
      DO I = 1, NTH
         EMAX(I) = 0.0D0
      END DO
      NF = 0

!$OMP PARALLEL NUM_THREADS(NTH) SHARED(EMAX,NF)
      CALL WORKER(OMP_GET_THREAD_NUM(), N, NITER, LOCKED, MTH,
     &            EMAX, NF)
!$OMP END PARALLEL

      CALL MPI_FINALIZE(IERR)

      DO I = 1, NTH
         WRITE(6,'(A,I3,A,ES12.4)') ' thread ', I-1,
     &       '  max |x-1| = ', EMAX(I)
      END DO
      IF (NF > 0) THEN
         WRITE(6,'(A,I6)') ' FAILED: cycles reporting errors = ', NF
         CALL FLUSH(6)
         CALL EXIT(1)
      END IF
      WRITE(6,'(A)') ' PASSED'
      CALL FLUSH(6)
      CALL EXIT(0)
      END

C-----------------------------------------------------------------------
C  WORKER: one MUMPS instance, NITER cycles, solution verified
C-----------------------------------------------------------------------
      SUBROUTINE WORKER(T, N, NITER, LOCKED, MTH, EMAX, NF)
      IMPLICIT NONE
      INCLUDE 'dmumps_struc.h'
C     T is the value of OMP_GET_THREAD_NUM(), i.e. omp_integer_kind (4
C     bytes) even when this file is compiled with -i8; a default INTEGER
C     dummy would read 8 bytes from the 4-byte actual argument.
      INTEGER(4) T
      INTEGER N, NITER, LOCKED, MTH, NF
      DOUBLE PRECISION EMAX(*)
      TYPE(DMUMPS_STRUC) ID
      INTEGER K, JOBFAIL, NFAIL
      DOUBLE PRECISION ERR, ETHIS

      ETHIS = 0.0D0
      JOBFAIL = 0
      NFAIL = 0
      DO K = 1, NITER
         IF (LOCKED == 2) THEN
!$OMP CRITICAL(MUMPS_CYCLE)
            CALL DO_CYCLE(T, K, ID, N, MTH, LOCKED, ERR, JOBFAIL)
!$OMP END CRITICAL(MUMPS_CYCLE)
         ELSE
            CALL DO_CYCLE(T, K, ID, N, MTH, LOCKED, ERR, JOBFAIL)
         END IF
         ETHIS = MAX(ETHIS, ERR)
         IF (JOBFAIL /= 0) EXIT
         IF (ERR > 1.0D-10) THEN
            CALL REPORT2(T, K, ERR)
            NFAIL = 1
            EXIT
         END IF
      END DO

      EMAX(T+1) = ETHIS
      IF (JOBFAIL + NFAIL > 0) THEN
!$OMP CRITICAL(COUNT)
         NF = NF + 1
!$OMP END CRITICAL(COUNT)
      END IF
      END

C-----------------------------------------------------------------------
C  DO_CYCLE: JOB = -1, 1, 2, 3, -2 for one instance, one time
C  JOB = -1 resets the structure, so N, NNZ, the matrix arrays and
C  ICNTL(16) are set afterwards, as in the MUMPS examples
C-----------------------------------------------------------------------
      SUBROUTINE DO_CYCLE(T, K, ID, N, MTH, LOCKED, ERR, JOBFAIL)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'
      INTEGER(4) T
      INTEGER K, N, MTH, LOCKED, JOBFAIL
      TYPE(DMUMPS_STRUC) ID
      DOUBLE PRECISION ERR
      INTEGER I, J, NNZ, ISTAT
      LOGICAL ALLOCED

      ERR = 0.0D0
      JOBFAIL = 0
      ALLOCED = .FALSE.

      ID%COMM = MPI_COMM_WORLD
      ID%SYM = 0
      ID%PAR = 1
      ID%JOB = -1
      CALL RUN_JOB(ID, LOCKED)
      IF (ID%INFOG(1) < 0) THEN
         CALL REPORT(T, K, -1, ID%INFOG(1))
         JOBFAIL = 1
         RETURN
      END IF

C     Quiet: ICNTL(4) = 0 level of printing, ICNTL(1) = error unit
      ID%ICNTL(4) = 0
      ID%ICNTL(1) = 6
      IF (MTH >= 0) ID%ICNTL(16) = MTH

      NNZ = 3 * N - 2
      ALLOCATE(ID%IRN(NNZ), ID%JCN(NNZ), ID%A(NNZ), ID%RHS(N))
      ALLOCED = .TRUE.

C     Tridiagonal matrix: 4 on the diagonal, -1 on both subdiagonals
      J = 0
      DO I = 1, N
         J = J + 1
         ID%IRN(J) = I
         ID%JCN(J) = I
         ID%A(J) = 4.0D0
         IF (I > 1) THEN
            J = J + 1
            ID%IRN(J) = I
            ID%JCN(J) = I - 1
            ID%A(J) = -1.0D0
         END IF
         IF (I < N) THEN
            J = J + 1
            ID%IRN(J) = I
            ID%JCN(J) = I + 1
            ID%A(J) = -1.0D0
         END IF
      END DO

      ID%N = N
      ID%NNZ = NNZ

      ID%JOB = 1
      CALL RUN_JOB(ID, LOCKED)
      IF (ID%INFOG(1) < 0) THEN
         CALL REPORT(T, K, 1, ID%INFOG(1))
         JOBFAIL = 1
         GOTO 900
      END IF

      ID%JOB = 2
      CALL RUN_JOB(ID, LOCKED)
      IF (ID%INFOG(1) < 0) THEN
         CALL REPORT(T, K, 2, ID%INFOG(1))
         JOBFAIL = 1
         GOTO 900
      END IF

C     b = A * 1, i.e. 3 at the ends and 2 inside; the solve overwrites RHS
      DO I = 1, N
         ID%RHS(I) = 2.0D0
      END DO
      ID%RHS(1) = 3.0D0
      ID%RHS(N) = 3.0D0

      ID%JOB = 3
      CALL RUN_JOB(ID, LOCKED)
      IF (ID%INFOG(1) < 0) THEN
         CALL REPORT(T, K, 3, ID%INFOG(1))
         JOBFAIL = 1
         GOTO 900
      END IF

      DO I = 1, N
         ERR = MAX(ERR, ABS(ID%RHS(I) - 1.0D0))
      END DO

      ID%JOB = -2
      CALL RUN_JOB(ID, LOCKED)
      IF (ID%INFOG(1) < 0) THEN
         CALL REPORT(T, K, -2, ID%INFOG(1))
         JOBFAIL = 1
      END IF

  900 CONTINUE
      IF (ALLOCED) THEN
         DEALLOCATE(ID%IRN, ID%JCN, ID%A, ID%RHS, STAT=ISTAT)
      END IF
      END

C-----------------------------------------------------------------------
C  RUN_JOB: one CALL DMUMPS; the per-call mutex is locked = 1 only
C-----------------------------------------------------------------------
      SUBROUTINE RUN_JOB(ID, LOCKED)
      IMPLICIT NONE
      INCLUDE 'dmumps_struc.h'
      TYPE(DMUMPS_STRUC) ID
      INTEGER LOCKED
      IF (LOCKED == 1) THEN
!$OMP CRITICAL(MUMPS_CALL)
         CALL DMUMPS(ID)
!$OMP END CRITICAL(MUMPS_CALL)
      ELSE
         CALL DMUMPS(ID)
      END IF
      END

C-----------------------------------------------------------------------
C  REPORT: a MUMPS error return, printed one line at a time
C-----------------------------------------------------------------------
      SUBROUTINE REPORT(T, K, JOB, INFO)
      IMPLICIT NONE
      INTEGER(4) T
      INTEGER K, JOB, INFO
!$OMP CRITICAL(IO)
      WRITE(6,'(A,I3,A,I4,A,I3,A,I9)') ' thread ', T, ' cycle ', K,
     &    ' JOB ', JOB, ' failed, INFOG(1) = ', INFO
!$OMP END CRITICAL(IO)
      END

C-----------------------------------------------------------------------
C  REPORT2: a wrong solution, printed one line at a time
C-----------------------------------------------------------------------
      SUBROUTINE REPORT2(T, K, ERR)
      IMPLICIT NONE
      INTEGER(4) T
      INTEGER K
      DOUBLE PRECISION ERR
!$OMP CRITICAL(IO)
      WRITE(6,'(A,I3,A,I4,A,ES12.4)') ' thread ', T, ' cycle ', K,
     &    ' wrong solution, max |x-1| = ', ERR
!$OMP END CRITICAL(IO)
      END

C-----------------------------------------------------------------------
C  GET_INT_ARG: reads the POS-th command line argument; leaves VALUE as
C  it is when the argument is absent or empty
C-----------------------------------------------------------------------
      SUBROUTINE GET_INT_ARG(POS, VALUE)
      IMPLICIT NONE
      INTEGER POS, VALUE
      CHARACTER(LEN=64) :: ARG
      IF (POS > COMMAND_ARGUMENT_COUNT()) RETURN
      CALL GET_COMMAND_ARGUMENT(POS, ARG)
      IF (LEN_TRIM(ARG) == 0) RETURN
      READ(ARG,*) VALUE
      END
