/*
 *
 *  This file is part of MUMPS 5.9.2, released
 *  on Fri Nov 13 10:00:55 UTC 2026
 *
 */
/*
 * C churn test for the persistent-worker dispatcher.
 *
 * 32 plain POSIX threads each own a DMUMPS_STRUC_C and run 20 complete
 * lifecycles (JOB = -1, 1, 2, 3, -2) of the same small SPD system.  Every
 * MUMPS call goes through the persistent worker installed by the library's
 * dispatcher, so the Fortran solver never sees two caller threads at once.
 *
 * The test matrix is the N x N tridiagonal A = tridiag(-1, 4, -1) and the
 * right-hand side is b = A * (1, ..., 1), hence the exact solution is
 * x = (1, ..., 1).  ICNTL(16) is left at 0, so MUMPS takes its own OpenMP
 * thread count from OMP_NUM_THREADS, as requested by the churn run.
 *
 * It takes no arguments, prints one summary line and exits non-zero on any
 * wrong solution or MUMPS error return.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mumps_compat.h"
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
# include <pthread.h>
#endif
#include "mpi.h"
#include "dmumps_c.h"

#define NTHREADS 32
#define NCYCLES  20
#define N        10
#define NNZ      (3 * N - 2)
#define JOB_INIT -1
#define JOB_END  -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1]

typedef struct
{
    int index;
    double max_err;
    int failures;
} thread_result;

static void
build_system(MUMPS_INT *irn, MUMPS_INT *jcn, double *a)
{
    int i, k = 0;

    for (i = 0; i < N; i++)
      {
        irn[k] = (MUMPS_INT) (i + 1);
        jcn[k] = (MUMPS_INT) (i + 1);
        a[k]   = 4.0;
        k++;
        if (i > 0)
          {
            irn[k] = (MUMPS_INT) (i + 1);
            jcn[k] = (MUMPS_INT) i;
            a[k]   = -1.0;
            k++;
          }
        if (i < N - 1)
          {
            irn[k] = (MUMPS_INT) (i + 1);
            jcn[k] = (MUMPS_INT) (i + 2);
            a[k]   = -1.0;
            k++;
          }
      }
}

static int
check_job(DMUMPS_STRUC_C *id, thread_result *r, int cycle, int job)
{
    if (id->infog[0] != 0 || id->info[0] != 0)
      {
        fprintf(stderr,
                "c_thread_safe_test: thread %d cycle %d JOB %d failed:"
                " INFOG(1)=%lld INFO(1)=%lld\n",
                r->index, cycle, job,
                (long long) id->infog[0], (long long) id->info[0]);
        r->failures++;
        return -1;
      }
    return 0;
}

#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
static void *
churn_worker(void *arg)
{
    thread_result *r = (thread_result *) arg;
    DMUMPS_STRUC_C id;
    MUMPS_INT irn[NNZ], jcn[NNZ];
    double a[NNZ], rhs[N];
    int cycle, i;

    build_system(irn, jcn, a);

    for (cycle = 0; cycle < NCYCLES; cycle++)
      {
        memset(&id, 0, sizeof(id));
        id.comm_fortran = USE_COMM_WORLD;
        id.par = 1;
        id.sym = 0;
        id.job = JOB_INIT;
        dmumps_c(&id);
        if (check_job(&id, r, cycle, JOB_INIT) != 0)
          {
            return NULL;
          }

        /* Quiet run; ICNTL(16) stays 0, i.e. OMP_NUM_THREADS decides. */
        id.ICNTL(1) = -1;
        id.ICNTL(2) = -1;
        id.ICNTL(3) = -1;
        id.ICNTL(4) = 0;
        id.ICNTL(16) = 0;

        id.n = N;
        id.nnz = (MUMPS_INT8) NNZ;
        id.irn = irn;
        id.jcn = jcn;
        id.a = a;
        id.rhs = rhs;

        id.job = 1;
        dmumps_c(&id);
        if (check_job(&id, r, cycle, 1) != 0)
          {
            return NULL;
          }

        id.job = 2;
        dmumps_c(&id);
        if (check_job(&id, r, cycle, 2) != 0)
          {
            return NULL;
          }

        /* b = A * ones: 3 at both ends, 2 in the interior. */
        for (i = 0; i < N; i++)
          {
            rhs[i] = 2.0;
          }
        rhs[0] = 3.0;
        rhs[N - 1] = 3.0;

        id.job = 3;
        dmumps_c(&id);
        if (check_job(&id, r, cycle, 3) != 0)
          {
            return NULL;
          }

        for (i = 0; i < N; i++)
          {
            double err = fabs(rhs[i] - 1.0);
            if (err > r->max_err)
              {
                r->max_err = err;
              }
            if (err > 1e-12)
              {
                fprintf(stderr,
                        "c_thread_safe_test: thread %d cycle %d wrong"
                        " solution at row %d: |x-1|=%.3e\n",
                        r->index, cycle, i + 1, err);
                r->failures++;
                return NULL;
              }
          }

        id.job = JOB_END;
        dmumps_c(&id);
        if (check_job(&id, r, cycle, JOB_END) != 0)
          {
            return NULL;
          }
      }

    return NULL;
}
#endif /* ! MUMPS_WIN32 && ! WITHOUT_PTHREAD */

int
main(int argc, char **argv)
{
#if defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD)
    (void) argc;
    (void) argv;
    printf("c_thread_safe_test: skipped (pthread support is not available)\n");
    return 0;
#else
    pthread_t tids[NTHREADS];
    thread_result results[NTHREADS];
    int i, created = 0, failures = 0;
    double max_err = 0.0;
    int ierr;

    (void) argc;
    (void) argv;
    /* MUMPS reads OMP_NUM_THREADS itself because ICNTL(16) is left at 0. */
    (void) getenv("OMP_NUM_THREADS");

    ierr = MPI_Init(&argc, &argv);
    if (ierr != 0)
      {
        fprintf(stderr, "c_thread_safe_test: MPI_Init failed (%d)\n", ierr);
        return 1;
      }

    for (i = 0; i < NTHREADS; i++)
      {
        results[i].index = i;
        results[i].max_err = 0.0;
        results[i].failures = 0;
        if (pthread_create(&tids[i], NULL, churn_worker, &results[i]) != 0)
          {
            fprintf(stderr,
                    "c_thread_safe_test: pthread_create failed for thread %d\n",
                    i);
            failures++;
            break;
          }
        created++;
      }

    for (i = 0; i < created; i++)
      {
        (void) pthread_join(tids[i], NULL);
        if (results[i].max_err > max_err)
          {
            max_err = results[i].max_err;
          }
        failures += results[i].failures;
      }

    ierr = MPI_Finalize();
    if (ierr != 0)
      {
        fprintf(stderr, "c_thread_safe_test: MPI_Finalize failed (%d)\n", ierr);
        failures++;
      }

    printf("threads=%d cycles=%d max|x-1|=%.3e failures=%d\n",
           NTHREADS, NCYCLES, max_err, failures);

    return (failures == 0) ? 0 : 1;
#endif
}
