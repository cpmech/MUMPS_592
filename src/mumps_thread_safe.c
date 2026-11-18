/*
 *
 *  This file is part of MUMPS 5.9.2, released
 *  on Fri Nov 13 10:00:55 UTC 2026
 *
 *
 *  Copyright 1991-2026 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
/*
 *  Modified 2026-11-18 by Dorival Pedroso: thread-safe entry points, per the
 *  accompanying cover letter.  The CeCILL-C notice above is unchanged.
 */
/*
 * Persistent single-worker dispatcher.
 *
 * The C entry points of every arithmetic ([sdcz]mumps_c) hand their whole
 * call, structure included, to one process-wide worker thread.  The request
 * lives on the caller's stack and the caller does not return until the
 * worker has set its "done" flag, so the frame stays alive for the worker.
 * Callers that arrive concurrently queue on the request mutex: exactly one
 * request is outstanding at a time, which is the intended single-worker
 * design.
 *
 * The worker is created lazily by the first dispatch.  If creation (or the
 * TLS key needed to detect re-entry) fails, the dispatcher becomes
 * permanently unavailable and says so once on stderr; it never pretends a
 * call ran.  There is deliberately no watchdog: if the worker dies, callers
 * wait forever.  That is the documented trade of running a solver on one
 * persistent thread.
 *
 * A call made from the worker itself runs inline, so a future nested entry
 * point cannot try to queue behind the worker that is already executing.
 *
 * Thread safety needs pthreads.  Windows builds and builds configured with
 * -DWITHOUT_PTHREAD run every call inline and report 0 from
 * mumps_is_thread_safe().  Real MPI builds also run inline: thread safety is
 * provided only for the sequential libmpiseq builds, because MPI requires
 * every call on the MPI_Init thread.  Which build this is comes from the
 * mpi.h that the build's include path selects: only libseq/mpi.h defines
 * MUMPS_MPI_H, a real mpi.h does not.
 */
#if defined(__linux__)
# define _GNU_SOURCE 1
#endif
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mumps_common.h"

#if defined(_OPENMP)
# if defined(__has_include)
#  if __has_include(<omp.h>)
#   include <omp.h>
#  else
/* Some compilers accept -fopenmp without shipping the runtime header. */
extern int omp_get_max_threads(void);
extern void omp_set_num_threads(int);
#  endif
# else
#  include <omp.h>
# endif
#endif

#if ! defined(MUMPS_MPI) && defined(MUMPS_MPI_H) && \
    ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
# include <pthread.h>
#endif

#ifndef MUMPS_THREAD_SAFE_VERSION
# define MUMPS_THREAD_SAFE_VERSION 1
#endif

#if ! defined(MUMPS_MPI) && defined(MUMPS_MPI_H) && \
    ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
typedef struct mumps_ts_request
{
    mumps_ts_fn fn;
    void *arg;
    int caller_omp_max;
    int done;
} mumps_ts_request;

static pthread_once_t mumps_ts_once = PTHREAD_ONCE_INIT;
static pthread_mutex_t mumps_ts_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t mumps_ts_work = PTHREAD_COND_INITIALIZER;
static pthread_cond_t mumps_ts_done = PTHREAD_COND_INITIALIZER;
static pthread_key_t mumps_ts_key;
static char mumps_ts_worker_marker[1];
static mumps_ts_request *volatile mumps_ts_slot = NULL;
/* 0 = not tried, 1 = ready, -1 = failed (permanently unavailable) */
static volatile int mumps_ts_worker_state = 0;
static int mumps_ts_failure_reported = 0;

static void
mumps_ts_report_failure(const char *what, int err)
{
    if (!mumps_ts_failure_reported)
      {
        mumps_ts_failure_reported = 1;
        fprintf(stderr,
                "MUMPS: thread-safe dispatcher unavailable: %s failed (%d)\n",
                what, err);
      }
}

static void mumps_ts_init(void);

/*
 * Run one request on the worker.  ICNTL(16) <= 0 means the caller keeps the
 * OpenMP thread count it configured with omp_set_num_threads(); apply that
 * count on the worker for the duration of the call and restore the worker's
 * own setting afterwards.  caller_omp_max == 0 leaves the worker alone.
 */
static void
mumps_ts_run(mumps_ts_fn fn, void *arg, int caller_omp_max)
{
#if defined(_OPENMP)
    if (caller_omp_max > 0)
      {
        int saved = omp_get_max_threads();
        omp_set_num_threads(caller_omp_max);
        fn(arg);
        omp_set_num_threads(saved);
      }
    else
      {
        fn(arg);
      }
#else
    /* Without OpenMP there is no per-thread setting to apply. */
    (void) caller_omp_max;
    fn(arg);
#endif
}

static void *
mumps_ts_main(void *unused)
{
    (void) unused;
    /*
     * Mark this thread as the dispatcher worker.  pthread_once guarantees
     * that the key exists before the worker is created.
     */
    (void) pthread_setspecific(mumps_ts_key, (void *) mumps_ts_worker_marker);
#if defined(__linux__)
    {
        /*
         * Linux limits the comm name to 15 characters plus the terminating
         * NUL, so the requested name cannot be stored in full.  Try it
         * verbatim first and fall back to the 15-character prefix.
         */
        static const char worker_name[] = "mumps-omp-master-0";
        static const char worker_name_short[] = "mumps-omp-maste";
        if (pthread_setname_np(pthread_self(), worker_name) != 0)
          {
            (void) pthread_setname_np(pthread_self(), worker_name_short);
          }
    }
#endif
    for (;;)
      {
        mumps_ts_request *req;

        pthread_mutex_lock(&mumps_ts_mutex);
        while (mumps_ts_slot == NULL)
          {
            pthread_cond_wait(&mumps_ts_work, &mumps_ts_mutex);
          }
        /*
         * Copy the request pointer out but leave the slot posted until the
         * work is done: other callers must keep queueing behind us.
         */
        req = mumps_ts_slot;
        pthread_mutex_unlock(&mumps_ts_mutex);

        mumps_ts_run(req->fn, req->arg, req->caller_omp_max);

        pthread_mutex_lock(&mumps_ts_mutex);
        req->done = 1;
        mumps_ts_slot = NULL;
        /* Broadcast: callers may be waiting for their own completion or for
         * the slot to become free. */
        pthread_cond_broadcast(&mumps_ts_done);
        pthread_mutex_unlock(&mumps_ts_mutex);
      }

    return NULL; /* not reached */
}

static void
mumps_ts_init(void)
{
    pthread_t worker;
    int err;

    err = pthread_key_create(&mumps_ts_key, NULL);
    if (err != 0)
      {
        mumps_ts_report_failure("pthread_key_create", err);
        mumps_ts_worker_state = -1;
        return;
      }

    err = pthread_create(&worker, NULL, mumps_ts_main, NULL);
    if (err != 0)
      {
        mumps_ts_report_failure("pthread_create", err);
        (void) pthread_key_delete(mumps_ts_key);
        mumps_ts_worker_state = -1;
        return;
      }

    (void) pthread_detach(worker);
    mumps_ts_worker_state = 1;
}
#endif /* ! MUMPS_MPI && MUMPS_MPI_H && ! MUMPS_WIN32 &&
          ! WITHOUT_PTHREAD */

int
mumps_ts_dispatch(int arith, void *struc, mumps_ts_fn fn, int caller_omp_max)
{
#if defined(MUMPS_MPI) || ! defined(MUMPS_MPI_H)
    /*
     * A real MPI build must keep every call on the thread that called
     * MPI_Init (MPI_THREAD_SINGLE/MPI_THREAD_FUNNELED).  The header that the
     * build's include path selects is what tells us which build this is:
     * libseq/mpi.h defines MUMPS_MPI_H and a real mpi.h does not, while
     * -DMUMPS_MPI remains an explicit override.  Run inline and never create
     * a worker.
     */
    (void) arith;
    (void) caller_omp_max;
    fn(struc);
    return 0;
#elif defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD)
    /* Thread safety needs pthreads: run inline and say so at run time. */
    (void) arith;
    (void) caller_omp_max;
    fn(struc);
    return 0;
#else
    int err;

    /*
     * arith is accepted for a future worker pool that will key the worker
     * on the arithmetic; this version uses one worker for all of them.
     */
    (void) arith;

    err = pthread_once(&mumps_ts_once, mumps_ts_init);
    if (err != 0)
      {
        mumps_ts_report_failure("pthread_once", err);
        mumps_ts_worker_state = -1;
        return -1;
      }
    if (mumps_ts_worker_state != 1)
      {
        return -1;
      }

    /*
     * A nested call made from the worker itself runs inline.  Queueing it
     * would deadlock on the single worker that is already busy with us.
     */
    if (pthread_getspecific(mumps_ts_key) != NULL)
      {
        mumps_ts_run(fn, struc, caller_omp_max);
        return 0;
      }

    {
        mumps_ts_request req;

        req.fn = fn;
        req.arg = struc;
        req.caller_omp_max = caller_omp_max;
        req.done = 0;

        pthread_mutex_lock(&mumps_ts_mutex);
        /* Wait for the worker to clear any previous request. */
        while (mumps_ts_slot != NULL)
          {
            pthread_cond_wait(&mumps_ts_done, &mumps_ts_mutex);
          }
        mumps_ts_slot = &req;
        pthread_cond_signal(&mumps_ts_work);
        /* Do not return before the worker has finished with our stack. */
        while (!req.done)
          {
            pthread_cond_wait(&mumps_ts_done, &mumps_ts_mutex);
          }
        pthread_mutex_unlock(&mumps_ts_mutex);
    }

    return 0;
#endif
}

int
MUMPS_CALL
mumps_ts_dispatch_f77(int *arith, void *struc, mumps_ts_fn fn,
                      int *caller_omp_max)
{
#if defined(MUMPS_MPI) || ! defined(MUMPS_MPI_H)
    /*
     * A real MPI build must keep every call on the thread that called
     * MPI_Init (MPI_THREAD_SINGLE/MPI_THREAD_FUNNELED).  The header that the
     * build's include path selects is what tells us which build this is:
     * libseq/mpi.h defines MUMPS_MPI_H and a real mpi.h does not, while
     * -DMUMPS_MPI remains an explicit override.  Run inline and never create
     * a worker.
     */
    (void) arith;
    (void) caller_omp_max;
    fn(struc);
    return 0;
#elif defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD)
    /* Thread safety needs pthreads: run inline and say so at run time. */
    (void) arith;
    (void) caller_omp_max;
    fn(struc);
    return 0;
#else
    return mumps_ts_dispatch(*arith, struc, fn, *caller_omp_max);
#endif
}

/*
 * Fortran-callable spelling of the same entry point.  F_SYMBOL applies the
 * compiler's decoration (a trailing underscore for gfortran/ifort with
 * -DAdd_), which is what a plain EXTERNAL call from [SDCZ]MUMPS expects;
 * calling mumps_ts_dispatch_f77 directly would leave that underscore off
 * and fail to link.  With no decoration the two names coincide, so the
 * alias is unnecessary and would be a redefinition.
 */
#if defined(UPPER) || defined(MUMPS_WIN32) || defined(Add_) || defined(Add__)
int
MUMPS_TS_DISPATCH_F77(int *arith, void *struc, mumps_ts_fn fn,
                      int *caller_omp_max)
{
    return mumps_ts_dispatch_f77(arith, struc, fn, caller_omp_max);
}
#endif

int
MUMPS_CALL
mumps_is_thread_safe(void)
{
#if ! defined(MUMPS_MPI) && defined(MUMPS_MPI_H) && \
    ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
    return 1;
#else
    /* Thread safety needs pthreads, and MPI calls stay on the MPI thread. */
    return 0;
#endif
}
