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
 *  ===========================================================================
 *  Modification for thread safety (Phase 1), 2026-11-17
 *  Author: Dorival Pedroso (Kleo completed the edit)
 *  Nature: added the persistent single-worker MUMPS dispatcher and the
 *          mumps_is_thread_safe capability symbol.
 *
 *  MUMPS_THREAD_SAFE_VERSION is 1 (see mumps_common.h as well).
 *  ===========================================================================
 */
/*
 * Persistent single-worker dispatcher.
 *
 * The C entry points of every arithmetic ([sdcz]mumps_c) hand their whole
 * call, structure included, to one process-wide worker thread.  The request
 * lives on the caller's stack and the caller does not return until the
 * worker has set its "done" flag, so the frame stays alive for the worker.
 * Callers that arrive concurrently queue on the request mutex: exactly one
 * request is outstanding at a time, which is the intended Phase 1 design.
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
 */
#if defined(__linux__)
# define _GNU_SOURCE 1
#endif
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include "mumps_common.h"

#ifndef MUMPS_THREAD_SAFE_VERSION
# define MUMPS_THREAD_SAFE_VERSION 1
#endif

typedef struct mumps_ts_request
{
    mumps_ts_fn fn;
    void *arg;
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

        req->fn(req->arg);

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

int
mumps_ts_dispatch(int arith, void *struc, mumps_ts_fn fn)
{
    int err;

    /*
     * arith is accepted for the Phase 2 worker pool, which will key the
     * worker on the arithmetic; Phase 1 uses one worker for all of them.
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
        fn(struc);
        return 0;
      }

    {
        mumps_ts_request req;

        req.fn = fn;
        req.arg = struc;
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
}

int
mumps_ts_dispatch_f77(int arith, void *struc, mumps_ts_fn fn)
{
    return mumps_ts_dispatch(arith, struc, fn);
}

/*
 * Fortran-callable spelling of the same entry point.  F_SYMBOL applies the
 * compiler's decoration (a trailing underscore for gfortran/ifort with
 * -DAdd_), which is what a plain EXTERNAL call from [SDCZ]MUMPS expects;
 * calling mumps_ts_dispatch_f77 directly would leave that underscore off
 * and fail to link.
 */
int
MUMPS_TS_DISPATCH_F77(int arith, void *struc, mumps_ts_fn fn)
{
    return mumps_ts_dispatch(arith, struc, fn);
}

int
mumps_is_thread_safe(void)
{
    return 1;
}
