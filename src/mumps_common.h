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
#ifndef MUMPS_COMMON_H
#define MUMPS_COMMON_H
#include "mumps_compat.h"
#include "mumps_c_types.h"
/**
 * F_SYMBOL is a macro that converts a couple (lower case symbol, upper
 * case symbol) into the symbol defined by the compiler convention.
 * Example: For MUMPS_XXX, first define
 *   #define MUMPS_XXX F_SYMBOL(xxx,XXX) and then use
 *   MUMPS_XXX in the code to get rid of any symbol convention annoyance.
 *
 * NB: We need to provide both upper and lower case versions because to our
 *     knowledge, there is no way to perform the conversion with CPP
 *     directives only.
 */
#if defined(UPPER) || defined(MUMPS_WIN32)
# define F_SYMBOL(lower_case,upper_case) MUMPS_##upper_case
#elif defined(Add_)
# define F_SYMBOL(lower_case,upper_case) mumps_##lower_case##_
#elif defined(Add__)
# define F_SYMBOL(lower_case,upper_case) mumps_##lower_case##__
#else
# define F_SYMBOL(lower_case,upper_case) mumps_##lower_case
#endif
/*
 * Serializes the [SDCZ]MUMPS_F77 JOB=-1 instance-table management.
 * The Fortran wrappers call these through plain EXTERNAL calls, so the
 * symbols follow the usual F_SYMBOL naming convention.
 */
#define MUMPS_F77_LOCK F_SYMBOL(f77_lock,F77_LOCK)
void MUMPS_CALL MUMPS_F77_LOCK(void);
#define MUMPS_F77_UNLOCK F_SYMBOL(f77_unlock,F77_UNLOCK)
void MUMPS_CALL MUMPS_F77_UNLOCK(void);
#define MUMPS_ICOPY_32TO64_64C_IP_C \
    F_SYMBOL(icopy_32to64_64c_ip_c,ICOPY_32TO64_64C_IP_C)
void MUMPS_CALL
MUMPS_ICOPY_32TO64_64C_IP_C(MUMPS_INT *inouttab, MUMPS_INT8 *sizetab);
#define MUMPS_ICOPY_64TO32_64C_IP_C \
    F_SYMBOL(icopy_64to32_64c_ip_c,ICOPY_64TO32_64C_IP_C)
void MUMPS_CALL
MUMPS_ICOPY_64to32_64C_IP_C(MUMPS_INT8 *inouttab, MUMPS_INT8 *sizetab);
#define MUMPS_MALLOC_C \
    F_SYMBOL(malloc_c,MALLOC_C)
void MUMPS_CALL
MUMPS_MALLOC_C(MUMPS_INT8 *address, MUMPS_INT8 *size);
#define MUMPS_FREE_C \
    F_SYMBOL(free_c,FREE_C)
void MUMPS_CALL
MUMPS_FREE_C(void *address);
#define MUMPS_RCOPY_32TO64_64C_IP_C \
    F_SYMBOL(rcopy_32to64_64c_ip_c,RCOPY_32TO64_64C_IP_C)
void MUMPS_CALL
MUMPS_RCOPY_32TO64_64C_IP_C(float *inouttab, MUMPS_INT8 *sizetab);
#define MUMPS_RCOPY_64TO32_64C_IP_C \
    F_SYMBOL(rcopy_64to32_64c_ip_c,RCOPY_64TO32_64C_IP_C)
void MUMPS_CALL
MUMPS_RCOPY_64to32_64C_IP_C(double *inouttab, MUMPS_INT8 *sizetab);
/*
 * Persistent single-worker dispatcher.  MUMPS_THREAD_SAFE_VERSION is 1.
 * mumps_ts_dispatch returns 0 once the worker has run the call and -1 when
 * the worker is unavailable; mumps_is_thread_safe is the capability symbol
 * clients can probe.  caller_omp_max is the OpenMP thread count to apply on
 * the worker for this call (0 leaves the worker's setting alone); it exists
 * so that ICNTL(16) <= 0 honours a caller's omp_set_num_threads().
 * mumps_ts_dispatch_f77 is the same entry point exposed for the C side, and
 * MUMPS_TS_DISPATCH_F77 follows the compiler symbol convention so that
 * [SDCZ]MUMPS can reach it with a plain EXTERNAL declaration from Fortran.
 *
 * MUMPS_TS_DISPATCH_ERROR (-99) is reported in INFO(1)/INFOG(1) when the
 * dispatcher cannot run the call at all.  It is deliberately outside the
 * documented negative INFO(1) values (-1, -8, -10, -13, -14, -60, -80, -81,
 * -800), in particular -9, which means that workarray S is too small.
 */
#define MUMPS_THREAD_SAFE_VERSION 1
#define MUMPS_TS_DISPATCH_ERROR (-99)
typedef void (*mumps_ts_fn)(void *);
int mumps_ts_dispatch(int arith, void *struc, mumps_ts_fn fn,
                      int caller_omp_max);
int MUMPS_CALL mumps_ts_dispatch_f77(int *arith, void *struc, mumps_ts_fn fn,
                                     int *caller_omp_max);
#define MUMPS_TS_DISPATCH_F77 \
    F_SYMBOL(ts_dispatch_f77,TS_DISPATCH_F77)
int MUMPS_TS_DISPATCH_F77(int *arith, void *struc, mumps_ts_fn fn,
                          int *caller_omp_max);
int MUMPS_CALL mumps_is_thread_safe(void);
#endif /* MUMPS_COMMON_H */
