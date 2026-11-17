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
 *  Modification for thread safety (Phase 0), 2026-11-17
 *  Author: Dorival Pedroso (Kleo completed the edit)
 *  Nature: removed the getter/setter declarations for the six process-global
 *          output pointers that [SDCZ]MUMPS_F77 no longer needs, and declared
 *          the process-wide lock protecting the F77 instance table.
 *  ===========================================================================
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
 * Persistent single-worker dispatcher (Phase 1).  MUMPS_THREAD_SAFE_VERSION
 * is 1.  mumps_ts_dispatch returns 0 once the worker has run the call and
 * -1 when the worker is unavailable; mumps_is_thread_safe is the capability
 * symbol clients can probe.  mumps_ts_dispatch_f77 is the same entry point
 * exposed for the C side, and MUMPS_TS_DISPATCH_F77 follows the compiler
 * symbol convention so that [SDCZ]MUMPS can reach it with a plain EXTERNAL
 * declaration from Fortran.
 */
#define MUMPS_THREAD_SAFE_VERSION 1
#define MUMPS_TS_DISPATCH_ERROR (-9)
typedef void (*mumps_ts_fn)(void *);
int mumps_ts_dispatch(int arith, void *struc, mumps_ts_fn fn);
int mumps_ts_dispatch_f77(int arith, void *struc, mumps_ts_fn fn);
#define MUMPS_TS_DISPATCH_F77 \
    F_SYMBOL(ts_dispatch_f77,TS_DISPATCH_F77)
int MUMPS_TS_DISPATCH_F77(int arith, void *struc, mumps_ts_fn fn);
int mumps_is_thread_safe(void);
#endif /* MUMPS_COMMON_H */
