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


#ifndef MUMPS_C_TYPES_H
#define MUMPS_C_TYPES_H

#include <stdint.h>
#include "mumps_compat.h"

/* mumps_int_def.h will define either MUMPS_INTSIZE32 (default)
   or MUMPS_INTSIZE64 (if compilation is with -DINTSIZE64 to
   match Fortran -i8 or equivalent option). This allows one to
   test from an external code whether MUMPS_INT is 64bits or not */
#include "mumps_int_def.h"

#ifdef MUMPS_INTSIZE64
#define MUMPS_INT int64_t
#else
#define MUMPS_INT int
#endif

#define MUMPS_INT8 int64_t

#define SMUMPS_COMPLEX float
#define SMUMPS_REAL float

#define DMUMPS_COMPLEX double
#define DMUMPS_REAL double

/* Complex datatypes */
typedef struct {float r,i;} mumps_complex;
typedef struct {double r,i;} mumps_double_complex;

#define CMUMPS_COMPLEX mumps_complex
#define CMUMPS_REAL float

#define ZMUMPS_COMPLEX mumps_double_complex
#define ZMUMPS_REAL double


#ifndef mumps_ftnlen
/* When passing a string, what is the type of the extra argument
 * passed by value ? */
# define mumps_ftnlen MUMPS_INT
#endif


#define MUMPS_ARITH_s 1
#define MUMPS_ARITH_d 2
#define MUMPS_ARITH_c 4
#define MUMPS_ARITH_z 8

#define MUMPS_ARITH_REAL   ( MUMPS_ARITH_s | MUMPS_ARITH_d )
#define MUMPS_ARITH_CMPLX  ( MUMPS_ARITH_c | MUMPS_ARITH_z )
#define MUMPS_ARITH_SINGLE ( MUMPS_ARITH_s | MUMPS_ARITH_c )
#define MUMPS_ARITH_DBL    ( MUMPS_ARITH_d | MUMPS_ARITH_z )

#define MUMPS_OFF_T MUMPS_INT8

/* Thread-safety capability (additive; the structures and the [sdcz]mumps_c
 * prototypes are unchanged).
 *   Distinct instances may be used concurrently from any threads, but one
 *   MUMPS_STRUC_C must not be accessed concurrently.
 *   1 = calls are dispatched internally to a persistent MUMPS worker.
 *   0 = calls run on the calling thread and the caller must serialize them:
 *       the library was built without pthread support, or for MPI.
 * Consecutive dispatch is active only when MUMPS is built with the
 * sequential MPI stub.  Source-level clients may test
 * MUMPS_THREAD_SAFE_VERSION at compile time and mumps_is_thread_safe() at
 * run time. */
#define MUMPS_THREAD_SAFE_VERSION 1
#ifdef __cplusplus
extern "C" {
#endif
int MUMPS_CALL mumps_is_thread_safe(void);
#ifdef __cplusplus
}
#endif

#endif /* MUMPS_C_TYPES_H */
