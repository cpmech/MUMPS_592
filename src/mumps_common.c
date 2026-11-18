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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
# include <pthread.h>
#endif
#include "mumps_common.h"
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
/*
 * One process-wide lock for the [SDCZ]MUMPS_F77 instance tables. It only
 * guards table scans/growth and instance creation/destruction; solver calls
 * themselves do not hold it.
 */
static pthread_mutex_t mumps_f77_mutex = PTHREAD_MUTEX_INITIALIZER;
void MUMPS_CALL
MUMPS_F77_LOCK(void)
{
    if (pthread_mutex_lock(&mumps_f77_mutex) != 0)
      {
        fprintf(stderr, "MUMPS: fatal error: mumps_f77_lock failed\n");
        abort();
      }
}
void MUMPS_CALL
MUMPS_F77_UNLOCK(void)
{
    if (pthread_mutex_unlock(&mumps_f77_mutex) != 0)
      {
        fprintf(stderr, "MUMPS: fatal error: mumps_f77_unlock failed\n");
        abort();
      }
}
#else
/* Thread safety needs pthreads; without them the instance-table lock is a
 * no-op (same convention as mumps_io.c and mumps_register_thread.c). */
void MUMPS_CALL
MUMPS_F77_LOCK(void)
{
}
void MUMPS_CALL
MUMPS_F77_UNLOCK(void)
{
}
#endif
void MUMPS_CALL
MUMPS_ICOPY_32TO64_64C_IP_C(MUMPS_INT *inouttab, MUMPS_INT8 *sizetab)
/* Copies in-place *sizetab int values starting at address inouttab
   into *sizetab int64_t values starting at the same address.
*/
{
   MUMPS_INT8 i8; /* signed integer needed for reversed loop below */
   for (i8=*sizetab-1; i8 >=0; i8--)
     {
       // See "strict aliasing" comment in MUMPS_COPY_32TO64_64C_IP_C
       /* outtab8[i8]=(MUMPS_INT8)intab4[i8]; */
       //((MUMPS_INT8 *)inouttab)[i8]=(MUMPS_INT8)inouttab[i8];
       MUMPS_INT8 tmpval = (MUMPS_INT8) inouttab[i8];
       memcpy( ((char*) inouttab) + i8 * sizeof(MUMPS_INT8), (char*) &tmpval, sizeof(MUMPS_INT8) ); 
     }
}
void MUMPS_CALL
MUMPS_ICOPY_64TO32_64C_IP_C(MUMPS_INT8 *inouttab, MUMPS_INT8 *sizetab)
/* Copies in-place *sizetab int64_t values starting at address inouttab
   into *sizetab int values starting at the same address */
{
   MUMPS_INT8 i8;
   for (i8=0; i8 < *sizetab; i8++)
     {
       // See "strict aliasing" comment in MUMPS_COPY_32TO64_64C_IP_C
       /*       outtab4[i8]=(MUMPS_INT)intab8[i8]; */
       //((MUMPS_INT *)inouttab)[i8]=(MUMPS_INT)inouttab[i8];
       MUMPS_INT tmpval = (MUMPS_INT) inouttab[i8];
       memcpy( ((char*) inouttab) + i8 * sizeof(MUMPS_INT), (char*) &tmpval, sizeof(MUMPS_INT) ); 
     }
}
void MUMPS_CALL
MUMPS_MALLOC_C(MUMPS_INT8 *address, MUMPS_INT8 *size)
{
  void * ptr;
  ptr=malloc(*size);
  *address=(MUMPS_INT8)(ptr);
}
void MUMPS_CALL
MUMPS_FREE_C(void *address)
{
  free(address);
}
void MUMPS_CALL
MUMPS_RCOPY_32TO64_64C_IP_C(float *inouttab, MUMPS_INT8 *sizetab)
/* Copies in-place *sizetab float values starting at address inouttab
   into *sizetab double values starting at the same address.
*/
{
   MUMPS_INT8 i8; /* signed integer needed for reversed loop below */
   for (i8=*sizetab-1; i8 >=0; i8--)
   {
       /* Initial implementation was undefined in C
	* the "strict aliasing" rule allows the compiler
	* to assume pointer of incompatible types as independant objects.
	* It can lead to undefined behaviors
	*/
       //((double *)inouttab)[i8]=(double)inouttab[i8];
       double tmpval = (double) inouttab[i8];
       memcpy( ((char*) inouttab) + i8 * sizeof(double), (char*) &tmpval, sizeof(double) ); 
   }
}
void MUMPS_CALL
MUMPS_RCOPY_64TO32_64C_IP_C(double *inouttab, MUMPS_INT8 *sizetab)
/* Copies in-place *sizetab double values starting at address inouttab
   into *sizetab float values starting at the same address */
{
   MUMPS_INT8 i8;
   for (i8=0; i8 < *sizetab; i8++)
     {
       // See "strict aliasing" comment in MUMPS_COPY_32TO64_64C_IP_C
       /*       outtab4[i8]=(MUMPS_INT)intab8[i8]; */
       //((float *)inouttab)[i8]=(float)inouttab[i8];
       float tmpval = (float) inouttab[i8];
       memcpy( ((char*) inouttab) + i8 * sizeof(float), (char*) &tmpval, sizeof(float) ); 
     }
}
