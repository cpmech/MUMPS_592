/*
 * L0-OMP solve probe: a 2-D 5-point Laplacian solved through the C interface
 * with MUMPS' own OpenMP threads (ICNTL(16) = 0, ICNTL(48) = 1 default), so the
 * tree-parallel solve path ([sdcz]sol_omp_m / [sdcz]sol_lr) is exercised.
 * Phase 0 L0-OMP/large-problem regression test: 250x250 five-point Laplacian,
 * 62 500 unknowns, 311 500 nonzeros, ICNTL(16)=0.
 * Prints INFOG(1), the factorization size and the relative residual, so two
 * libraries can be compared on the same run.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "dmumps_c.h"

#define NX 250
#define NY 250
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1]

int main(int argc, char **argv)
{
    DMUMPS_STRUC_C id;
    MUMPS_INT n = (MUMPS_INT) NX * NY;
    MUMPS_INT8 nnz;
    MUMPS_INT *irn, *jcn;
    double *a, *rhs;
    double res, res0, xmax = 0.0, xsum = 0.0;
    int i, j, k, myid, ierr, nt, mb;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    nt = (argc > 1) ? atoi(argv[1]) : 1;
    mb = (argc > 2) ? atoi(argv[2]) : 0;

    irn = (MUMPS_INT *) malloc(sizeof(MUMPS_INT) * (size_t) (5 * n));
    jcn = (MUMPS_INT *) malloc(sizeof(MUMPS_INT) * (size_t) (5 * n));
    a = (double *) malloc(sizeof(double) * (size_t) (5 * n));
    rhs = (double *) malloc(sizeof(double) * (size_t) n);

    k = 0;
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            int row = j * NX + i + 1;   /* 1-based */
            irn[k] = row; jcn[k] = row; a[k] = 4.0; k++;
            if (i > 0)      { irn[k] = row; jcn[k] = row - 1; a[k] = -1.0; k++; }
            if (i < NX - 1) { irn[k] = row; jcn[k] = row + 1; a[k] = -1.0; k++; }
            if (j > 0)      { irn[k] = row; jcn[k] = row - NX; a[k] = -1.0; k++; }
            if (j < NY - 1) { irn[k] = row; jcn[k] = row + NX; a[k] = -1.0; k++; }
        }
    }
    nnz = k;
    for (i = 0; i < n; i++) rhs[i] = 1.0;

    id.comm_fortran = USE_COMM_WORLD;
    id.par = 1; id.sym = 0; id.job = JOB_INIT;
    dmumps_c(&id);
    id.n = n; id.nnz = nnz; id.irn = irn; id.jcn = jcn; id.a = a; id.rhs = rhs;
    id.ICNTL(1) = -1; id.ICNTL(2) = -1; id.ICNTL(3) = -1; id.ICNTL(4) = 0;
    id.ICNTL(16) = nt;      /* MUMPS internal threads: 0 = OpenMP default */
    if (mb > 0) id.ICNTL(23) = mb;   /* working memory in MB: no dynamic CB */

    id.job = 1; dmumps_c(&id);   /* analysis */
    if (id.infog[0] < 0) { printf("ANALYSIS FAILED INFOG(1)=%d\n", (int) id.infog[0]); return 2; }
    id.job = 2; dmumps_c(&id);   /* factorization */
    if (id.infog[0] < 0) { printf("FACTORIZATION FAILED INFOG(1)=%d\n", (int) id.infog[0]); return 3; }
    id.job = 3; dmumps_c(&id);   /* solve */
    if (id.infog[0] < 0) { printf("SOLVE FAILED INFOG(1)=%d\n", (int) id.infog[0]); return 4; }

    /* residual: b - A x, with b = 1 and x = id.rhs (the solve is in place) */
    res = 0.0;
    {
        double *ax = (double *) calloc((size_t) n, sizeof(double));
        for (k = 0; k < nnz; k++) ax[irn[k] - 1] += a[k] * id.rhs[jcn[k] - 1];
        for (i = 0; i < n; i++) {
            double r = fabs(1.0 - ax[i]);
            if (r > res) res = r;
        }
        free(ax);
    }
    for (i = 0; i < n; i++) {
        if (fabs(id.rhs[i]) > xmax) xmax = fabs(id.rhs[i]);
        xsum += id.rhs[i];
    }
    printf("n=%d nnz=%lld nthreads=%d INFOG(1)=%d res_inf=%.6e xmax=%.10f xsum=%.10f\n",
           (int) n, (long long) nnz, nt, (int) id.infog[0], res, xmax, xsum);

    id.job = JOB_END; dmumps_c(&id);
    MPI_Finalize();
    free(irn); free(jcn); free(a); free(rhs);
    return (res < 1e-10) ? 0 : 5;
}
