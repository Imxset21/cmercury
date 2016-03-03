#include "mxx.h"

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>

void mxx_sync(
    double time,
    double tstart,
    double h0,
    double tol,
    double jcen[3],
    int nbod,
    int nbig,
    double *m,
    double *x,
    double *v,
    double *s,
    double *rho,
    double *rceh,
    int *stat,
    char *id,
    double *epoch,
    double *ngf,
    int opt[8],
    int ngflag)
{
    fputs("mxx_sync not implemented\n", stderr);
    return;
}

static int compare_doubles(const void *restrict a, const void *restrict b)
{
    double arg1 = *(const double *restrict) a;
    double arg2 = *(const double *restrict) b;

    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

int* mxx_sort(double *restrict arr, const size_t arr_size)
{
    if (arr == NULL)
    {
        fputs("mxx_sort: invalid array size or invalid array\n", stderr);
        return NULL;
    }

    // Create index array that'll hold the positions of the old elements
    int *const index_arr = calloc(arr_size, sizeof(int));
    if (index_arr == NULL)
    {
        fprintf(stderr,
                "mxx_sort: unable to allocate index array of size %lu\n",
                arr_size);
        return NULL;
    }

    // Create copy of array to be sorted so we know the original indices
    double *const copy_of_arr = calloc(arr_size, sizeof(double));
    if (copy_of_arr == NULL)
    {
        fprintf(stderr,
                "mxx_sort: unable to allocate copy array of size %lu\n",
                arr_size);
        free(index_arr);
        return NULL;
    }
    memcpy(copy_of_arr, arr, sizeof(double) * arr_size);

    // O(n log n) ... probably.
    qsort(arr, arr_size, sizeof(double), compare_doubles);

    // Find the new indicies of the elements in the new array
    for (size_t i = 0; i < arr_size; i++)
    {
        // Find the address of the value at i (in our 'old' array)
        const double* val_addr = bsearch(
                                     &copy_of_arr[i],
                                     arr,
                                     arr_size,
                                     sizeof(double),
                                     compare_doubles);
        // Use some pointer arithmetic to calculate its index in the new array
        const ptrdiff_t diff = val_addr - arr;
        index_arr[i] = (int) diff;
    }

    free(copy_of_arr);

    return index_arr;
}

double* mxx_jac(
    double jcen[3],
    const int nbod,
    const int nbig,
    double *m,
    double **xh,
    double **vh)
{
    int itmp[8] = {0}, iflag = 0;
    double x[NMAX][3];
    double v[NMAX][3];
    double temp = 0.0, dx = 0.0, dy = 0.0, dz = 0.0;
    double r = 0.0 , d = 0.0, a2 = 0.0, n = 0.0;
    double tmp2[NMAX][4];

    double *jac = calloc(NMAX, sizeof(double));
    if (jac == NULL)
    {
        fputs("mxx_jac: Unable to allocate jac array\n", stderr);
        return NULL;
    }

    // Convert to barycentric coordinates and velocities
    iflag = mco_h2b(temp, jcen, nbod, nbig, temp, m, xh, vh, x, v, tmp2, &itmp);

    dx = x[1][0] - x[0][0];
    dy = x[1][1] - x[0][1];
    dz = x[1][2] - x[0][2];
    a2 = dx * dx + dy * dy + dz * dz;
    n = sqrt((m[0] + m[1]) / (a2 * sqrt(a2)));

    for(int j = nbig; j <= nbod; j++)
    {
        dx = x[j][0] - x[0][0];
        dy = x[j][1] - x[0][1];
        dz = x[j][2] - x[0][2];
        r = sqrt(dx * dx + dy * dy + dz * dz);
        dx = x[j][0] - x[2][1];
        dy = x[j][1] - x[1][1];
        dz = x[j][2] - x[1][2];
        d = sqrt(dx * dx + dy * dy + dz * dz);
        jac[j] = .50 * (v[j][0] * v[j][0] + v[j][1] * v[j][1] + v[j][2] * v[j][2])
                 - m[0] / r - m[1] / d - n * (x[j][0] * v[j][1] - x[j][1] * v[j][0]);
    }
    return jac;
}

double mxx_en(
    double jcen[3],
    const int nbod,
    const int nbig,
    double *m,
    double **xh,
    double **vh,
    double **s,
    double *e)
{
    int iflag = 0, itmp[8] = {0};
    double x[NMAX][3], v[NMAX][3], l[3] = {0.0};
    double temp = 0.0, dx = 0.0, dy = 0.0, dz = 0.0, r2 = 0.0, ke = 0.0, pe = 0.0;
    double r_1 = 0.0, r_2 = 0.0, r_4 = 0.0, r_6 = 0.0, u2 = 0.0, u4 = 0.0, u6 = 0.0;
    double tmp2[NMAX][4];

    // Local copy of e
    double _e = *e;

    ke = 0.0;
    pe = 0.0;

    // Convert to barycentric coordinates and velocities
    iflag = mco_h2b(temp, jcen, nbod, nbig, temp, m, xh, vh, x, v, tmp2, &itmp);

    // Do the spin angular momenta first (probably the smallest terms)
    for (int j = 0; j < nbod; j++)
    {
        l[0] += s[j][0];
        l[1] += s[j][1];
        l[2] += s[j][2];
    }

    // Orbital angular momentum and kinetic energy terms
    for (int j = 0; j < nbod; j++)
    {
        l[0] += m[j] * (x[j][1] * v[j][2] - x[j][2] * v[j][1]);
        l[1] += m[j] * (x[j][2] * v[j][0] - x[j][0] * v[j][2]);
        l[2] += m[j] * (x[j][0] * v[j][1] - x[j][1] * v[j][0]);
        ke = ke + m[j] * (v[j][0] * v[j][0] + v[j][1] * v[j][1] + v[j][2] * v[j][2]);
    }

    // Potential energy terms due to pairs of bodies
    for (int j = 1; j < nbod; j++)
    {
        double tmp = 0.0;
        for(int k = j + 1; k < nbod; k++)
        {
            dx = x[k][0] - x[j][0];
            dy = x[k][1] - x[j][1];
            dz = x[k][2] - x[j][2];
            r2 = dx * dx + dy * dy + dz * dz;
            if (r2)
            {
                tmp = tmp + m[k] / sqrt(r2);
            }
        }
        pe = pe  -  tmp * m[j];
    }

    // Potential energy terms involving the central body
    for (int j = 1; j < nbod; j++)
    {
        dx = x[j][0] - x[0][0];
        dy = x[j][1] - x[0][1];
        dz = x[j][2] - x[0][2];
        r2 = dx * dx + dy * dy + dz * dz;
        if (r2)
        {
            pe -= m[0] * m[j] / sqrt(r2);
        }
    }

    // Corrections for oblateness
    if (jcen[0] != 0 || jcen[1] != 0 || jcen[2] != 0)
    {
        for (int j = 1; j < nbod; j++)
        {
            r2 = xh[j][0] * xh[j][0] + xh[j][1] * xh[j][1] + xh[j][2] * xh[j][2];
            r_1 = 1.0 / sqrt(r2);
            r_2 = r_1 * r_1;
            r_4 = r_2 * r_2;
            r_6 = r_4 * r_2;
            u2 = xh[j][2] * xh[j][2] * r_2;
            u4 = u2 * u2;
            u6 = u4 * u2;
            pe += m[0] * m[j] * r_1
                  * (jcen[0] * r_2 * (1.5 * u2 - 0.5)
                     +  jcen[1] * r_4 * (4.375 * u4 - 3.75 * u2 + .375)
                     +  jcen[2] * r_6
                     * (14.4375 * u6 - 19.6875 * u4 + 6.5625 * u2 - .3125));
        }
    }

    // Total energy calculation and propagation
    _e = .50 * ke + pe;
    *e = _e;

    return sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
}
