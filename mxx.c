#include "mxx.h"

#include "config.h"

#include "mio.h"
#include "mco.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <tgmath.h>

void cmxx_sync(
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
    char **id,
    double *epoch,
    double *ngf,
    int opt[8],
    int ngflag)
{
    fputs("cmxx_sync not implemented\n", stderr);
    return;
}

static int compare_doubles(const void *restrict a, const void *restrict b)
{
    double arg1 = *(const double * restrict) a;
    double arg2 = *(const double * restrict) b;

    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

int* cmxx_sort(double *restrict arr, const size_t arr_size)
{
    if (arr == NULL)
    {
        fputs("cmxx_sort: invalid array size or invalid array\n", stderr);
        return NULL;
    }

    // Create index array that'll hold the positions of the old elements
    int *const index_arr = calloc(arr_size, sizeof(int));
    if (index_arr == NULL)
    {
        fprintf(stderr,
                "cmxx_sort: unable to allocate index array of size %lu\n",
                arr_size);
        return NULL;
    }

    // Create copy of array to be sorted so we know the original indices
    double *const copy_of_arr = calloc(arr_size, sizeof(double));
    if (copy_of_arr == NULL)
    {
        fprintf(stderr,
                "cmxx_sort: unable to allocate copy array of size %lu\n",
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

double* cmxx_jac(
    double jcen[3],
    const int nbod,
    const int nbig,
    double *m,
    const double **restrict xh,
    const double **restrict vh)
{
    int itmp[8] = {0};
    double x[NMAX][3];
    double v[NMAX][3];

    // Make double** compatible alias of x and v for passing to other functions
    double *x_ptr[3];
    for (int i = 0; i < 3; i++)
    {
        x_ptr[i] = x[i];
    }
    double *v_ptr[3];
    for (int i = 0; i < 3; i++)
    {
        v_ptr[i] = v[i];
    }

    double temp = 0.0, dx = 0.0, dy = 0.0, dz = 0.0;
    double r = 0.0 , d = 0.0, a2 = 0.0, n = 0.0;
    double tmp2[NMAX][4];

    double *jac = calloc(NMAX, sizeof(double));
    if (jac == NULL)
    {
        fputs("cmxx_jac: Unable to allocate jac array\n", stderr);
        return NULL;
    }

    // Convert to barycentric coordinates and velocities
    cmco_h2b(nbod, m, xh, vh, x_ptr, v_ptr);

    dx = x[1][0] - x[0][0];
    dy = x[1][1] - x[0][1];
    dz = x[1][2] - x[0][2];
    a2 = dx * dx + dy * dy + dz * dz;
    n = sqrt((m[0] + m[1]) / (a2 * sqrt(a2)));

    for (int j = nbig; j <= nbod; j++)
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

double cmxx_en(
    double jcen[3],
    const int nbod,
    const int nbig,
    double *m,
    const double **restrict xh,
    const double **restrict vh,
    double **s,
    double *e)
{
    int iflag = 0, itmp[8] = {0};
    double l[3] = {0.0};
    double temp = 0.0, dx = 0.0, dy = 0.0, dz = 0.0, r2 = 0.0, ke = 0.0, pe = 0.0;
    double r_1 = 0.0, r_2 = 0.0, r_4 = 0.0, r_6 = 0.0, u2 = 0.0, u4 = 0.0, u6 = 0.0;
    double tmp2[NMAX][4];

    double x[NMAX][3], v[NMAX][3];
    // Make double** compatible alias of x and v for passing to other functions
    double *x_ptr[3];
    for (int i = 0; i < 3; i++)
    {
        x_ptr[i] = x[i];
    }
    double *v_ptr[3];
    for (int i = 0; i < 3; i++)
    {
        v_ptr[i] = v[i];
    }    

    // Local copy of e
    double _e = *e;

    ke = 0.0;
    pe = 0.0;

    // Convert to barycentric coordinates and velocities
    cmco_h2b(nbod, m, xh, vh, x_ptr, v_ptr);

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

void cmxx_elim(
    int *nbod,
    int *nbig,
    double *m,
    double **x,
    double **v,
    double **s,
    double *rho,
    double *rceh,
    double **ngf,
    double *stat,
    char **id,
    char mem[NMESS],
    int lmem[NMESS],
    FILE *outfile)
{
    int elim[NMAX] = {0};

    // Find out how many objects are to be removed
    int nelim = 0;
    int nbigelim = 0;
    for(int j = 1; j < *nbod; j++)
    {
        if (stat[j] < 0)
        {
            elim[nelim] = j;
            nelim++;
            if (j <= *nbig)
            {
                nbigelim++;
            }
        }
    }
    elim[nelim] = *nbod + 1;

    // Eliminate unwanted objects
    for(int k = 0; k < nelim; k++)
    {
        for(int j = elim[k] - k;  j < elim[k + 1] - k; j++)
        {
            const int l = j + k;
            x[j][0] = x[l][0];
            x[j][1] = x[l][1];
            x[j][2] = x[l][2];
            v[j][0] = v[l][0];
            v[j][1] = v[l][1];
            v[j][2] = v[l][2];
            m[j]   = m[l];
            s[j][0] = s[l][0];
            s[j][1] = s[l][1];
            s[j][2] = s[l][2];
            rho[j] = rho[l];
            rceh[j] = rceh[l];
            stat[j] = stat[l];
            strncpy(id[j], id[l], sizeof(char) * 8);
            ngf[j][0] = ngf[l][0];
            ngf[j][1] = ngf[l][1];
            ngf[j][2] = ngf[l][2];
            ngf[j][3] = ngf[l][3];
        }
    }
    // Update total number of bodies and number of Big bodies
    *nbod -= nelim;
    *nbig -= nbigelim;

    // If no massive bodies remain, stop the integration
    if (*nbig < 1)
    {
        puts("Simulation concluded - no massive bodies remain.\n");
        // TODO: Write simulation results to outfile
        // open (23,file=outfile,status='old',access='append',err=10)
        // write (23,'(2a)') mem(81)(1:lmem(81)),mem(124)(1:lmem(124))
        fclose(outfile);
        exit(EXIT_SUCCESS);
    }
}

int cmxx_ejec(
    double time,
    double tstart,
    double rmax,
    double en[3],
    double am[3],
    double jcen[3],
    int i0,
    int nbod,
    int nbig,
    double *m,
    const double **restrict x,
    const double **restrict v,
    double **s,
    int *stat,
    char **id,
    int opt[8],
    FILE *outfile,
    char mem[NMESS],
    char lmem[NMESS])
{
    int year = 0, month = 0;
    double r2 = 0.0, rmax2 = 0.0, t1 = 0.0, e = 0.0, l = 0.0;
    // char flost[38];
    // char tstring[6];

    if (i0 <= 0)
    {
        i0 = 2;
    }
    int ejflag = 0;
    rmax2 = rmax * rmax;

    // Calculate initial energy and angular momentum
    l = cmxx_en(jcen, nbod, nbig, m, x, v, s, &e);

    // Flag each object which is ejected, and set its mass to zero
    for(int j = i0; j < nbod; j++)
    {
        r2 = x[j][0] * x[j][0] + x[j][1] * x[j][1] + x[j][2] * x[j][2];
        if (r2 > rmax2)
        {
            ejflag = 1;
            stat[j] = -3;
            m[j] = 0.0;
            s[j][0] = 0.0;
            s[j][1] = 0.0;
            s[j][2] = 0.0;
        }
        // TODO: Write message to information file
        // open (23,file=outfile,status='old',access='append',err=20)
        if (opt[2] == 1)
        {
            cmio_jd2y(time, &year, &month, &t1);
            // flost = "(1x,a8,a,i10,1x,i2,1x,f8.5)";
            // write (23,flost) id(j),mem(68)(1:lmem(68)),year,month,t1
        } else {
            if (opt[2] == 3)
            {
                t1 = (time - tstart) / 365.250;;
                // tstring = mem[1];
                // flost = "(1x,a8,a,f18.7,a)";
            } else {
                if (opt[2] == 0)
                {
                    t1 = time;
                }
                if (opt[2] == 2)
                {
                    t1 = time - tstart;
                }
                // tstring = mem[0];
                // flost = '(1x,a8,a,f18.5,a)';
            }
            // write (23,flost) id(j),mem(68)(1:lmem(68)),t1,tstring;
        }
        // close(23);
    }

    // If ejections occurred, update ELOST and LLOST
    if (ejflag)
    {
        am[1] = cmxx_en(jcen, nbod, nbig, m, x, v, s, &en[1]);
        en[2] = en[2] + (e - en[1]);
        am[2] = am[2] + (l - am[1]);
    }

    return ejflag;
}
