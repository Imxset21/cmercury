#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "orbel.h"

double orbel_zget(double q)
{
    int iflag = 0;
    double x, tmp, orbel_zget_out;

    if (q < 0.0)
    {
        iflag = 1;
        q = -q;
    }

    // If q is smaller than 0.001
    if (q < 1.0e-3)
    {
        orbel_zget_out = q * (1.0 - (q * q / 3.0) * (1.0 - q * q));
    } else {
        x = 0.50 * (3.0 * q + sqrt(9.0 * (pow(q, 2)) + 4.0));
        tmp = pow(x, (1.0 / 3.0));
        orbel_zget_out = tmp - 1.0 / tmp;
    }

    if (iflag == 1)
    {
        orbel_zget_out = -orbel_zget_out;
        q = -q;
    }

    return orbel_zget_out;
}

double orbel_flon(double e, double capn)
{
    int iflag = 0, IMAX = 10;
    double orbel_flon_out;
    double a, b, sq, biga, bigb;
    double x, x2;
    double f, fp, dx;
    double diff;
    double a0, a1;
    double b1;
    double a11 = 156.0, a9 = 17160.0, a7 = 1235520.0;
    double a5 = 51891840.0, a3 = 1037836800.0;
    double b11 = 11.0 * a11, b9 = 9.0 * a9, b7 = 7.0 * a7;
    double b5 = 5.0 * a5, b3 = 3.0 * a3;

    // Function to solve "Kepler's eqn" for F (here called x) for given e
    // and CAPN. Only good for smallish CAPN.
    if (capn < 0.0)
    {
        iflag = 1;
        capn = -capn;
    }

    a1 = 6227020800.0 * (1.0 - 1.0 / e);
    a0 = -6227020800.0 * capn / e;
    b1 = a1;

    // Set iflag nonzero if capn < 0., in which case solve for -capn
    // and change the sign of the final answer for F.
    // Begin with a reasonable guess based on solving the cubic for small F

    a = 6.0 * (e - 1.0) / e;
    b = -6.0 * capn / e;
    sq = sqrt(0.25 * b * b + a * a * a / 27.0);
    biga = pow((-0.5 * b + sq), 0.3333333333333333);
    bigb = pow(-(0.5 * b + sq), 0.3333333333333333);
    x = biga + bigb;
    orbel_flon_out = x;

    // If capn is tiny (or zero) no need to go further than cubic even for e =1
    if (capn < TINY)
    {
        goto orbel_flon_normal_ret;
    }

    for (int i = 1; i < IMAX; i++)
    {
        x2 = x * x;
        f = a0 + x * (a1 + x2 * (a3 + x2 * (a5 + x2 * (a7 + x2 * (a9 + x2 * (a11 + x2))))));
        fp = b1 + x2 * (b3 + x2 * (b5 + x2 * (b7 + x2 * (b9 + x2 * (b11 + 13.0 * x2)))));
        dx = -f / fp;

        // TODO: Convert 'format(1x, i3, 3(2x,1p1e22.15))' to something C-like
        orbel_flon_out = x + dx;

        // If we have converged here there's no point in going on
        if (fabs(dx) < TINY)
        {
            goto orbel_flon_normal_ret;
        }

        x = orbel_flon_out;
    }

    // Abnormal return here - we've gone thru the loop IMAX times without convergence
    if (iflag == 1)
    {
        orbel_flon_out = -orbel_flon_out;
        capn = -capn;
    }
    // Print out diagnostic information regarding non-convergence state
    fputs("FLON : RETURNING WITHOUT COMPLETE CONVERGENCE\n", stderr);
    diff = e * sinh(orbel_flon_out) - orbel_flon_out - capn;
    fprintf(stderr,
            "N: %f, F: %f, ecc * sinh(F) - F - N: %f\n",
            capn, orbel_flon_out, diff);
    return orbel_flon_out;

    // Normal return here, but check if capn was originally negative
orbel_flon_normal_ret:
    if (iflag == 1)
    {
        orbel_flon_out = -orbel_flon_out;
        capn = -capn;
    }

    return orbel_flon_out;
}
