#include "mio.h"

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int** mio_spl(
    const size_t len,
    const char *str,
    int nsub)
{
    char c = ' ';
    nsub = 0;
    size_t j = 0;

    // Allocate master array
    int **delimit = calloc(100, sizeof(int*));
    if (delimit == NULL)
    {
        fputs("mio_spl: Could not allocate delimiter array.\n", stderr);
        return NULL;
    }

    // Allocate sub-arrays
    for (size_t i = 0; i < 100; i++)
    {
        delimit[i] = calloc(2, sizeof(int));
        if (delimit[i] == NULL)
        {
            fprintf(stderr,
                    "mio_spl: Could not allocate delimiter sub-array[%lu].\n",
                    i);
            return NULL;
        }
    }

ten:
    j++;
    if (j > len)
    {
        goto ninety_nine;
    }
    c = str[j];
    if (c == ' ' || c == '=')
    {
        goto ten;
    }

    // Find the end of string
    size_t k = j;
twenty:
    k++;
    if (k > len)
    {
        goto thirty;
    }
    c = str[k];
    if (c != ' ' && c != '=')
    {
        goto twenty;
    }

    // Store details for this string
thirty:
    nsub = nsub + 1;
    delimit[nsub][0] = j;
    delimit[nsub][1] = k - 1;

    if (k < len)
    {
        j = k;
        goto ten;
    }
ninety_nine:
    return delimit;
}

char* mio_re2c(const double x, const double xmin, const double xmax)
{
    double z = 0.0;

    char *mio_out = calloc(8, sizeof(char));
    if (mio_out == NULL)
    {
        fputs("Unable to allocate mio_out buffer\n", stderr);
        return NULL;
    }
    memset(mio_out, ' ', sizeof(char) * 8);

    const double y = (x - xmin) / (xmax - xmin);

    if (y >= 1)
    {
        for (int j = 0; j < 8; j++)
        {
            mio_out[j] = (char) 255;
        }
    }
    else if (y > 0)
    {
        z = y;
        for(int j = 0; j < 8; j++)
        {
            z = remainder(z, 1.0) * 224.0;
            mio_out[j] = (char) ((int) z + 32);
        }
    }

    return mio_out;
}

int mio_out(
    double time,
    double jcen[3],
    double rcen,
    double rmax,
    int nbod,
    int nbig,
    double *m,
    double **xh,
    double **vh,
    double **s,
    double *rho,
    int *stat,
    char **id,
    int opt[8],
    int opflag,
    int algor,
    FILE *outfile)
{
    int k, len, nchar;
    double rhocgs, k_2, rfac, rcen_2, fr, fv, theta, phi, vtheta, vphi;
    char header[80];
    char c[NMAX][80];
    char mio_fl2c_out[8], mio_re2c_out[8];
    char fout[5] = {'(', 'a', ' ', ' ', ')'};

    rhocgs = AU * AU * AU * K2 / MSUN;
    k_2 = 1.0 / K2;
    rcen_2 = 1.0 / (rcen * rcen);

    // Scaling factor (maximum possible range) for distances
    rfac = log10(rmax / rcen);

    // Create the format list, FOUT, used when outputting the orbital elements
    if (opt[3] == 1) { nchar = 2; }
    if (opt[3] == 2) { nchar = 4; }
    if (opt[3] == 3) { nchar = 7; }
    len = 3 + 6 * nchar;
    // if (len.lt.10) write (fout(3:3),'(i1)') len;
    // if (len.ge.10) write (fout(3:4),'(i2)') len;

    // TODO:  Open the orbital-elements output file
    // open (21, file=outfile, status='old', access='append', err=10)


    // SPECIAL  OUTPUT  PROCEDURE
    // If this is a new integration or a complete output is required (e.g. because
    // the number of objects has changed), then output object details & parameters.
    if (opflag == -1 || opflag == 1)
    {
        // Compose a header line with time, number of objects and relevant parameters
        // TODO: Compose header in mio_out
        /*
        header(1:8)   = mio_fl2(time);
        header(9:16)  = mio_re2(dble(nbig - 1),   0.d0, 11239423.99d0);
        header(12:19) = mio_re2(dble(nbod - nbig),0.d0, 11239423.99d0);
        header(15:22) = mio_fl2(m(1) * k_2);
        header(23:30) = mio_fl2(jcen(1) * rcen_2);
        header(31:38) = mio_fl2(jcen(2) * rcen_2 * rcen_2);
        header(39:46) = mio_fl2(jcen(3) * rcen_2 * rcen_2 * rcen_2);
        header(47:54) = mio_fl2(rcen);
        header(55:62) = mio_fl2(rmax);
        */

        // For each object, compress its index number, name, mass, spin
        // components and density (some of these need to be converted to normal
        // units).
        for(int k = 1; k < nbod; k++)
        {
            // TODO: Compress values and write to file in mio_out
            /*
            c(k)(1:8) = mio_re2(dble(k - 1), 0.d0, 11239423.99d0);
            c(k)(4:11) = id(k);
            c(k)(12:19) = mio_fl2(m(k) * k_2);
            c(k)(20:27) = mio_fl2(s(1,k) * k_2);
            c(k)(28:35) = mio_fl2(s(2,k) * k_2);
            c(k)(36:43) = mio_fl2(s(3,k) * k_2);
            c(k)(44:51) = mio_fl2(rho(k) / rhocgs);
            */
        }
        // TODO: Write compressed output to file
        // write (21,'(a1,a2,i2,a62,i1)') char(12),'6a',algor,header(1:62), opt(4)
        for (int k = 1; k < nbod; k++)
        {
            // write (21,'(a51)') c(k)(1:51);
        }
    }

    // NORMAL  OUTPUT  PROCEDURE

    // TODO: Compose a header line containing the time and number of objects
    /*
    header(1:8)   = mio_fl2(time);
    header(9:16)  = mio_re2(dble(nbig - 1),    0.d0, 11239423.990);
    header(12:19) = mio_re2(dble(nbod - nbig), 0.d0, 11239423.990);
    */

    // TODO: Calculate output variables for each body & convert to compressed format
    for(int k = 1; k < nbod; k++)
    {
        /*
        mco_x2ov (rcen,rmax,m(1),m(k),xh(1,k),xh(2,k),xh(3,k),
                  vh(1,k),vh(2,k),vh(3,k),fr,theta,phi,fv,vtheta,vphi);
        // Calculate object's index number and output variables
        c(k)(1:8) = mio_re2(dble(k - 1), 0.d0, 11239423.99d0);
        c(k)(4:11)                 = mio_re2(fr,     0.d0, rfac);
        c(k)(4+  nchar:11+  nchar) = mio_re2(theta,  0.d0, PI);
        c(k)(4+2*nchar:11+2*nchar) = mio_re2(phi,    0.d0, TWOPI);
        c(k)(4+3*nchar:11+3*nchar) = mio_re2(fv,     0.d0, 1.d0);
        c(k)(4+4*nchar:11+4*nchar) = mio_re2(vtheta, 0.d0, PI);
        c(k)(4+5*nchar:11+5*nchar) = mio_re2(vphi,   0.d0, TWOPI);
        */
    }

    // TODO: Write compressed output to file
    // write (21,'(a1,a2,a14)') char(12),'6b',header(1:14);
    for(int k = 1; k < nbod; k++)
    {
        // write (21,fout) c(k)(1:len);
    }

    fclose(outfile);
    opflag = 0;
    return opflag;
}

int mio_log(
    const double time,
    const double tstart,
    const double en[3],
    const double am[3],
    const int opt[8],
    char **mem,
    int lmem[NMESS])
{
    int year = 0, month = 0;
    char flog[38] = {0};
    char tstring[6] = {0};

    if (opt[2] == 0 || opt[2] == 2)
    {
        strncpy(tstring, mem[0], sizeof(tstring));
        snprintf(flog, sizeof(flog), "(1x,a,f14.1,a,2(a,1p1e12.5))");
    }
    else if (opt[2] == 1)
    {
        snprintf(flog, sizeof(flog), "(1x,a,i10,1x,i2,1x,f4.1,2(a,1p1e12.5))");
    }
    else
    {
        strncpy(tstring, mem[1], sizeof(tstring));
        snprintf(flog, sizeof(flog), "(1x,a,f14.3,a,2(a,1p1e12.5))");
    }

    double tmp0 = 0.0;
    double tmp1 = 0.0;
    if (en[0]) { tmp0 = (en[1] + en[2] - en[0]) / abs(en[0]); }
    if (am[0]) { tmp1 = (am[1] + am[2] - am[0]) / abs(am[0]); }

    double t1 = 0.0;
    if (opt[2] == 1)
    {
        mio_jd2y(time, &year, &month, &t1);
        // write (*,flog) mem(64)(1:lmem(64)), year, month, t1,
        // mem(65)(1:lmem(65)), tmp0,mem(66)(1:lmem(66)), tmp1
    }
    else
    {
        if (opt[2] == 0) { t1 = time; }
        if (opt[2] == 2) { t1 = time - tstart; }
        if (opt[2] == 3) { t1 = (time - tstart) / 365.250; }
        //write (*,flog) mem(63)(1:lmem(63)), t1, tstring,
        //    mem(65)(1:lmem(65)), tmp0, mem(66)(1:lmem(66)), tmp1
    }

    return 0;
}

void mio_jd2y(double jd0, int *year, int *month, double *day)
{
    int i, a, b, c, d, e, g;
    double jd, f, temp, x, y, z;

    if (jd0 <= 0)
    {
        goto fifty;
    }

    jd = jd0 + 0.50;
    i = SIGN((int) fabs(jd), jd);
    f = jd - 1.0 * i;

    // If on or after 15th October 1582
    if (i > 2299160)
    {
        temp = (1.0 * i - 1867216.250) / 36524.250;
        a = SIGN((int) fabs(temp), temp );
        temp = .250 * a;
        b = i + 1 + a - SIGN((int) fabs(temp), temp);
    }
    else
    {
        b = i;
    }

    c = b + 1524;
    temp = (1.0 * c - 122.10) / 365.250;
    d = SIGN((int) fabs(temp), temp );
    temp = 365.250 * d;
    e = SIGN((int) fabs(temp), temp );
    temp = (c - e) / 30.60010;
    g = SIGN((int) fabs(temp), temp );

    temp = 30.60010 * g;
    *day = 1.0 * (c - e) + f - 1.0 * SIGN((int) fabs(temp), temp );

    if (g <= 13) { *month = g - 1; }
    if (g > 13) { *month = g - 13; }

    if (*month > 2) { *year = d - 4716; }
    if (*month <= 2) { *year = d - 4715; }

    if (*day > 32)
    {
        *day -= 32;
        *month += 1;
    }

    if (*month > 12)
    {
        *month -= 12;
        *year += 1;
    }
    return;

fifty:

    // Algorithm for negative Julian day numbers (Duffett-Smith doesn't work)
    x = jd0 - 2232101.5;
    f = x - (int) x;
    if (f < 0) { f = f + 1.0; }
    y = (int) (remainder(x, 1461.0) + 1461.0);
    z = (int) (remainder(y, 365.250));
    *month = (int) ((z + 0.50) / 30.610);
    *day = (int) (z + 1.50 - 30.610 * (double)(*month)) + f;
    *month = remainder(*month + 2, 12) + 1;

    *year = 1399 + (int) (x / 365.250);
    if (x < 0) { year = year - 1; }
    if (*month < 3) { year = year + 1; }
}

char* mio_fl2c(const double x)
{
    char *mio_out = calloc(8, sizeof(char));
    if (mio_out == NULL)
    {
        fputs("mio_fl2c: unable to allocate memory for float conversion.\n", stderr);
        return NULL;
    }
    // Write double to string
    int result = 0;
#ifdef __STDC_LIB_EXT1__
    result = snprintf_s(mio_out, sizeof(char) * 8, "%lf", x);
#else
    result = snprintf(mio_out, sizeof(char) * 8, "%lf", x);
#endif

    if (result < 0)
    {
        fprintf(stderr, "mio_fl2c: error while writing %lf to buffer.\n", x);
        return NULL;
    }

    return mio_out;
}

void mio_err(
    const int unit,
    const char *s1,
    const char *s2,
    const char *s3,
    const char *s4)
{
    fputs("ERROR: Program terminated. See information file for details.\n", stderr);
    fprintf(stderr, "%d: %s\n%s\n%s\n%s", unit, s1, s2, s3, s4);
    exit(EXIT_FAILURE);
}

int mio_ce(
    double time,
    double tstart,
    double rcen,
    double rmax,
    int nbod,
    int nbig,
    double *m,
    int *stat,
    char **id,
    int nclo,
    int *iclo,
    int *jclo,
    int opt[8],
    int stopflag,
    double *tclo,
    double *dclo,
    double **ixvclo,
    double **jxvclo,
    char **mem,
    int lmem[NMESS],
    FILE* outfile,
    int nstored,
    int ceflush)
{
    int k, year, month;
    double tmp0, t1, rfac, fr, fv, theta, phi, vtheta, vphi;
    static char c[200][80];
    char fstop[38];
    char *mio_fl2c_out, *mio_re2c_out;
    char tstring[6];

    // Scaling factor (maximum possible range) for distances
    rfac = log10(rmax / rcen);

    // TODO: Store details of each new close-encounter minimum
    for (int k = 0; k < nclo; k++)
    {
        nstored++;
        /*
        c(nstored)(1:8)   = mio_fl2c(tclo[k]);
        c(nstored)(9:16)  = mio_re2c((double)(iclo[k]-1),0.0,11239423.990);
        c(nstored)(12:19) = mio_re2c((double)(jclo[k]-1),0.0,11239423.990);
        c(nstored)(15:22) = mio_fl2c(dclo[k]);

        mco_x2ov (rcen,rmax,m(1),0.d0,ixvclo(1,k),ixvclo(2,k),
         ixvclo(3,k),ixvclo(4,k),ixvclo(5,k),ixvclo(6,k),fr,theta,phi,
                       fv,vtheta,vphi);
        c(nstored)(23:30) = mio_re2c (fr    , 0.0, rfac);
        c(nstored)(27:34) = mio_re2c (theta , 0.0, PI);
        c(nstored)(31:38) = mio_re2c (phi   , 0.0, TWOPI);
        c(nstored)(35:42) = mio_re2c (fv    , 0.0, 1.0);
        c(nstored)(39:46) = mio_re2c (vtheta, 0.0, PI);
        c(nstored)(43:50) = mio_re2c (vphi  , 0.0, TWOPI);

        mco_x2ov (rcen,rmax,m(1),0.d0,jxvclo(1,k),jxvclo(2,k),
                  jxvclo(3,k),jxvclo(4,k),jxvclo(5,k),jxvclo(6,k),fr,theta,phi,
                  fv,vtheta,vphi);
        c(nstored)(47:54) = mio_re2c (fr    , 0.0, rfac);
        c(nstored)(51:58) = mio_re2c (theta , 0.0, PI);
        c(nstored)(55:62) = mio_re2c (phi   , 0.0, TWOPI);
        c(nstored)(59:66) = mio_re2c (fv    , 0.0, 1.0);
        c(nstored)(63:74) = mio_re2c (vtheta, 0.0, PI);
        c(nstored)(67:78) = mio_re2c (vphi  , 0.0, TWOPI);
        */
    }

    // If required, output the stored close encounter details
    if (nstored >= 100 || ceflush == 0)
    {
        // open (22, file=outfile(2), status='old', access='append',err=10);
        for(int k = 0; k < nstored; k++)
        {
            //write (22,'(a1,a2,a70)') char(12),'6b',c(k)(1:70)
        }
        // close (22)
        nstored = 0;
    }

    // If new encounter minima have occurred, decide whether to stop integration
    stopflag = 0;
    if (opt[0] == 1 && nclo > 0)
    {
        // open (23, file=outfile(3), status='old', access='append',err=20)
        // If time style is Gregorian date then...
        tmp0 = tclo[0];
        if (opt[2] == 1)
        {
            snprintf(fstop, sizeof(fstop), "(5a,/,9x,a,i10,1x,i2,1x,f4.1)");
            mio_jd2y(tmp0, &year, &month, &t1);
            // write (23,fstop) mem(121)(1:lmem(121)),mem(126)
            //      (1:lmem(126)),id(iclo(1)),',',id(jclo(1)),
            //      mem(71)(1:lmem(71)),year,month,t1
        }
        // Otherwise...
        else 
        {
            if (opt[2] == 3)
            {
                // tstring = mem(2);
                snprintf(fstop, sizeof(fstop), "(5a,/,9x,a,f14.3,a)");
                t1 = (tmp0 - tstart) / 365.250;
            }
            else
            {
                // tstring = mem(1);
                snprintf(fstop, sizeof(fstop), "(5a,/,9x,a,f14.1,a)");
                if (opt[2] == 0) { t1 = tmp0; }
                if (opt[2] == 2) { t1 = tmp0 - tstart; }
            }
            // write (23,fstop) mem(121)(1:lmem(121)),mem(126)
            // %      (1:lmem(126)),id(iclo(1)),',',id(jclo(1)),
            // %      mem(71)(1:lmem(71)),t1,tstring
        }
        stopflag = 1;
        //close(23);
    }

    return 0;
}

double mio_c2re(
    char *restrict c,
    const double xmin,
    const double xmax,
    const size_t nchar __attribute__((unused)))
{
    double y = 0;

    int result = EOF;
#ifdef __STDC_LIB_EXT1__
    result = sscanf_s(c, "%lf", &y);
#else
    result = sscanf(c, "%lf", &y);
#endif
    if (result == EOF || result != 1)
    {
        fprintf(stderr, "mio_c2re: Unable to parse double in c: %s", c);
        return 0xDEADBEEF;
    }
    return xmin + y * (xmax - xmin);
}
