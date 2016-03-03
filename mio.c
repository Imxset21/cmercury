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
    } else if (y > 0) {
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
