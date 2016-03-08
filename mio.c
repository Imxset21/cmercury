#include "config.h"

#include "mio.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <tgmath.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>

//! Simulation algorithm choice
enum SIMUL_ALGO
{
    BS_ALGO = 0,
    MVS_ALGO,
    BS2_ALGO,
    HYBRID_ALGO,
    RADAU_ALGO
};

//! Simulation precision (not to be confused with floating/double precision!)
enum SIMUL_PRECISION
{
    LOW_PRECISION = 0,
    MEDIUM_PRECISION,
    HIGH_PRECISION,
};

//! Convenience object for passing simulation parameters
typedef struct simul_params
{
    //! Simulation algorithm
    enum SIMUL_ALGO algo;
    //! Start time (days)
    double simul_start;
    //! Stop time (days)
    double stop_time;
    //! Output interval (days)
    double output_interval;
    //! Timestep (days)
    unsigned int timestep;
    //! Accuracy parameter 
    double accuracy_param;
    //! Stop integration after a close encounter
    bool encounter_stop;
    //! Allow collisions to occur
    bool allow_collisions;
    //! Include collisional fragmentation
    bool collisional_fragmentation;
    //! Express time in days or years
    bool time_in_years;
    //! Express time relative to integration start time
    bool time_relative_to_start;
    //! Output precision
    enum SIMUL_PRECISION precision;
    //! Include user-defined force
    bool use_user_force;
    //! Ejection distance (AU)
    double ejection_distance_AU;
    //! Radius of central body (AU)
    double central_body_radius_AU;
    //! Central mass (solar)
    double central_mass_solar;
    double central_j2;
    double central_j4;
    double central_j6;
} *simul_params_t;

/**
 * @brief Function pointer to filtering function for parameter acceptance
 * @param[in] val_ptr Pointer to value to test
 * @returns true if the value is accepted, false otherwise
 */
typedef bool (*accept_val_t)(void* val_ptr);

char* cmio_re2c(const double x, const double xmin, const double xmax)
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

char* cmio_fl2c(const double x __attribute__((unused)))
{
    fputs("cmio_fl2c not implemented.\n", stderr);
    return NULL;
}

int cmio_out(
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
    int len, nchar;
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

int cmio_log(
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
    if (en[0]) { tmp0 = (en[1] + en[2] - en[0]) / fabs(en[0]); }
    if (am[0]) { tmp1 = (am[1] + am[2] - am[0]) / fabs(am[0]); }

    double t1 = 0.0;
    if (opt[2] == 1)
    {
        cmio_jd2y(time, &year, &month, &t1);
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

void cmio_jd2y(double jd0, int *year, int *month, double *day)
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

void cmio_err(void)
{
    fputs("ERROR: Program terminated. See information file for details.\n", stderr);
    exit(EXIT_FAILURE);
}

int cmio_ce(
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
    int year, month;
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
            cmio_jd2y(tmp0, &year, &month, &t1);
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

double cmio_c2re(
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

/**
 * @brief Reads the messages.in file and writes its contents to the given vars
 * @param[out] lines NMESS-sized array of 80-character arrays for file lines
 * @returns Number of lines read, or -1 if reading the file failed
 */
static int cmio_in_read_messages_in(char (*lines)[NMESS][80])
{
    if (lines == NULL)
    {
        fputs("cmio_read_messages_in: lines array is NULL\n", stderr);
        return -1;
    }
    
    // Open file for reading
    const char *min_fname = "message.in";
    FILE* message_in = fopen(min_fname, "r");
    if (message_in == NULL)
    {
        fprintf(stderr,
                "cmio_read_messages_in: unable to open %s: %s\n", min_fname,
                (errno == ENOENT) ?
                "file doesn't exist." :
                "an error occurred while reading the file.");
        return -1;
    }

    // Read lines one by one into the provided buffers
    int linum = 0;
    while (fgets(lines[0][linum], sizeof(lines[0][linum]), message_in) != NULL)
    {
        linum++;
    }

    // Check if we reached the end of the file or some other error occurred
    if (!feof(message_in) && (ferror(message_in)))
    {
        perror("fgets()");
        fprintf(stderr,"fgets() failed in at line #%d\n", linum);
        linum = -1;
    }
    fclose(message_in);

    return linum;
}

/**
 * @brief Gets the address of the '=' delimited substring in a string
 * @param[in] str String to search for substring in
 * @param[out] sub_str_len Length of substring
 * @returns Pointer to start of the substring, or NULL if an error occurred
 */
static char* get_postfix_substring(const char* str, size_t *sub_str_len)
{
    if (str == NULL)
    {
        fputs("get_postfix_substring: string is null\n", stderr);
        return NULL;
    } else if (sub_str_len == NULL) {
        fputs("get_postfix_substring: sub_str_len is null\n", stderr);
        return NULL;
    }
    
    const size_t str_len = strlen(str);
    size_t seen_equals_index = 0;
    bool seen_equals = false;
    // Scan the string
    for (size_t i = 0; i < str_len; i++)
    {
        if (str[i] == '=')
        {
            seen_equals_index = i;
            seen_equals = true;
            continue;
        }

        if (seen_equals)
        {
            *sub_str_len += 1;
        }
    }

    if (!seen_equals)
    {
        fputs("get_postfix_substring: couldn't locate substring start\n", stderr);
        return NULL;
    }

    char *substring = calloc(*sub_str_len, sizeof(char));
    if (substring == NULL)
    {
        fputs("get_postfix_substring: couldn't allocate substring\n", stderr);
        return NULL;
    }

    // Copy substring, ignoring non-alphanumeric characters (but including .)
    size_t i = 0;
    for (size_t j = seen_equals_index; j < str_len; j++)
    {
        if (isalnum(str[j]) || str[j] == '.' || str[j] == '+' || str[j] == '-')
        {
            substring[i] = str[j];
            i++;
        }
    }

    return substring;
}

/**
 * @brief Extracts a double floating-point value from a parameter file line
 * @param[in] buff String to slice and extract value from
 * @param[in] accept_fun Function to test if value is to be accepted
 * @param[out] val Pointer to value to set
 * @returns 0 is successful, -1 otherwise
 */
static int extract_double_param(
    const char* buff,
    accept_val_t accept_fun,
    double* val)
{
    if (buff == NULL)
    {
        fputs("extract_double_param: buff cannot be null\n", stderr);
        return -1;
    } else if (val == NULL) {
        fputs("extract_double_param: val cannot be null\n", stderr);
        return -1;
    }
    
    // Attempt to parse start time string
    size_t str_len = 0;
    char* sub_str = get_postfix_substring(buff, &str_len);
    if (sub_str == NULL)
    {
        fputs("read_params: failed to parse stime string\n", stderr);
        return -1;
    }
    // Convert value to double
    double temp = 0.0;
    if (sscanf(sub_str, "%lf", &temp) != 1)
    {
        fputs("read_params: failed to convert stime value\n", stderr);
        free(sub_str);
        return -1;
    } else if (accept_fun) {
        // If accept function is non-NULL, run it and test parameter value
        if (!accept_fun(&temp))
        {
            fprintf(stderr,
                    "read_params: illegal value for parameter: %lf\n",
                    temp);
            free(sub_str);
            return -1;
        }
    }
    *val = temp;
    free(sub_str);
    return 0;
}

/**
 * @brief Extracts an integer value from a parameter file line
 * @param[in] buff String to slice and extract value from
 * @param[in] accept_fun Function to test if value is to be accepted
 * @param[out] val Pointer to value to set
 * @returns 0 is successful, -1 otherwise
 */
static int extract_int_param(
    const char* buff,
    accept_val_t accept_fun,
    int* val)
{
    if (buff == NULL)
    {
        fputs("extract_int_param: buff cannot be null\n", stderr);
        return -1;
    } else if (val == NULL) {
        fputs("extract_int_param: val cannot be null\n", stderr);
        return -1;
    }
    
    // Attempt to parse start time string
    size_t str_len = 0;
    char* sub_str = get_postfix_substring(buff, &str_len);
    if (sub_str == NULL)
    {
        fputs("extract_int_param: failed to parse stime string\n", stderr);
        return -1;
    }
    // Convert value to int
    int temp = 0.0;
    if (sscanf(sub_str, "%d", &temp) != 1)
    {
        fputs("extract_int_param: failed to convert stime value\n", stderr);
        free(sub_str);
        return -1;
    }  else if (accept_fun) {
        // If accept function is non-NULL, run it and test parameter value
        if (!accept_fun(&temp))
        {
            fprintf(stderr,
                    "extract_int_param: illegal value for parameter: %i\n",
                    temp);
            free(sub_str);
            return -1;
        }
    }
    *val = temp;
    free(sub_str);
    return 0;
}

/**
 * @brief Extracts a yes/no value from a parameter file line as a boolean
 * @param[in] buff String to slice and extract value from
 * @param[out] val Pointer to value to set
 * @returns 0 is successful, -1 otherwise
 */
static int extract_yesno_param(const char* buff, bool* val)
{
    if (buff == NULL)
    {
        fputs("extract_yesno_param: buff cannot be null\n", stderr);
        return -1;
    } else if (val == NULL) {
        fputs("extract_yesno_param: val cannot be null\n", stderr);
        return -1;
    }
    
    // Attempt to parse start time string
    size_t str_len = 0;
    char* sub_str = get_postfix_substring(buff, &str_len);
    if (sub_str == NULL)
    {
        fputs("extract_yesno_param: failed to parse stime string\n", stderr);
        return -1;
    }

    // Canary for if neither "yes" or "no" are found ;)
    int tristate = -1;
    // Detect if "yes" or "no" is present in the string
    if (strstr(sub_str, "yes") || strstr(sub_str, "YES"))
    {
        tristate = 1;
    } else if (strstr(sub_str, "no") || strstr(sub_str, "NO")) {
        tristate = 0;
    }
    
    if (tristate == -1)
    {
        fputs("extract_yesno_param: ", stderr);
        return -1;
    }

    *val = (tristate) ? true : false;
    free(sub_str);
    return 0;
}

static enum SIMUL_ALGO algo_str_to_enum(const char* buff)
{
    // Attempt to parse algorithm string
    size_t algo_str_len = 0;
    const char* str = get_postfix_substring(buff, &algo_str_len);
    if (str == NULL)
    {
        fputs("read_params: failed to parse algo type\n", stderr);
        return -1;
    }

    if (strncmp("mvs", str, sizeof("mvs")) == 0)
    {
        return MVS_ALGO;
    } else if (strncmp("bs2", str, sizeof("bs2")) == 0) {
        return BS2_ALGO;
    } else if (strncmp("hybrid", str, sizeof("hybrid")) == 0) {
        return HYBRID_ALGO;
    } else if (strncmp("radau", str, sizeof("radau")) == 0) {
        return RADAU_ALGO;
    } else if (strncmp("bs", str, sizeof("bs")) == 0) {
        return BS_ALGO;
    } else {
        return -1;
    }
}

//! Function that only accepts a double value if it is nonnegative and finite
static bool accept_nonnegfin_double(void* val)
{
    const double *dval = (const double*) val;
    return (val && isfinite(*dval) && *dval >= 0.00);
}

//! Function that only accepts an int value if it is positive and finite
static bool accept_posfin_int(void* val)
{
    const int *ival = (const int*) val;
    return (val && isfinite(*ival) && *ival > 0);
}

//! Function that only accepts a double value if it is positive and finite
static bool accept_posfin_double(void* val)
{
    const double *dval = (const double*) val;
    return (val && isfinite(*dval) && *dval > 0);
}

/**
 * @brief Reads files.in file and places result in input params
 * @param[out] infiles Input filenames
 * @param[out] outfiles Output filenames
 * @param[out] dumpfiles Dump filenames
 * @returns 0 if successful, or -1 if reading the file failed
 */
static int cmio_in_read_files_in(
    char (*infiles)[3][80],
    char (*outfiles)[3][80],
    char (*dumpfiles)[4][80])
{
    if (infiles == NULL || outfiles == NULL || dumpfiles == NULL)
    {
        fputs("cmio_read_messages_in: one of the arrays is NULL\n", stderr);
        return -1;
    }

    // Open the file
    const char *fin_fname = "files.in";
    FILE* files_in = fopen(fin_fname, "r");
    if (files_in == NULL)
    {
        fprintf(stderr,
                "cmio_read_files_in: unable to open %s: %s\n", fin_fname,
                (errno == ENOENT) ?
                "file doesn't exist." :
                "an error occurred while reading the file.");
        return -1;
    }

    // Read the following filenames:
    // *3* input files (*.in)
    // *3* output files (*.out)
    // *4* dump files (*.dmp)
    int success = 0;
    for (size_t i = 0; i < 3 + 3 + 4; i++)
    {
        // Note: have to do (i % j) because of each array's indexing
        void* r = NULL;
        if (i < 3)
        {
            r = fgets(infiles[0][i], sizeof(infiles[0][i]), files_in);
        } else if (i < 6) {
            r = fgets(outfiles[0][i % 3], sizeof(outfiles[0][i % 3]), files_in);
        } else {
            r = fgets(dumpfiles[0][i % 4], sizeof(dumpfiles[0][i % 4]), files_in);
        }

        // If reading the file failed for any reason
        if (r == NULL)
        {
            perror("cmio_read_files_in: fgets()");
            fprintf(stderr, "cmio_read_files_in: fgets() failed at line #%lu\n", i);
            success = -1;
            break;
        } else {
            // Hackery to 'remove' leading whitespace and trailing newline
            int result = 0;
            if (i < 3)
            {
                result = sscanf(infiles[0][i], " %s\n", infiles[0][i]);
            } else if (i < 6) {
                result = sscanf(outfiles[0][i % 3], " %s\n", outfiles[0][i % 3]);
            } else {
                result = sscanf(dumpfiles[0][i % 4], " %s\n", dumpfiles[0][i % 4]);
            }
            // Check if our sscanf trick didn't work... oops
            if (result != 1)
            {
                fprintf(stderr,
                        "cmio_read_files_in: sscanf failed for line #%lu\n",
                        i);
                success = -1;
                break;
            }
        }
    }
    fclose(files_in);

    return success;
}

/**
 * @brief Attempts to opens the info file, given status of restart file
 * @param[in] info_filename Filename of the information file
 * @param[in] restart_filename Filename of the restart file
 * @returns FILE pointer to info file, NULL if an error occurred
 */
static FILE* cmio_in_open_info(const char (*info_filename)[80],
                               const char (*restart_filename)[80])
{
    // Try to open both the restart file or the info file
    FILE *restart_fp = fopen(restart_filename[0], "r");
    FILE *info_fp = fopen(info_filename[0], "r");

    // If BOTH files exist OR if NEITHER file exists, proceed normally
    // If one file exists but the other doesn't, then we're in trouble
    if ((restart_fp == NULL && info_fp) || (restart_fp && info_fp == NULL))
    {
        if (restart_fp)
        {
            fputs("cio_min: Restart file exists, but info doesn't!\n", stderr);
            fclose(restart_fp);
        }
        if (info_fp)
        {
            fputs("cio_min: Info file exists, but restart doesn't!\n", stderr);
            fclose(info_fp);
        }
        return NULL;
    }

    // Make sure restart file is closed if it exists
    if (restart_fp)
    {
        fclose(restart_fp);
        restart_fp = NULL;
    }
    // Re-open info file in append mode after closing its read-mode cousin
    if (info_fp)
    {
        fclose(info_fp);
    }
    info_fp = fopen(info_filename[0], "a");

    return info_fp;
}

/**
 * @brief Read parameter file values and updates the simulation data object
 * @param[in] param_fname File name of the parameter file
 * @param[in,out] simul_data Simulation data object to be updated
 * @returns 0 if successful, -1 otherwise
 */
int cmio_read_params(const char (*param_fname)[80],
                     cmercury_simul_data_t simul_data)
{
    if (simul_data == NULL)
    {
        fputs("cmio_read_params: invalid simulation data object", stderr);
        return -1;
    }
    if (param_fname == NULL)
    {
        fputs("cmio_read_params: invalid information file name", stderr);
        return -1;
    }

    // Try opening file
    FILE *param_file = fopen(param_fname[0], "r");
    if (param_file == NULL)
    {
        fprintf(stderr, "cmio_read_params: unable to read file %s\n", param_fname[0]);
        return -1;
    }

    char buff[150] = {0};

    // Read lines one by one into the provided buffers
    int linum = 0;
    while (fgets(buff, sizeof(buff), param_file) != NULL)
    {
        // Ignore comment lines
        if (buff[0] == ')')
        {
            linum++;
            continue;
        }

        // TODO:
        
        linum++;
    }

    // Check if we reached the end of the file or some other error occurred
    if (!feof(param_file) && (ferror(param_file)))
    {
        perror("fgets()");
        fprintf(stderr,"fgets() failed in at line #%d\n", linum);
        linum = -1;
    }
    fclose(param_file);

    return linum;
}

/**
 * @brief Reads parameter file and returns a pointer to a parameter object
 * @param[in] fname Name of the parameter file
 * @returns Pointer to parameter object, or NULL if parameter parsing failed
 */
static simul_params_t read_params(const char* fname)
{
    if (fname == NULL)
    {
        fputs("read_params: invalid file name\n", stderr);
        return NULL;        
    }
    FILE* simul_fp = fopen(fname, "r");
    if (simul_fp == NULL)
    {
        fprintf(stderr, "read_params: unable to read file '%s': ", fname);
        perror("");
        return NULL;
    }
    
    simul_params_t params = calloc(1, sizeof(struct simul_params));
    if (params == NULL)
    {
        fputs("read_params: unable to allocate parameter object\n", stderr);
        return NULL;
    }

    char buff[150] = {0};
    size_t linum = 0, absolute_linum = 1;
    while (fgets(buff, sizeof(buff), simul_fp) != NULL)
    {
        // Skip 'comment' lines without incrementing linum
        if (buff[0] == ')' || (buff[0] == ' ' && buff[1] == '<'))
        {
            memset(buff, 0, sizeof(buff));
            absolute_linum++;
            continue;
        }

        // Relative number determines parameter set/changed, not text
        switch (linum)
        {
        // Algorithm choice (encoded as a string)
        case 0:
        {
            // Attempt to convert algorithm string into enum denoting type
            int algo_i = -1;
            if ((algo_i = algo_str_to_enum(buff)) < 0)
            {
                fputs("read_params: invalid algorithm choice\n", stderr);
                goto fail_read_params;
            }
            params->algo = (enum SIMUL_ALGO) algo_i;
            printf("read_params: algo is of type %i\n", algo_i);
            break;
        }
        // Start time (in days)
        case 1:
        {
            double stime = 0.0;
            if (extract_double_param(buff, accept_nonnegfin_double, &stime))
            {
                fputs("read_params: failed to read start time value\n", stderr);
                goto fail_read_params;
            }
            params->simul_start = stime;
            printf("read_params: start_time is %lf\n", stime);
            break;
        }
        // Stop time (in days)
        case 2:
        {
            double stime = 0.0;
            if (extract_double_param(buff, accept_nonnegfin_double, &stime))
            {
                fputs("read_params: failed to read stop time value\n", stderr);
                goto fail_read_params;
            }
            if (stime < params->simul_start)
            {
                fprintf(stderr,
                        "read_params: stop time (%lf) is in the past from start time (%lf)\n",
                        stime, params->simul_start);
                goto fail_read_params;
            }
            params->stop_time = stime;
            printf("read_params: stop_time is %lf\n", stime);
            break;
        }
        // Output interval (days)
        case 3:
        {
            double out_interval = 0.0;
            if (extract_double_param(buff, accept_nonnegfin_double, &out_interval))
            {
                fputs("read_params: failed to read stop time value\n", stderr);
                goto fail_read_params;
            }
            params->output_interval = out_interval;
            printf("read_params: output_interval is %lf\n", out_interval);
            break;
        }
        // Time step (days)
        case 4:
        {
            int timestep = 0;
            if (extract_int_param(buff, accept_posfin_int, &timestep))
            {
                fputs("read_params: failed to read stop time value\n", stderr);
                goto fail_read_params;
            }
            params->timestep = timestep; 
            printf("read_params: output_interval is %d\n", timestep);
            break;
        }
        // Accuracy Parameter
        case 5:
        {
            double acc = 0.0;
            if (extract_double_param(buff, accept_posfin_double, &acc))
            {
                fputs("read_params: failed to read accuracy parameter\n", stderr);
                goto fail_read_params;
            }
            params->accuracy_param = acc;
            printf("read_params: accuracy parameter is %e\n", acc);
            break;
        }
        // Option to stop integration after a close encounter
        case 6:
        {
            bool stop_simul = false;
            if (extract_yesno_param(buff, &stop_simul))
            {
                fputs("read_params: failed to read integration option\n", stderr);
                goto fail_read_params;
            }
            params->encounter_stop = stop_simul;
            printf("read_params: %s? %s\n",
                   "stop integration after a close encounter",
                   stop_simul ? "yes" : "no");
            break;
        }
        // Option to allow collisions to occur
        case 7:
        {
            bool allow_coll = false;
            if (extract_yesno_param(buff, &allow_coll))
            {
                fputs("read_params: failed to read integration option\n", stderr);
                goto fail_read_params;
            }
            params->allow_collisions = allow_coll;
            printf("read_params: %s? %s\n",
                   "allow collisions to occur",
                   allow_coll ? "yes" : "no");
            break;
        }
        // Option to include collisional fragmentation
        case 8:
        {
            bool coll_frag = false;
            if (extract_yesno_param(buff, &coll_frag))
            {
                fputs("read_params: failed to read integration option\n", stderr);
                goto fail_read_params;
            }
            params->collisional_fragmentation = coll_frag;
            printf("read_params: %s? %s\n",
                   "include collisional fragmentation",
                   coll_frag ? "yes" : "no");
            break;
        }
        // Option to express time in days or years
        case 9:
        {
            // Attempt to parse time format string
            size_t str_len = 0;
            char* sub_str = get_postfix_substring(buff, &str_len);
            if (sub_str == NULL)
            {
                fputs("read_params: failed to parse time format\n", stderr);
                goto fail_read_params;
            }
            // Detect if "years" or "days" is present in the string
            if (strstr(sub_str, "years"))
            {
                params->time_in_years = true;
            } else if (strstr(sub_str, "days")) {
                params->time_in_years = false;
            } else {
                fprintf(stderr,
                        "read_params: invalid value for time format: %s\n",
                        sub_str);
                goto fail_read_params;
            }
            printf("read_params: %s? %s\n",
                   "express time in days or years",
                   params->time_in_years ? "years" : "days");
            break;
        }
        // Option to express time relative to integration start time
        case 10:
        {
            bool relative = false;
            if (extract_yesno_param(buff, &relative))
            {
                fputs("read_params: failed to read integration option\n", stderr);
                goto fail_read_params;
            }
            params->time_relative_to_start = relative;
            printf("read_params: %s? %s\n",
                   "express time relative to integration start time",
                   relative ? "yes" : "no");
            break;
        }
        // Option for output precision
        case 11:
        {
            // Attempt to parse output precision string
            size_t str_len = 0;
            char* sub_str = get_postfix_substring(buff, &str_len);
            if (sub_str == NULL)
            {
                fputs("read_params: failed to parse output precision\n", stderr);
                goto fail_read_params;
            }
            // Determine precision
            enum SIMUL_PRECISION precision;
            if (strstr(sub_str, "low"))
            {
                precision = LOW_PRECISION;
            } else if (strstr(sub_str, "medium")) {
                precision = MEDIUM_PRECISION;
            } else if (strstr(sub_str, "high")) {
                precision = HIGH_PRECISION;
            } else {
                fprintf(stderr,
                        "read_params: invalid value for output precision: %s\n",
                        sub_str);
                goto fail_read_params;
            }
            params->precision = precision;
            printf("read_params: precision? %d\n", precision);
            break;
        }
        // Option to include relativity in integration
        case 12:
        {
            printf("read_params: '%s' not yet implemented, ignored\n",
                   "Include relativity in integration");
            break;
        }
        // Option to include user-defined force
        case 13:
        {
            bool user_forces = false;
            if (extract_yesno_param(buff, &user_forces))
            {
                fputs("read_params: failed to read integration option\n", stderr);
                goto fail_read_params;
            }
            params->use_user_force = user_forces;
            printf("read_params: %s? %s\n",
                   "express time relative to integration start time",
                   user_forces ? "yes" : "no");
            break;
        }
        // Ejection distance
        case 14:
        {
            double dist = 0.0;
            if (extract_double_param(buff, accept_posfin_double, &dist))
            {
                fputs("read_params: failed to read ejection distance\n", stderr);
                goto fail_read_params;
            }
            params->ejection_distance_AU = dist;
            printf("read_params: ejection distance is %e\n", dist);
            break;
        }
        // Central body radius
        case 15:
        {
            double radius = 0.0;
            if (extract_double_param(buff, accept_posfin_double, &radius))
            {
                fputs("read_params: failed to read central body radius\n", stderr);
                goto fail_read_params;
            }
            params->central_body_radius_AU = radius;
            printf("read_params: central body radius is %e\n", radius);
            break;
        }
        // Central body radius
        case 16:
        {
            double c_mass = 0.0;
            if (extract_double_param(buff, accept_posfin_double, &c_mass))
            {
                fputs("read_params: failed to read central body mass\n", stderr);
                goto fail_read_params;
            }
            params->central_mass_solar = c_mass;
            printf("read_params: central body mass is %e\n", c_mass);
            break;
        }
        // Central body J2
        case 17:
        {
            double j2 = 0.0;
            if (extract_double_param(buff, accept_nonnegfin_double, &j2))
            {
                fputs("read_params: failed to read central J2\n", stderr);
                goto fail_read_params;
            }
            params->central_j2 = j2;
            printf("read_params: central body J2 is %e\n", j2);
            break;
        }
        // Central body J4
        case 18:
        {
            double j4 = 0.0;
            if (extract_double_param(buff, accept_nonnegfin_double, &j4))
            {
                fputs("read_params: failed to read central J4\n", stderr);
                goto fail_read_params;
            }
            params->central_j4 = j4;
            printf("read_params: central body J4 is %e\n", j4);
            break;
        }
        // Central body J6
        case 19:
        {
            double j6 = 0.0;
            if (extract_double_param(buff, accept_nonnegfin_double, &j6))
            {
                fputs("read_params: failed to read central J6\n", stderr);
                goto fail_read_params;
            }
            params->central_j6 = j6;
            printf("read_params: central body J6 is %e\n", j6);
            break;
        }
        // Option to set hybrid integrator changeover
        case 20:
        {
            printf("read_params: '%s' not yet implemented, ignored\n",
                   "Hybrid integrator changeover");
            break;
        }
        // Option to set number of timesteps between data dumps
        case 21:
        {
            printf("read_params: '%s' not yet implemented, ignored\n",
                   "Set number of timesteps between data dumps");
            break;
        }
        // Option to set number of timesteps between data dumps
        case 22:
        {
            printf("read_params: '%s' not yet implemented, ignored\n",
                   "Set number of timesteps between periodic effects");
            break;
        }
        default:
        {
            fprintf(stderr,
                    "read_params: unrecognized parameter, skipping line %lu\n",
                    absolute_linum);
            break;
        }
        }

        absolute_linum++;
        linum++;
    }

    // Test reason for reaching NULL
    if (!feof(simul_fp) && ferror(simul_fp))
    {
        perror("read_params: fgets()");
fail_read_params:
        fprintf(stderr,
                "read_params: failed reading param file at line %lu\n",
                absolute_linum);
        free(params);
        fclose(simul_fp);
        return NULL;
    }

    fclose(simul_fp);
    return params;
}

int cmio_in(cmercury_simul_data_t simul_data)
{
    if (simul_data == NULL)
    {
        fputs("cmio_in: unable to read simulation data object\n", stderr);
        return -1;
    }

    int itmp,jtmp,informat,lim[10][2],nsub,year,month,lineno;
    double q,a,e,i,p,n,l,temp,tmp2,tmp3,rhocgs,t1,tmp4,tmp5,tmp6;
    bool test,oldflag,flag1,flag2;
    char c1;
    char c3[3] = {0}, alg[60][3];
    char filename[80] = {0}, c80[80] = {0};
    char string[150] = {0};
    
    // FIXME: I have no clue with this is decided at runtime and not as a macro
    rhocgs = AU * AU * AU * K2 / MSUN;

    // Read messages in
    char mem[NMESS][80];
    if (cmio_in_read_messages_in(&mem) < 0)
    {
        fputs("cmio_in: unable to read messages file.\n", stderr);
        return -1;
    }

    // Read filenames in
    char outfiles[3][80], infiles[3][80], dumpfiles[4][80];
    if (cmio_in_read_files_in(&infiles, &outfiles, &dumpfiles))
    {
        fputs("cmio_in: unable to read files file.\n", stderr);
        return -1;
    }

    // Try to open the info file, if it exists, or create it if it doesn't
    FILE *info_fp = cmio_in_open_info((const char (*)[80]) &outfiles[2],
                                      (const char (*)[80]) &dumpfiles[3]);
    if (info_fp == NULL)
    {
        fprintf(stderr, "cmio_in: failed to open info file %s\n", outfiles[2]);
        return -1;
    }

    // Read integration parameters
    simul_params_t params = read_params(infiles[2]);
    if (params == NULL)
    {
        fputs("cmio_in: failed to parse parameter file\n", stderr);
        fclose(info_fp);
        return -1;
    }
    free(params);

    fclose(info_fp);
    return 0;
}
