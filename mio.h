/**
 * @file mio.h
 * @date Time-stamp: <2016-03-03 11:56:03 pedro>
 * @author Pedro Rittner, John E. Chambers, Gregory Tabak
 * @brief String formatting and I/O functions for cmercury
 */
#ifndef MIO_H
#define MIO_H

#include "config.h"

#include <stddef.h>
#include <stdio.h>

/**
 * @author John E. Chambers
 * Given a character string STRING, of length LEN bytes, the routine finds 
 * the beginnings and ends of NSUB substrings present in the original, and 
 * delimited by spaces. The positions of the extremes of each substring are 
 * returned in the array DELIMIT.
 * Substrings are those which are separated by spaces or the = symbol.
 */
extern int** cmio_spl(const size_t len, const char *str, int nsub)
    __attribute__((nonnull(2)));

/**
 * @author John E. Chambers
 * @brief Converts a double X, where XMIN <= X < XMAX, into an ASCII string
 * @returns ASCII string of 8 characters, using the new format compression
 *
 * X is first converted to base 224, and then each base 224 digit is converted 
 * to an ASCII character, such that 0 -> character 32, 1 -> character 33...
 * and 223 -> character 255.
 *
 * ASCII characters 0 - 31 (CTRL characters) are not used, because they
 * cause problems when using some operating systems.
 */
extern char* cmio_re2c(const double x, const double xmin, const double xmax);

/**
 * @author John E. Chambers
 * @brief Converts a double X into a 8-character ASCII string
 * @returns Converted string as a char pointer
 * 
 * N.B. X must lie in the range -1.e112 < X < 1.e112
 */
extern char* cmio_fl2c(const double x);

/**
 * @author John E. Chambers
 *
 * Writes output variables for each object to an output file. Each variable
 * is scaled between the minimum and maximum possible values and then
 * written in a compressed format using ASCII characters.
 * The output variables are:
 *  r = the radial distance
 *  theta = polar angle
 *  phi = azimuthal angle
 *  fv = 1 / [1 + 2(ke/be)^2], where be and ke are the object's binding and
 *                             kineti* energies. (Note that 0 < fv < 1).
 *  vtheta = polar angle of velocity vector
 *  vphi = azimuthal angle of the velocity vector
 *
 * If this is the first output (OPFLAG = -1), or the first output since the 
 * number of the objects or their masses have changed (OPFLAG = 1), then 
 * the names, masses and spin components of all the objects are also output.
 *
 * @returns 0 if successful, -1 otherwise
 * N.B. Each object's distance must lie between RCEN < R < RMAX
 */
extern int cmio_out(
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
    __attribute__((nonnull(7, 8, 9, 10, 11, 12, 13)));

/**
 * @author John E. Chambers
 * @brief Writes a progress report to the log file
 * Writes to stdout if you are running Mercury interactively
 */
extern int cmio_log(
    const double time,
    const double tstart,
    const double en[3],
    const double am[3],
    const int opt[8],
    char **mem,
    int lmem[NMESS])
    __attribute__(());

/**
 * @author John E. Chambers
 * @date 7 July 1999
 * @brief Converts from Julian day number to Julian/Gregorian Calendar dates
 *
 * Assumes the dates are those used by the English calendar.
 *
 * Algorithm taken from `Practical Astronomy with your calculator' (1988)
 * by Peter Duffett-Smith, 3rd edition, C.U.P.
 *
 * Algorithm for negative Julian day numbers (Julian calendar assumed) by
 * J. E. Chambers.
 *
 * N.B. The output date is with respect to the Julian Calendar on or before
 * ===  4th October 1582, and with respect to the Gregorian Calendar on or 
 *      after 15th October 1582.
 */
extern void cmio_jd2y(double jd0, int *year, int *month, double *day);


/**
 * @author John E. Chambers
 * @brief Writes out an error message and terminates Mercury.
 * @returns This function never returns
 */
extern void cmio_err(
    const int unit,
    const char *s1,
    const char *s2,
    const char *s3,
    const char *s4)
    __attribute__((noreturn));

/**
 * @author John E. Chambers
 * @brief Writes details of close encounter minima to an output file
 * 
 * It also decides how to continue the integration depending upon the
 * close-encounter option chosen by the user. Close encounter details are stored
 * until either 100 have been accumulated, or a data dump is done, at which
 * point the stored encounter details are also output.
 *
 * For each encounter, the routine outputs the time and distance of closest
 * approach, the identities of the objects involved, and the output
 * variables of the objects at this time. The output variables are:
 * expressed as
 *  r = the radial distance
 *  theta = polar angle
 *  phi = azimuthal angle
 *  fv = 1 / [1 + 2(ke/be)^2], where be and ke are the object's binding and
 *                             kineti* energies. (Note that 0 < fv < 1).
 *  vtheta = polar angle of velocity vector
 *  vphi = azimuthal angle of the velocity vector
 */
extern int cmio_ce(
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
    int ceflush);

/**
 * @author John E. Chambers
 * @brief Converts an ASCII string into a double X, where XMIN <= X < XMAX
 * @date 1 July 1999
 */
extern double cmio_c2re(
    char *restrict c,
    const double xmin,
    const double xmax,
    const size_t nchar)
    __attribute__((nonnull(1)));

#endif /* MIO_H */
