/**
 * @file mio.h
 * @date Time-stamp: <2016-03-03 11:56:03 pedro>
 * @author Pedro Rittner, John E. Chambers, Gregory Tabak
 * @brief String formatting functions for cmercury
 */

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
extern int** mio_spl(const size_t len, const char *str, int nsub)
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
extern char* mio_re2c(const double x, const double xmin, const double xmax);

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
extern int mio_out(
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
