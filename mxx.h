/**
 * @file mxx.h
 * @brief 
 * @date Time-stamp: <2016-03-02 21:47:43 pedro>
 * @author Pedro Rittner, John E. Chambers, Gregory Tabak
 * @copyright LGPL v3
 */
#ifndef MXX_H
#define MXX_H

/**
 * Synchronizes the epochs of NBIG Big bodies (having a common epoch) and
 * NBOD-NBIG Small bodies (possibly having differing epochs), for an integration
 * using MERCURY.  The Small bodies are picked up in order starting with the one
 * with epoch furthest from the time, TSTART, at which the main integration will
 * begin producing output.  N.B. The synchronization integrations use Everhart's
 * RA15 routine.
 * @author John E. Chambers 
 */
extern void mxx_sync(
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
    int ngflag);

/**
 * @brief Sorts an array in-place, returning an array of new indices of elements
 */
extern int* mxx_sort(double *arr, int arr_size);

/**
 * @brief Calculates the Jacobi constant for massless particles
 * @author John E. Chambers, ErikSoft
 * @date 2 March 2001
 * @returns jac array
 * This assumes that *there are only 2 massive bodies (including the
 * central body) moving on circular orbits.
 * N.B. All coordinates and velocities must be heliocentric!!
 */
extern double* mxx_jac(
    double jcen[3],
    const int nbod,
    const int nbig,
    double *m,
    double **xh,
    double **vh)
    __attribute__((nonnull(4, 5, 6)));

/**
 * @author John E. Chambers
 * @brief Calculates the total energy and angular-momentum for a system
 * System has Masses M, coordinates X, velocities V and spin angular momenta S.
 * N.B. All coordinates and velocities must be with respect to the central body
 * @returns Total energy of the system
 */
extern double mxx_en(
    double jcen[3],
    const int nbod,
    const int nbig,
    double *m,
    double **xh,
    double **vh,
    double **s,
    double *e)
    __attribute__((nonnull(4, 5, 6, 7, 8)));

#endif /* MXX_H */
