/**
 * @file danby.h
 * @brief Solutions (partial or otherwise) for Danby-type drift equations
 * @date Time-stamp: <2016-03-01 16:22:24 pedro>
 * @author Pedro Rittner, John E. Chambers, Gregory Tabak
 * @copyright LGPL v3
 */
#ifndef DANBY_H
#define DANBY_H

//! Convergence criteria for danby
#ifndef DANBYAC
#define DANBYAC 1.0e-14
#endif
#ifndef DANBYB
#define DANBYB 1.0e-13
#endif

//! Struct for a Jacobi coordinate
struct jacobi_coord
{
    double x;
    double y;
    double z;
};

/**
  * @brief Calculates the danby-type drift for one particle
  * @param nbod number of massive bodies (int scalar)
  * @param mass mass of bodies (real array)
  * @param inital_pos initial position in jacobi coord 
  * @param vbles initial position in jacobi coord 
  * @param mu ???
  * @param dt time step
  * @returns integer (zero for successful step)
  * @authors Hal Levison, Martin Duncan 
  * @date 2/10/93
  *
  * This function does the danby-type drift for one particle, using 
  * appropriate vbles and redoing a drift if the accuracy is too poor 
  * (as indicated by a non-zero return value).
  */
extern int drift_one(int nbod,
                     double* mass,
                     struct jacobi_coord *initial_pos,
                     struct jacobi_coord *vbles,
                     double mu,
                     double dt)
    __attribute__((nonnull(2, 3, 4)));

/**
  * @brief Function for the calculation of stumpff functions
  * @see Danby p.172  equations 6.9.15
  * @param x argument
  * @params {c0,c1,c2,c3} c's from p171-172 (real scalors)
  * @author Hal Levison
  * @date 31/3/98
  */
extern void drift_kepu_stumpff(double x,
                               double *c0,
                               double *c1,
                               double *c2,
                               double *c3)
    __attribute__((nonnull(2, 3, 4, 5)));

#endif /* DANBY_H */
