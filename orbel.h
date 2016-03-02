/**
 * @file orbel.h
 * @brief Solutions (partial or otherwise) for Kepler's equations
 * @date Time-stamp: <2016-03-01 16:20:41 pedro>
 * @author Pedro Rittner, John E. Chambers, Gregory Tabak
 * @copyright LGPL v3
 */
#ifndef ORBEL_H
#define ORBEL_H

#include "config.h"

/**
 * @brief Solves the equivalent of Kepler's eqn. for a parabola given Q (Fitz. notation.)
 * @param q parabola mean anomaly as a real scalar
 * @returns anomaly as a real scalar
 * @date May 11, 1992
 * @author M. Duncan 
 * @see p. 70-72 of Fitzpatrick's book "Principles of Celestial Mechanics"
 *
 * For a parabola we can solve analytically.
 * Corrected for negative Q and use power series for small Q.
 */
extern double orbel_zget(double q);

/**
 * @brief Solves Kepler's eqn. for hyperbola using hybrid approach.  
 * @param e eccentricity anomaly. (real scalar)
 * @param capn hyperbola mean anomaly. (real scalar)
 * @returns eccentric anomaly. (real scalar)
 * @author M. Duncan 
 * @date May 26, 1992
 *
 * Uses power series for N in terms of F and Newton's method.
 * ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
 */
extern double orbel_flon(double e, double capn);


/**
 * @brief  Solves Kepler's eqn. for hyperbola using hybrid approach.  
 * @param e eccentricity anomaly. (real scalar)
 * @param n hyperbola mean anomaly. (real scalar)
 * @returns eccentric anomaly. (real scalar)
 * @author M. Duncan 
 * @date May 26,1992.
 * For abs(N) < 0.636*ecc -0.6 , use FLON.
 * For larger N, uses FGET.
 */
extern double orbel_fhybrid(double e, double n);

/**
 * @brief  Solves Kepler's eqn. for hyperbola using hybrid approach.  
 * @param e eccentricity anomaly. (real scalar)
 * @param capn hyperbola mean anomaly. (real scalar)
 * @returns eccentric anomaly. (real scalar)
 * @authors M. Duncan, JEC
 * @date May 11, 1992.
 * @see pp. 70-72 of Fitzpatrick's book "Principles of Cel. Mech. ".
 * @see Quartic convergence from Danby's book.
 */
extern double orbel_fget(double e, double capn);

#endif /* UTILS_H */
