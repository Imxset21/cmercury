/**
 * @file orbel.h
 * @brief Celestial orbit functions for cmercury
 * @date Time-stamp: <2016-03-01 14:43:08 pedro>
 * @author Pedro Rittner, John E. Chambers, Gregory Tabak
 * @copyright LGPL v3
 */
#ifndef ORBEL_H
#define ORBEL_H

//! A really small number
#ifndef TINY
#define TINY 4.0e-15
#endif

//! Max number of planets, including the Sun 
#ifndef NPLMAX
#define NPLMAX 202  
#endif

//! Max number of test particles
#ifndef NTPMAX
#define NTPMAX 2000
#endif 

//! Size of the test particle status flag
#ifndef NSTAT
#define NSTAT 3
#endif

//! Convergence criteria for danby
#ifndef DANBYAC
#define DANBYAC 1.0e-14
#endif
#ifndef DANBYB
#define DANBYB 1.0e-13
#endif

//! Loop limits in the Laguerre attempts
#ifndef NLAG1
#define NLAG1 50
#endif
#ifndef NLAG2
#define NLAG2 400
#endif

//! Numerical constant for Pi
#ifndef PI
#define PI 3.141592653589793
#endif
//! Numerical constant for 2*Pi
#ifndef TWOPI
#define TWOPI 2 * PI
#endif
//! Numerical constant for Pi/2
#ifndef PIBY2
#define PIBY2 PI / 2.0
#endif
//! Numerical constant for 180/Pi
#ifndef DEGRAD
#define DEGRAD 180.0 / PI
#endif

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

#endif /* UTILS_H */
