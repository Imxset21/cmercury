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

//! Loop limits in the Laguerre attempts
#ifndef NLAG1
#define NLAG1 50
#endif
#ifndef NLAG2
#define NLAG2 400
#endif

//! Struct for a Jacobi coordinate
struct jacobi_coord
{
    double x;
    double y;
    double z;
};

//! Struct for cvals from pp170
struct cvals
{
    double c1;
    double c2;
    double c3;
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
extern int drift_one(
    const int nbod,
    double* mass,
    struct jacobi_coord *initial_pos,
    struct jacobi_coord *vbles,
    const double mu,
    const double dt)
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
                               struct cvals *cvals)
    __attribute__((nonnull(2, 3)));

/**
 * @brief Returns the real root of cubic often found in solving kepler problem in universal variables.
 * @param dt time step (real scalar)
 * @param r0 Distance between `Sun' and paritcle
 * @param mu reduced mass of system (real scalar)
 * @param alpha Twice the binding energy (real scalar)
 * @param u  Vel. dot radial vector (real scalar)
 * @param s solution of cubi*eqn for the universal variable
 * @returns flag (= 0 if O.K.) (integer)
 * @author Martin Duncan  
 * @date March 12/93
 */
extern int drift_kepu_p3solve(
    const double dt,
    const double r0,
    const double mu,
    const double alpha,
    const double u,
    double *s)
    __attribute__((nonnull(6)));

/**
 * @brief Solves Kepler's eqn. in universal variables using NEWTON'S METHOD
 * @param[in,out] s universal variable
 * @param[in] dt time step (real scalor)
 * @param[in] r0 Distance between `Sun' and particle (real scalar)
 * @param[in] mu Reduced mass of system (real scalar)
 * @param[in] alpha energy (real scalar)
 * @param[in] u angular momentun (real scalar)
 * @param[in,out] fp f' from p170 (real scaloa)
 * @param[in,out] cvals c's from p171-172 (real scalors)
 * @returns 0 if converged; -1 otherwise
 * @author Hal Levison  
 * @date 31/3/98
 */
extern int drift_kepu_new(
    double *s,
    const double dt,
    const double r0,
    const double mu,
    const double alpha,
    const double u,
    double *fp,
    struct cvals *cvals)
    __attribute__((nonnull(1, 7, 8)));

/**
 * @brief Solves kepler's equation in universal variables using LAGUERRE'S METHOD
 * @param[in,out] s pointer to universal variable
 * @param[in] dt time step (real scalar)
 * @param[in] r0 Distance between `Sun' and paritcle (real scalar)
 * @param[in] mu Reduced mass of system (real scalar)
 * @param[in] alpha energy (real scalor)
 * @param[in] u angular momentun  (real scalar)
 * @param[in,out] fp f' from p170 (real scalars)
 * @param[in,out] cvals c's from p171-172 (real scalars)
 * @returns 0 if converged; -1 otherwise
 * @author Hal Levison  
 * @date 4/21/93
 */
extern int drift_kepu_lag(
    double *s,
    const double dt,
    const double r0,
    const double mu,
    const double alpha,
    const double u,
    double *fp,
    struct cvals *cvals)
    __attribute__((nonnull(1, 7, 8)));

/**
 * @brief Initial guess for solving kepler's equation using universal variables
 * @param[in] dt time step (real scalor)
 * @param[in] r0 Distance between `Sun' and particle (real scalor)
 * @param[in] mu Reduced mass of system (real scalor)
 * @param[in] alpha energy (real scalor)
 * @param[in] u angular momentun  (real scalor)
 * @returns initial guess for the value of universal variable
 * @author Hal Levison, Martin Duncan
 * @date 8/6/98
 */
extern double drift_kepu_guess(
    const double dt,
    const double r0,
    const double mu,
    const double alpha,
    const double u,
    double s);

/**
 * @brief Returns the value of the function f of which we are trying to find the root
 *        in universal variables.
 * @param[in] dt time step (real scalar)
 * @param[in] r0 Distance between `Sun' and particle (real scalar)
 * @param[in] mu Reduced mass of system (real scalar)
 * @param[in] alpha Twice the binding energy (real scalar)
 * @param[in] u Vel. dot radial vector (real scalar)
 * @param[in] s Approx. root of f 
 * @returns 0 if O.K., -1 otherwise
 * @author Martin Duncan  
 * @date 12/93
 */
extern double drift_kepu_fchk(
    const double dt,
    const double r0,
    const double mu,
    const double alpha,
    const double u,
    const double s);

/**
 * @brief Solves kepler's equation using universal variables.
 * @param[in] dt time step (real scalor)
 * @param[in] r0            ==>  Distance between `Sun' and paritcle (real scalor)
 * @param[in] mu            ==>  Reduced mass of system (real scalor)
 * @param[in] alpha         ==>  energy (real scalor)
 * @param[in] u             ==>  angular momentun  (real scalor)
 * @param[in,out] fp  f' from p170 (real scalors)
 * @param[in,out] cvals  c's from p171-172 (real scalors)
 * @returns 0 if converged; !=0 if not
 * @author Hal Levison
 * @date 2/3/93
 */
extern int drift_kepu(
    const double dt,
    const double r0,
    const double mu,
    const double alpha,
    const double u,
    double *fp,
    struct cvals *cvals)
    __attribute__((nonnull(6, 7)));

/**
 * @brief  Subroutine for solving kepler's equation in difference form for an
 * ellipse, given SMALL dm and SMALL eccentricity.
 * @see drift_dan
 * WARNING - BUILT FOR SPEED : DOES NOT CHECK HOW WELL THE ORIGINAL
 * EQUATION IS SOLVED! (CAN DO THAT IN THE CALLING ROUTINE BY
 * CHECKING HOW CLOSE (x - ec*s +es*(1.-c) - dm) IS TO ZERO.
 * @param[in] dm increment in mean anomaly M (real*8 scalar)
 * @param[in] es ecc. times sin of E_0 (real*8 scalar)
 * @param[in] ec ecc. times cos of E_0 (real*8 scalar)
 * @param[in,out] x solution to Kepler's difference eqn (real*8 scalar)
 * @param[in,out] sx sin of x (real*8 scalar)
 * @param[in,out] cx cos of x (real*8 scalar)
 */
extern void drift_kepmd(
    const double dm,
    const double es,
    const double ec,
    double *x,
    double *sx,
    double *cx);

/**
 * @brief Performs the Danby and decides which vbles to use
 * @param[in] nbod number of massive bodies (int scalar)
 * @param[in] mass mass of bodies (real array)
 * @param[in] dt0 time step
 * @param[in,out] position initial, final position in jacobi coord (real scalars)
 * @param[in,out] vbles initial, final position in jacobi coord (real scalars)
 * @returns integer flag - zero if satisfactory, non-zero if nonconvergence
 * @authors Hal Levison, Martin Duncan  
 * @date April 6/93 - MD adds dt and keeps dt0 unchanged
 */
extern int drift_dan(
    const int nbod,
    double* mass,
    struct jacobi_coord *position,
    struct jacobi_coord *vbles,
    const double mu,
    const double dt)
    __attribute__((nonnull(2, 3, 4)));

#endif /* DANBY_H */
