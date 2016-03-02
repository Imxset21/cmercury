#ifndef CONFIG_H
#define CONFIG_H

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

//! Macro for DSIGN in FORTRAN77: If B >= 0 then the result is ABS(A), else it is -ABS(A). 
#define DSIGN(A, B) B >= 0 ? fabs(A) : -fabs(A)

#endif /* CONFIG_H */
