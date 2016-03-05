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

//! Maximum number of bodies
#ifndef NMAX
#define NMAX 2000
#endif

//! Maximum number of close-encounter minima monitored simultaneously
#ifndef CMAX
#define CMAX 50
#endif

//! Maximum number of messages in message.in
#ifndef NMESS
#define NMESS 200
#endif

//! An implausibly large number
#ifndef HUGE_VAL
#define HUGE_VAL 9.9e29
#endif

//! Maximum number of files that can be open at the same time
#ifndef NFILES
#define NFILES 50
#endif

//! Gaussian gravitational constant squared
#ifndef K2
#define K2 2.959122082855911e-4
#endif

//! Astronomical unit in cm
#ifndef AU
#define AU 1.4959787e13
#endif

//! Mass of the Sun in g
#ifndef MSUN
#define MSUN 1.9891e33
#endif

//! Macro for DSIGN in FORTRAN77: If B >= 0 then the result is ABS(A), else it is -ABS(A). 
#define SIGN(A, B) B >= 0 ? fabs(A) : -fabs(A)

#endif /* CONFIG_H */
