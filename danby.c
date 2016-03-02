#include "danby.h"
#include "config.h"
#include <math.h>

int drift_one(const int nbod,
              double* mass,
              struct jacobi_coord *initial_pos,
              struct jacobi_coord *vbles,
              const double mu,
              const double dt)
{
    int iflg = drift_dan(mu, initial_pos, vbles, dt);

    if (iflg)
    {
        for (int i = 1; i <= 10; i++)
        {
            double dttmp = dt / 10.0;
            iflg = drift_dan(mu, initial_pos, vbles, dttmp);

            if (iflg)
            {
                return iflg;
            }
        }
    }

    return 0;
}

void drift_kepu_stumpff(double x,
                        double *c0,
                        struct cvals *cvals)
{
    int n = 0;
    double xm = 0.1, x2, x3, x4, x5, x6;

    while (fabs(x) > xm)
    {
        n = n + 1;
        x = x * .25;
    }

    x2 = x  * x;
    x3 = x  * x2;
    x4 = x2 * x2;
    x5 = x2 * x3;
    x6 = x3 * x3;

    cvals->c2 = 1.147074559772972 - 11 * x6 - 2.087675698786810 - 9 * x5
        + 2.755731922398589 - 7 * x4  - 2.480158730158730 - 5 * x3
        + 1.388888888888889 - 3 * x2  - 4.166666666666667 - 2 * x + .50;

    cvals->c3 = 7.647163731819816 - 13 * x6 - 1.605904383682161 - 10 * x5
        + 2.505210838544172 - 8 * x4  - 2.755731922398589 - 6 * x3
        + 1.984126984126984 - 4 * x2  - 8.333333333333333 - 3 * x
        + 1.666666666666667 - 1;

    cvals->c1 = 1. - x * cvals->c3;
    *c0 = 1. - x * cvals->c2;

    if (n)
    {
        for (int i = n; i > 1; i--)
        {
            cvals->c3 = (cvals->c2 + *c0 * cvals->c3) * .25;
            cvals->c2 = cvals->c1 * cvals->c1 * 0.5;
            cvals->c1 = *c0 * cvals->c1;
            *c0 = 2.0 * *c0 * *c0 - 1.0;
            x = x * 4.0;
        }
    }

    return;
}

int drift_kepu_p3solve(const double dt,
                       const double r0,
                       const double mu,
                       const double alpha,
                       const double u,
                       double *s)
{
    double denom,a0,a1,a2,q,r,sq2;

	denom = (mu - alpha*r0)/6.0;
	a2 = 0.5*u/denom;
	a1 = r0/denom;
	a0 =-dt/denom;

	q = (a1 - a2*a2/3.0)/3.0;
	r = (a1*a2 -3.0*a0)/6.0 - (pow(a2,3))/27.0;
	sq2 = pow(q, 3) + pow(r, 2);

	if (sq2 >= 0.0)
    {
        double sq = sqrt(sq2), p2, p1;

        if ((r+sq) <= 0.0)
        {
            p1 = pow(-(-(r + sq)), (1.0/3.0));
        } else {
            p1 = pow((r + sq), (1.0/3.0));
        }
        if ((r-sq) <= 0.0)
        {
            p2 =  pow(-(-(r - sq)), (1.0/3.0));
        } else {
            p2 = pow((r - sq), (1.0/3.0));
        }

        *s = p1 + p2 - a2/3.0;
        return 0;
    }
    *s = 0;
    return -1;
}

int drift_kepu_new(double *s,
                   const double dt,
                   const double r0,
                   const double mu,
                   const double alpha,
                   const double u,
                   double *fp,
                   struct cvals *cvals)
{
    // Local variables for input/output
    double _s = *s, _fp = *fp;

    for (int nc = 0; nc <= 6; nc++)
    {
        const double s2 = _s * _s;
        const double x = s2 * alpha;
        double c0 = 0.0;
        drift_kepu_stumpff(x, &c0, cvals);
        cvals->c1 = cvals->c1 * _s;
        cvals->c2 = cvals->c2 * s2;
        cvals->c3 = cvals->c3 *_s * s2;
        const double f = r0 * cvals->c1 + u * cvals->c2 + mu * cvals->c3 - dt;
        _fp = r0 * c0 + u * cvals->c1 + mu * cvals->c2;
        const double fpp = (mu - r0 * alpha) * cvals->c1 + u * c0;
        const double fppp = (mu - r0 * alpha) * c0 - u * alpha * cvals->c1;
        double ds = - f / _fp;
        ds = - f / (_fp + .5 * ds * fpp);
        ds = -f / (_fp + .5 * ds * fpp + ds * ds * fppp * .1666666666666667);
        _s = _s + ds;
        const double fdt = f / dt;
        // quartic convergence, newton's method succeeded
        if (fdt * fdt < DANBYB * DANBYB)
        {
            *s = _s, *fp = _fp; 
            return 0;
        }
    }

    // newton's method failed
    *s = _s, *fp = _fp;
    return -1;
}

int drift_kepu_lag(double *s,
                   const double dt,
                   const double r0,
                   const double mu,
                   const double alpha,
                   const double u,
                   double *fp,
                   struct cvals *cvals)
{
    int ncmax;
    double _s = *s, _fp = *fp;

    // const int NTMP = NLAG2 + 1;
    // To get close approch needed to take lots of iterations if alpha < 0
    if(alpha < 0.0)
    {
        ncmax = NLAG2;
    } else {
        ncmax = NLAG2;
    }

    const double ln = 5.0;
    // Start laguere's method
    for (int nc = 0; nc <= ncmax; nc++)
    {
        const double x = _s*_s*alpha;
        double c0 = 0.0;
        drift_kepu_stumpff(x, &c0, cvals);
        cvals->c1 = cvals->c1*_s;
        cvals->c2 = cvals->c2*_s*_s;
        cvals->c3 = cvals->c3*_s*_s*_s;
        const double f = r0*cvals->c1 + u*cvals->c2 + mu*cvals->c3 - dt;
        _fp = r0*c0 + u*cvals->c1 + mu*cvals->c2;
        const double fpp = (-40.0*alpha + mu)*cvals->c1 + u*c0;
        const double ds = -ln*f/(_fp + DSIGN(1.0, _fp) * sqrt(fabs((ln - 1.0) *
                                                     (ln - 1.0) * _fp * _fp) - (ln - 1.0) * ln * f * fpp));
        _s = _s + ds;

        const double fdt = f / dt;

        // quartic convergence, Laguerre's method succeeded
        if (fdt * fdt < DANBYB*DANBYB)
        {
            *s = _s, *fp = _fp;
            return 0;
        }

    }

    *s = _s, *fp = _fp;
    return 2;
}
