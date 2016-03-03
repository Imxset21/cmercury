#include "danby.h"
#include "config.h"
#include <math.h>

int drift_one(
    struct jacobi_coord *initial_pos,
    struct jacobi_coord *vbles,
    const double mu,
    const double dt)
{
    int iflg = drift_dan(initial_pos, vbles, mu, dt);

    if (iflg)
    {
        for (int i = 1; i <= 10; i++)
        {
            double dttmp = dt / 10.0;
            iflg = drift_dan(initial_pos, vbles, mu, dttmp);

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

    while (fabs(x) > 0.1)
    {
        n++;
        x = x * .25;
    }

    const double x2 = x  * x;
    const double x3 = x  * x2;
    const double x4 = x2 * x2;
    const double x5 = x2 * x3;
    const double x6 = x3 * x3;

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
            *c0 = 2.0 * *c0 **c0 - 1.0;
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
    const double denom = (mu - alpha * r0) / 6.0;
    const double a2 = 0.5 * u / denom;
    const double a1 = r0 / denom;
    const double a0 = -dt / denom;

    const double q = (a1 - a2 * a2 / 3.0) / 3.0;
    const double r = (a1 * a2 - 3.0 * a0) / 6.0 - (pow(a2, 3)) / 27.0;
    const double sq2 = pow(q, 3) + pow(r, 2);

    if (sq2 >= 0.0)
    {
        const double sq = sqrt(sq2);
        double p2 = 0.0, p1 = 0.0;

        if ((r + sq) <= 0.0)
        {
            p1 = pow(-(-(r + sq)), (1.0 / 3.0));
        } else {
            p1 = pow((r + sq), (1.0 / 3.0));
        }

        if ((r - sq) <= 0.0)
        {
            p2 =  pow(-(-(r - sq)), (1.0 / 3.0));
        } else {
            p2 = pow((r - sq), (1.0 / 3.0));
        }

        *s = p1 + p2 - a2 / 3.0;
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
        cvals->c3 = cvals->c3 * _s * s2;
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
    double _s = *s, _fp = *fp;

    // const int NTMP = NLAG2 + 1;
    // To get close approch needed to take lots of iterations if alpha < 0
    const int ncmax = (alpha < 0.0) ? NLAG2 : NLAG2;

    const double ln = 5.0;

    // Start laguere's method
    for (int nc = 0; nc <= ncmax; nc++)
    {
        const double x = _s * _s * alpha;
        double c0 = 0.0;
        drift_kepu_stumpff(x, &c0, cvals);
        cvals->c1 = cvals->c1 * _s;
        cvals->c2 = cvals->c2 * _s * _s;
        cvals->c3 = cvals->c3 * _s * _s * _s;
        const double f = r0 * cvals->c1 + u * cvals->c2 + mu * cvals->c3 - dt;
        _fp = r0 * c0 + u * cvals->c1 + mu * cvals->c2;
        const double fpp = (-40.0 * alpha + mu) * cvals->c1 + u * c0;
        const double ds = -ln * f / (_fp + DSIGN(1.0, _fp) * sqrt(fabs((ln - 1.0) *
                                     (ln - 1.0) * _fp * _fp) - (ln - 1.0) * ln * f * fpp));
        _s = _s + ds;

        const double fdt = f / dt;

        // quartic convergence, Laguerre's method succeeded
        if (fdt * fdt < DANBYB * DANBYB)
        {
            *s = _s, *fp = _fp;
            return 0;
        }

    }

    *s = _s, *fp = _fp;
    return 2;
}

static inline void mco_sine(double *x, double *sx, double *cx)
{
    *x = (*x > 0) ? remainder(*x, TWOPI) : remainder(*x, TWOPI) + TWOPI;
    *cx = cos(*x);
    *sx = (*x > PI) ? -sqrt(1.0 - (*cx **cx)) : sqrt(1.0 - (*cx **cx));
}

double drift_kepu_guess(
    const double dt,
    const double r0,
    const double mu,
    const double alpha,
    const double u,
    double s)
{
    if (alpha > 0.0)
    {
        // Find initial guess for elliptic motion
        if ( dt / r0 <= 0.4)
        {
            s = dt / r0 - (dt * dt * u) / (2.0 * r0 * r0 * r0);
            return s;
        } else {
            const double a = mu / alpha;
            const double en = sqrt(mu / (a * a * a));
            const double ec = 1.0 - r0 / a;
            const double es = u / (en * a * a);
            const double e = sqrt(ec * ec + es * es);
            double y = en * dt - es;

            double sy = 0.0, cy = 0.0;
            mco_sine(&y, &sy, &cy);

            const double sigma = DSIGN(1.0, (es * cy + ec * sy));
            const double x = y + sigma * 0.85 * e;
            s = x / sqrt(alpha);
        }
    } else {
        // Find initial guess for hyperbolic motion.
        const int iflg = drift_kepu_p3solve(dt, r0, mu, alpha, u, &s);

        if (iflg)
        {
            s = dt / r0;
        }
    }

    return s;
}

double drift_kepu_fchk(
    const double dt,
    const double r0,
    const double mu,
    const double alpha,
    const double u,
    const double s)
{
    double c0 = 0.0;
    struct cvals my_cvals = {.c1 = 0.0, .c2 = 0.0, .c3 = 0.0};
    const double x = s * s * alpha;
    drift_kepu_stumpff(x, &c0, &my_cvals);
    my_cvals.c1 = my_cvals.c1 * s;
    my_cvals.c2 = my_cvals.c2 * s * s;
    my_cvals.c3 = my_cvals.c3 * s * s * s;
    const double f = r0 * my_cvals.c1 + u * my_cvals.c2 + mu * my_cvals.c3 - dt;
    return f;
}

int drift_kepu(
    const double dt,
    const double r0,
    const double mu,
    const double alpha,
    const double u,
    double *fp,
    struct cvals *cvals)
{
    int iflg = 0;
    double s = 0.0;

    // Store initial guess for possible use later in
    // laguerre's method, in case newton's method fails.
    s = drift_kepu_guess(dt, r0, mu, alpha, u, s);
    const double st = s;

    iflg = drift_kepu_new(&s, dt, r0, mu, alpha, u, fp, cvals);

    if (iflg)
    {
        const double fo = drift_kepu_fchk(dt, r0, mu, alpha, u, st);
        const double fn = drift_kepu_fchk(dt, r0, mu, alpha, u, s);

        if (fabs(fo) < fabs(fn))
        {
            s = st;
        }

        iflg = drift_kepu_lag(&s, dt, r0, mu, alpha, u, fp, cvals);
    }

    return iflg;
}

void drift_kepmd(
    const double dm,
    const double es,
    const double ec,
    double *x,
    double *sx,
    double *cx)
{
    const double A0 = 39916800.0, A1 = 6652800.0, A2 = 332640.0,
                 A3 = 7920.0, A4 = 110.0;
    double dx;
    double y;

    double _x = *x, _sx  = *sx, _cx = *cx;

    // Calc initial guess for root
    const double fac1 = 1.0 / (1.0 - ec);
    const double q = fac1 * dm;
    const double fac2 = es * es * fac1 - ec / 3.0;
    _x = q * (1.0 - 0.5 * fac1 * q * (es - q * fac2));

    // Excellent approx. to sin and cos of x for small x.
    y = _x * _x;
    _sx = _x * (A0 - y * (A1 - y * (A2 - y * (A3 - y * (A4 - y))))) / A0;
    _cx = sqrt(1.0 - _sx * _sx);

    // Compute better value for the root using quartic Newton method
    const double f = _x - ec * _sx + es * (1. - _cx) - dm;
    const double fp = 1. - ec * _cx + es * _sx;
    const double fpp = ec * _sx + es * _cx;
    const double fppp = ec * _cx - es * _sx;
    dx = -f / fp;
    dx = -f / (fp + 0.5 * dx * fpp);
    dx = -f / (fp + 0.5 * dx * fpp + 0.16666666666666666 * dx * dx * fppp);
    _x = _x + dx;;

    // Excellent approx. to sin and cos of x for small x.
    y = _x * _x;
    _sx = _x * (A0 - y * (A1 - y * (A2 - y * (A3 - y * (A4 - y))))) / A0;
    _cx = sqrt(1.0 - _sx * _sx);

    *x = _x, *sx = _sx, *cx = _cx;
    return;
}

int drift_dan(
    struct jacobi_coord *pos,
    struct jacobi_coord *vbles,
    const double mu,
    const double dt0)
{
    int iflg;

    struct cvals my_cvals = {0.0, 0.0, 0.0};
    double x,y,z,vx,vy,vz;
    double f,g,fdot;
    double gdot;
    double u,alpha,fp,r0,v0s;
    double a,asq,en;
    double dm,ec,es,esq,xkep;
    double fchk,s,c;

    // Set dt = dt0 to be sure timestep is not altered while solving for new coords.
	double dt = dt0;
	iflg = 0;
    r0 = sqrt(pos->x*pos->x + pos->y*pos->y + pos->z*pos->z);
    v0s = vbles->x*vbles->x + vbles->y*vbles->y + vbles->z*vbles->z;
    u = pos->x*vbles->x + pos->y*vbles->y + pos->z*vbles->z;
    alpha = 2.0*mu/r0 - v0s;
        
	if (alpha > 0.0)
    {
        a = mu/alpha;
        asq = a*a;
        en = sqrt(mu/(a*asq));
        ec = 1.0 - r0/a;
        es = u/(en*asq);
        esq = ec*ec + es*es;
        dm = dt*en - (int) (dt * en / TWOPI) * TWOPI;
        dt = dm/en;
        if ((dm*dm > 0.16) || (esq > 0.36))
        {
            goto calc_drift_kepyu;
        }

        if (esq*dm*dm < 0.0016)
        {
           drift_kepmd(dm, es, ec, &xkep, &s, &c);
	       fchk = (xkep - ec*s +es*(1.-c) - dm);

	       if(fchk*fchk > DANBYB)
           {
               iflg = 1;
               return iflg;
           }

           fp = 1. - ec*c + es*s;
           f = (a/r0) * (c-1.) + 1.;
           g = dt + (s-xkep)/en;
           fdot = - (a/(r0*fp))*en*s;
           gdot = (c-1.)/fp + 1.;

           x = pos->x*f + vbles->x*g;
           y = pos->y*f + vbles->y*g;
           z = pos->z*f + vbles->z*g;
           vx = pos->x*fdot + vbles->x*gdot;
           vy = pos->y*fdot + vbles->y*gdot;
           vz = pos->z*fdot + vbles->z*gdot;

           pos->x = x;
           pos->y = y;
           pos->z = z;
           vbles->x = vx;
           vbles->y = vy;
           vbles->z = vz;

	       iflg = 0;
	       return iflg;

        }
    }
             
calc_drift_kepyu:
    drift_kepu(dt,r0,mu,alpha,u, &fp, &my_cvals);

    if (iflg == 0)
    {
        f = 1.0 - (mu/r0)*my_cvals.c2;
        g = dt - mu*my_cvals.c3;
        fdot = -(mu/(fp*r0))*my_cvals.c1;
        gdot = 1. - (mu/fp)*my_cvals.c2;

        x = pos->x*f + vbles->x*g;
        y = pos->y*f + vbles->y*g;
        z = pos->z*f + vbles->z*g;
        vx = pos->x*fdot + vbles->x*gdot;
        vy = pos->y*fdot + vbles->y*gdot;
        vz = pos->z*fdot + vbles->z*gdot;

        pos->x = x;
        pos->y = y;
        pos->z = z;
        vbles->x = vx;
        vbles->y = vy;
        vbles->z = vz;
	}
    
    return iflg;
}
