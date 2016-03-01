#include "danby.h"

#include <math.h>

int drift_one(int nbod,
              double* mass,
              struct jacobi_coord *initial_pos,
              struct jacobi_coord *vbles,
              double mu,
              double dt)
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
                        double *c1,
                        double *c2,
                        double *c3)
{
    int n = 0;
    double xm = 0.1, x2, x3, x4, x5, x6;
    double _c0 = *c0, _c1 = *c1, _c2 = *c2, _c3 = *c3;

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

    _c2 = 1.147074559772972 - 11 * x6 - 2.087675698786810 - 9 * x5
        + 2.755731922398589 - 7 * x4  - 2.480158730158730 - 5 * x3
        + 1.388888888888889 - 3 * x2  - 4.166666666666667 - 2 * x + .50;

    _c3 = 7.647163731819816 - 13 * x6 - 1.605904383682161 - 10 * x5
        + 2.505210838544172 - 8 * x4  - 2.755731922398589 - 6 * x3
        + 1.984126984126984 - 4 * x2  - 8.333333333333333 - 3 * x
        + 1.666666666666667 - 1;

    _c1 = 1. - x * _c3;
    _c0 = 1. - x * _c2;

    if (n)
    {
        for (int i = n; i > 1; i--)
        {
            _c3 = (_c2 + _c0 * _c3) * .25;
            _c2 = _c1 * _c1 * 0.5;
            _c1 = _c0 * _c1;
            _c0 = 2.0 * _c0 * _c0 - 1.0;
            x = x * 4.0;
        }
    }

    // Update values passed by caller
    *c0 = _c0;
    *c1 = _c1;
    *c2 = _c2;
    *c3 = _c3;

    return;
}
