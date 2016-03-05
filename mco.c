#include "mco.h"

void cmco_h2b(
    const int nbod,
    const double *restrict m,
    const double **restrict xh,
    const double **restrict vh,
    double **restrict x,
    double **restrict v)
{
    double mtot = 0.0;
    x[0][0] = 0.0;
    x[0][1] = 0.0;
    x[0][2] = 0.0;
    v[0][0] = 0.0;
    v[0][1] = 0.0;
    v[0][2] = 0.0;

    // Calculate coordinates and velocities of the central body
    for (int j = 1; j < nbod; j++)
    {
        mtot += m[j];
        x[0][0] += m[j] * xh[j][0];
        x[0][1] += m[j] * xh[j][1];
        x[0][2] += m[j] * xh[j][2];
        v[0][0] += m[j] * vh[j][0];
        v[0][1] += m[j] * vh[j][1];
        v[0][2] += m[j] * vh[j][2];
    }

    const double temp = -1.0 / (mtot + m[0]);
    x[0][0] *= temp;
    x[0][1] *= temp;
    x[0][2] *= temp;
    v[0][0] *= temp;
    v[0][1] *= temp;
    v[0][2] *= temp;

    // Calculate the barycentric coordinates and velocities
    for (int j = 1; j < nbod; j++)
    {
        x[j][0] += x[0][0];
        x[j][1] += x[0][1];
        x[j][2] += x[0][2];
        v[j][0] += v[0][0];
        v[j][1] += v[0][1];
        v[j][2] += v[0][2];
    }
}
