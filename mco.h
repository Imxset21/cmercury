/**
 * @file mco.h
 * @date Time-stamp: <2016-03-03 23:39:24 pedro>
 * @author Pedro Rittner, John E. Chambers, Gregory Tabak
 * @brief 
 */
#ifndef MCO_H
#define MCO_H

/**
 * @author John E. Chambers
 * @date 2 March 2001
 * @brief Converts coordinates with respect to the central body to barycentric
 * 
 */
extern void cmco_h2b(
    const int nbod,
    const double *restrict m,
    const double **restrict xh,
    const double **restrict vh,
    double **restrict x,
    double **restrict v)
    __attribute__((nonnull(2, 3, 4, 5, 6)));

#endif /* MCO_H */
