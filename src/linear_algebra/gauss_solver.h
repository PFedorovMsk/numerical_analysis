#ifndef GAUSS_SOLVER_H
#define GAUSS_SOLVER_H

#include "src/matrix.h"


class GaussEstimation
{
public:
    GaussEstimation();

    static void compute(Matrix &mat);
    static void compute(Matrix &mat, Vector &vec);
};


class GaussSolver
{
public:
    GaussSolver();

    static Vector solve(const Matrix &lhs, const Vector &rhs);
};


#endif // GAUSS_SOLVER_H
