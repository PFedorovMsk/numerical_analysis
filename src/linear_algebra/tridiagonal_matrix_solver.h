#ifndef TRIDIAGONAL_MATRIX_SOLVER_H
#define TRIDIAGONAL_MATRIX_SOLVER_H

#include "src/matrix.h"


class TridiagonalMmatrixSolver
{
public:
    TridiagonalMmatrixSolver();

    static bool solve(const Matrix &lhs, const Vector &rhs, Vector &x);
    static bool solve(const Vector &a, const Vector &b, const Vector &c, const Vector &rhs,
                      Vector &x);
};


#endif // TRIDIAGONAL_MATRIX_SOLVER_H
