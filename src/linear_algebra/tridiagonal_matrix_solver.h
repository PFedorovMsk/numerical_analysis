#ifndef TRIDIAGONAL_MATRIX_SOLVER_H
#define TRIDIAGONAL_MATRIX_SOLVER_H

#include "src/matrix.h"


class TridiagonalMmatrixSolver
{
public:
    TridiagonalMmatrixSolver();

    static Vector solve(const Matrix &lhs, const Vector &rhs);
};


#endif // TRIDIAGONAL_MATRIX_SOLVER_H
