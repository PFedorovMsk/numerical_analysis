#ifndef GAUSS_SEIDEL_SOLVER_H
#define GAUSS_SEIDEL_SOLVER_H

#include "src/constants.h"
#include "src/matrix.h"


class GaussSeidelSolver
{
public:
    GaussSeidelSolver();

    static Vector solve(const Matrix &lhs, const Vector &rhs, int &iterations,
                        double tol = Const::EPS);

    inline static bool removeZeroesFromDiagonal(Matrix &lhs, Vector &rhs);
};


#endif // GAUSS_SEIDEL_SOLVER_H
