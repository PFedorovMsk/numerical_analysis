#ifndef JACOBI_SOLVER_H
#define JACOBI_SOLVER_H

#include "src/constants.h"
#include "src/matrix.h"


class JacobiSolver
{
public:
    JacobiSolver();

    static Vector solve(const Matrix &lhs, const Vector &rhs, int &iterations, double tol = Const::EPS);

    static bool removeZeroesFromDiagonal(Matrix &lhs, Vector &rhs);
};


#endif // JACOBI_SOLVER_H
