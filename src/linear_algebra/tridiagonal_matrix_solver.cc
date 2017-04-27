#include "tridiagonal_matrix_solver.h"


TridiagonalMmatrixSolver::TridiagonalMmatrixSolver()
{
}

Vector TridiagonalMmatrixSolver::solve(const Matrix &lhs, const Vector &rhs)
{
    int size = rhs.size();
    assert(lhs.rows() == size);
    assert(lhs.rows() == size);

    Vector a = Vector::Zero(size);
    Vector b = Vector::Zero(size);

    a(0) = -lhs(0, 1) / lhs(0, 0);
    b(0) = rhs(0) / lhs(0, 0);

    for (int i = 1; i < size - 1; ++i) {
        double denominator = lhs(i, i) + lhs(i, i - 1) * a(i - 1);
        a(i)               = -lhs(i, i + 1) / denominator;
        b(i)               = (rhs(i) - lhs(i, i - 1) * b(i - 1)) / denominator;
    }
    a(size - 1) = 0.0;
    b(size - 1) = (rhs(size - 1) - lhs(size - 1, size - 2) * b(size - 2)) /
                  (lhs(size - 1, size - 1) + lhs(size - 1, size - 2) * a(size - 2));

    Vector x    = Vector::Zero(size);
    x(size - 1) = b(size - 1);
    for (int i = size - 1; i > 0; --i) {
        x(i - 1) = a(i - 1) * x(i) + b(i - 1);
    }

    return x;
}
