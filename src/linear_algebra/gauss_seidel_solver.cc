#include "gauss_seidel_solver.h"
#include "jacobi_solver.h"
#include <iostream>


GaussSeidelSolver::GaussSeidelSolver()
{
}

Vector GaussSeidelSolver::solve(const Matrix &lhs, const Vector &rhs, int &iterations, double tol)
{
    int size = rhs.size();
    assert(lhs.rows() == size);
    assert(lhs.rows() == size);

    Matrix a = lhs;
    Vector b = rhs;

    if (!JacobiSolver::removeZeroesFromDiagonal(a, b)) {
        return Vector(0);
    }

    Matrix alpha = Matrix::Zero(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i != j) {
                alpha(i, j) = -a(i, j) / a(i, i);
            }
        }
    }
    Vector beta = Vector::Zero(size);
    for (int i = 0; i < size; ++i) {
        beta(i) = b(i) / a(i, i);
    }

    Vector x         = Vector::Zero(size);
    double normAlpha = std::fabs(alpha.norm());
    iterations       = 0;

    while (true) {
        ++iterations;
        Vector prevX = x;

        // Jacobi method:
        // x_new = beta + alpha * x_old;

        // Gauss-Seidel method:
        for (int i = 0; i < size; ++i) {
            x(i) = beta(i);
            for (int j = 0; j < size; ++j) {
                x(i) = x(i) + alpha(i, j) * x(j);
            }
        }

        double tk = std::fabs((x - prevX).norm() * normAlpha / (1.0 - normAlpha));
        if (tk <= tol) {
            return x;
        }
        if (iterations >= 1000) {
            std::cout << "Warning: exit after 1000 iterations" << std::endl;
            return x;
        }
    }
}
