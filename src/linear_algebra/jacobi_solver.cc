#include "jacobi_solver.h"
#include <iostream>


JacobiSolver::JacobiSolver()
{
}


Vector JacobiSolver::solve(const Matrix &lhs, const Vector &rhs, int &iterations, double tol)
{
    int size = rhs.size();
    assert(lhs.rows() == size);
    assert(lhs.rows() == size);

    Matrix a = lhs;
    Vector b = rhs;

    if (!removeZeroesFromDiagonal(a, b)) {
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
        x            = beta + alpha * x;
        double tk    = std::fabs((x - prevX).norm() * normAlpha / (1.0 - normAlpha));
        if (tk <= tol) {
            return x;
        }
        if (iterations >= 1000) {
            std::cout << "Warning: exit after 1000 iterations" << std::endl;
            return x;
        }
    }
}

bool JacobiSolver::removeZeroesFromDiagonal(Matrix &lhs, Vector &rhs)
{
    int size = rhs.size();
    for (int i = 0; i < size; ++i) {
        if (std::fabs(lhs(i, i)) <= Const::EPS) {
            bool swaped = false;
            for (int k = i + 1; k < size; ++k) {
                if (std::fabs(lhs(k, i)) > Const::EPS) {
                    lhs.row(k).swap(lhs.row(i));
                    rhs.row(k).swap(rhs.row(i));
                    swaped = true;
                    break;
                }
            }
            if (!swaped) {
                std::cout << "Error: matrix is singular" << std::endl;
                return false;
            }
        }
    }
    return true;
}
