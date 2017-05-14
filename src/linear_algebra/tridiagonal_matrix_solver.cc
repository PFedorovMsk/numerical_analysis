#include "tridiagonal_matrix_solver.h"


TridiagonalMmatrixSolver::TridiagonalMmatrixSolver()
{
}

bool TridiagonalMmatrixSolver::solve(const Matrix &lhs, const Vector &rhs, Vector &x)
{
    int size = rhs.size();

    Vector a = Vector::Zero(size);
    Vector b = Vector::Zero(size);
    Vector c = Vector::Zero(size);

    a(0) = 0.0;
    b(0) = lhs(0, 0);
    c(0) = lhs(0, 1);

    for (int i = 1; i < size - 1; ++i) {
        a(i) = lhs(i, i - 1);
        b(i) = lhs(i, i);
        c(i) = lhs(i, i + 1);
    }

    a(size - 1) = lhs(size - 1, size - 2);
    b(size - 1) = lhs(size - 1, size - 1);
    c(size - 1) = 0.0;

    return solve(a, b, c, rhs, x);
}

bool TridiagonalMmatrixSolver::solve(const Vector &a, const Vector &b, const Vector &c,
                                     const Vector &rhs, Vector &x)
{
    int size = rhs.size();
    assert(a.size() == size);
    assert(b.size() == size);
    assert(c.size() == size);

    // is correct?
    for (int i = 0; i < size; ++i) {
        if (std::fabs(a(i)) < 0 || std::fabs(b(i)) <= 0 || std::fabs(c(i)) < 0) {
            return false;
        }
    }
    for (int i = 1; i < size - 1; ++i) {
        if (std::fabs(b(i)) < std::fabs(a(i)) + std::fabs(c(i))) {
            return false;
        }
    }
    if (std::fabs(b(0)) < std::fabs(c(0)) || std::fabs(b(size - 1)) < std::fabs(a(size - 1))) {
        return false;
    }
    // end

    Vector A = Vector::Zero(size);
    Vector B = Vector::Zero(size);

    A(0) = -c(0) / b(0);
    B(0) = rhs(0) / b(0);

    for (int i = 1; i < size - 1; ++i) {
        double denominator = b(i) + a(i) * A(i - 1);
        A(i)               = -c(i) / denominator;
        B(i)               = (rhs(i) - a(i) * B(i - 1)) / denominator;
    }
    A(size - 1) = 0.0;
    B(size - 1) =
        (rhs(size - 1) - a(size - 1) * B(size - 2)) / (b(size - 1) + a(size - 1) * A(size - 2));

    x           = Vector::Zero(size);
    x(size - 1) = B(size - 1);
    for (int i = size - 1; i > 0; --i) {
        x(i - 1) = A(i - 1) * x(i) + B(i - 1);
    }

    return true;
}
