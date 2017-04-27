#include "lu_decomposition.h"
#include "src/constants.h"
#include <iostream>

LUDecomposition::LUDecomposition()
{
}


Vector LUDecomposition::solve(const Matrix &lhs, const Vector &rhs, Matrix &lu)
{
    int size = rhs.size();
    assert(lhs.rows() == size);
    assert(lhs.rows() == size);

    // decompose
    double maxA, absA;
    int    imax;

    Vector p(size); // permutation
    for (int i = 0; i < size; ++i) {
        p(i) = i;
    }

    lu = lhs;
    for (int i = 0; i < size; ++i) {
        maxA = 0.0;
        imax = i;
        for (int k = i; k < size; ++k) {
            if ((absA = fabs(lhs(k, i))) > maxA) {
                maxA = absA;
                imax = k;
            }
        }
        if (imax != i) {
            p.row(i).swap(p.row(imax));
            lu.row(i).swap(lu.row(imax));
        }
        for (int j = i + 1; j < size; ++j) {
            lu(j, i) = lu(j, i) / lu(i, i);
            for (int k = i + 1; k < size; ++k) {
                lu(j, k) = lu(j, k) - lu(j, i) * lu(i, k);
            }
        }
    }

    // solve
    Vector x = Vector::Zero(size);
    for (int i = 0; i < size; ++i) {
        x(i) = rhs(int(p(i)));
        for (int k = 0; k < i; ++k) {
            x(i) = x(i) - lu(i, k) * x(k);
        }
    }
    for (int i = size - 1; i >= 0; --i) {
        for (int k = i + 1; k < size; ++k) {
            x(i) = x(i) - lu(i, k) * x(k);
        }
        x(i) = x(i) / lu(i, i);
    }

    return x;
}
