#include "gauss_solver.h"
#include "src/constants.h"
#include <iostream>


GaussEstimation::GaussEstimation()
{
}

void GaussEstimation::compute(Matrix &mat)
{
    int m = std::min(mat.rows(), mat.cols());

    for (int k = 0; k < m; ++k) {
        double maxValue = mat(k, k);
        int    maxIndex = k;
        for (int i = k + 1; i < m; ++i) {
            double tmp = std::fabs(mat(i, k));
            if (tmp > maxValue) {
                maxValue = tmp;
                maxIndex = i;
            }
        }
        if (std::fabs(mat(maxIndex, k)) < Const::EPS) {
            std::cout << "Error: matrix is singular" << std::endl;
            return;
        }
        mat.row(k).swap(mat.row(maxIndex));

        for (int i = k + 1; i < m; ++i) {
            double f = mat(i, k) / mat(k, k);
            for (int j = k + 1; j < mat.cols(); ++j) {
                mat(i, j) = mat(i, j) - f * mat(k, j);
            }
            mat(i, k) = 0.0;
        }
    }
}

void GaussEstimation::compute(Matrix &mat, Vector &vec)
{
    Matrix m;
    MakeBlockMatrix(mat, vec, m, false);
    compute(m);
    mat = Matrix(m.block(0, 0, mat.rows(), mat.cols()));
    vec = Vector(m.block(0, mat.cols(), vec.rows(), 1));
}


Vector GaussSolver::solve(const Matrix &lhs, const Vector &rhs)
{
    int size = rhs.size();
    assert(lhs.rows() == size);
    assert(lhs.rows() == size);

    Matrix mat;
    MakeBlockMatrix(lhs, rhs, mat, false);
    GaussEstimation::compute(mat);

    Vector x = Vector::Zero(size);
    for (int i = size - 1; i >= 0; --i) {
        x(i) = mat(i, size) / mat(i, i);
        for (int k = i - 1; k >= 0; --k) {
            mat(k, size) = mat(k, size) - mat(k, i) * x(i);
        }
    }

    return x;
}
