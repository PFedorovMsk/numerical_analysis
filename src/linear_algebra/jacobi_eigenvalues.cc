#include "jacobi_eigenvalues.h"
#include "src/math.h"


JacobiEigenvalues::JacobiEigenvalues()
{
}

inline double MaxNondiagonalElement(const Matrix &mat, int &row, int &col)
{
    row        = 1;
    col        = 0;
    double max = fabs(mat(row, col));

    for (int j = 1; j < mat.cols(); ++j) {
        for (int i = 0; i < j; ++i) {
            double value = fabs(mat(i, j));
            if (value > max) {
                value = max;
                row   = i;
                col   = j;
            }
        }
    }

    return max;
}

inline double SumSqr(const Matrix &mat)
{
    double sum = 0;
    for (int j = 1; j < mat.cols(); ++j) {
        for (int i = 0; i < j; ++i) {
            sum += pow(mat(i, j), 2);
        }
    }
    return sqrt(sum);
}


void JacobiEigenvalues::compute(const Matrix &mat, Vector &eigenValues, Matrix &eigenVectors, int &iterations,
                                double &error, double tol)
{
    Matrix a = mat;
    assert(a == a.transpose());
    int size = a.rows();

    iterations   = 0;
    error        = 2 * tol;
    eigenVectors = Matrix::Identity(size, size);

    while (true) {
        ++iterations;

        int i, j;
        MaxNondiagonalElement(a, i, j);
        double phi = 0.25 * Const::PI;
        if (fabs(a(i, i) - a(j, j)) > Const::EPS) {
            phi = 0.5 * atan(2 * a(i, j) / (a(i, i) - a(j, j)));
        }
        Matrix u = Matrix::Identity(size, size);
        u(i, i) = u(j, j) = cos(phi);
        u(i, j) = sin(phi);
        u(j, i) = -u(i, j);

        a     = u.transpose() * a * u;
        error = SumSqr(a);

        eigenVectors *= u;

        if (error < tol) {
            eigenValues = Vector(a.diagonal());
            return;
        }
    }
}
