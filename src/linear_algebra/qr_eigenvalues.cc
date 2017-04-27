#include "qr_eigenvalues.h"
#include <complex>

using Complex = std::complex<double>;
using CPair   = std::pair<Complex, Complex>;


QREigenvalues::QREigenvalues()
{
}

inline Matrix ComputeMinor(const Matrix &mat, int d)
{
    Matrix res = Matrix::Zero(mat.rows(), mat.cols());
    for (int i = 0; i < d; ++i) {
        res(i, i) = 1.0;
    }
    for (int i = d; i < mat.rows(); ++i) {
        for (int j = d; j < mat.cols(); j++) {
            res(i, j) = mat(i, j);
        }
    }
    return res;
}

void QREigenvalues::householder(Matrix &mat, Matrix &Q, Matrix &R)
{

    int m = mat.rows();
    int n = mat.cols();

    std::vector<Matrix> qv(m); // array of factor Q1, Q2, ... Qm
    Matrix              z(mat);
    Matrix              z1;

    for (int k = 0; k < n && k < m - 1; k++) {
        z1       = ComputeMinor(z, k);
        Vector x = z1.col(k);

        double a = std::fabs(x.norm());
        if (mat(k, k) > 0) {
            a = -a;
        }

        Vector e(x.size());
        for (int i = 0; i < e.size(); ++i) {
            e(i) = (i == k) ? 1 : 0;
        }
        e = x + a * e;
        e = e / std::fabs(e.norm());

        qv[size_t(k)] = Matrix::Identity(e.size(), e.size()) - 2 * e * e.transpose();
        z             = qv[size_t(k)] * z1;
    }

    Q = qv[0];
    for (int i = 1; i < n && i < m - 1; ++i) {
        z1 = qv[size_t(i)] * Q;
        Q  = z1;
    }

    R = Q * mat;
    Q = Matrix(Q.transpose());
}

inline CPair SolveSquareEquation(Complex a, Complex b, Complex c)
{
    Complex four(4, 0);
    Complex two(2, 0);
    Complex d  = std::sqrt(b * b - four * a * c);
    Complex x1 = (-b + d) / (two * a);
    Complex x2 = (-b - d) / (two * a);
    return CPair(x1, x2);
}

void QREigenvalues::compute(const Matrix &mat, CVector &eigenValues, int &iterations, double &error, double tol)
{
    eigenValues.resize(mat.rows());

    iterations = 0;
    Matrix a, q, r;
    Matrix a0 = mat;

    std::vector<CPair> roots0(a0.rows() - 1), roots1(a0.rows() - 1);

    while (true) {
        householder(a0, q, r);
        a  = q * r;
        a0 = r * q;
        ++iterations;

        roots0 = roots1;
        for (int j = 0; j < roots1.size(); ++j) {
            Complex A(1, 0);
            Complex B(a(j, j) + a(j + 1, j + 1), 0);
            Complex C(a(j, j) * a(j + 1, j + 1) - a(j, j + 1) * a(j + 1, j), 0);
            CPair   cp = SolveSquareEquation(A, B, C);
            roots1[j]  = cp;
        }

        double sum = 0;
        for (int m = 0; m < a.rows(); ++m) {
            for (int l = m + 1; l < a.rows(); ++l) {
                sum += pow(a(l, m), 2);
            }
        }
        error = sqrt(sum);

        for (int j = 0; j < a.rows(); ++j) {
            bool complexSwaped = false;
            if (j < roots1.size()) {
                Complex x01 = roots0[j].first;
                Complex x02 = roots0[j].second;
                Complex x11 = roots1[j].first;
                Complex x12 = roots1[j].second;
                double  e   = std::max(std::abs(x11 - x01), std::abs(x12 - x02));
                if (e <= tol) {
                    eigenValues(j)     = x11;
                    eigenValues(j + 1) = x12;
                    error              = std::min(e, error);
                    complexSwaped      = true;
                }
            }
            if (!complexSwaped) {
                Complex prev   = eigenValues[j];
                eigenValues[j] = Complex(a(j, j), 0);
                double e       = std::abs(eigenValues[j] - prev);
                if (e <= tol) {
                    error = std::min(e, error);
                }
            }
        }

        if (error <= tol) {
            return;
        }
    }
}
