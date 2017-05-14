#include "parabolic_equation.h"

#include <cassert>


ParabolicEquation::ParabolicEquation(double _a, double _b, double _c, Function2 *const _f,
                                     double _alpha, double _beta, Function *const _phi1,
                                     double _gamma, double _delta, Function *const _phi2,
                                     Function *const _psi, double _xMax, double _tMax)
    : a(_a)
    , b(_b)
    , c(_c)
    , alpha(_alpha)
    , beta(_beta)
    , gamma(_gamma)
    , delta(_delta)
    , xMax(_xMax)
    , tMax(_tMax)
    , phi1(*_phi1)
    , phi2(*_phi2)
    , psi(*_psi)
    , f(*_f)
{
    assert(a > 0);
    assert(xMax > 0);
    assert(tMax > 0);
}

void ParabolicEquation::computeSteps(double &dx, double &dt, int &countbyX, int &countbyT) const
{
    assert(dx > 0);
    assert(dx <= 0.25 * xMax);
    assert(dt > 0);
    assert(dt <= 0.25 * tMax);

    if (dt > 0.5 * dx * dx / a) {
        dt = 0.5 * dx * dx / a;
    }

    countbyT = int(tMax / dt) + 1;
    countbyX = int(xMax / dx) + 1;

    dt = tMax / (countbyT - 1);
    dx = xMax / (countbyX - 1);
}

void ParabolicEquation::computeError(const Matrix &u, Function2 *const analitic, double dx,
                                     double dt, double &error)
{
    Matrix a = Matrix::Zero(u.rows(), u.cols());
    for (int k = 0; k < u.rows(); ++k) {
        for (int j = 0; j < u.cols(); ++j) {
            a(k, j) = (*analitic)(j * dx, k * dt);
        }
    }
    error = (u - a).norm();
}

void ParabolicEquation::solveExplicit(double &dx, double &dt, Matrix &u) const
{
    int J, K;
    computeSteps(dx, dt, J, K);

    u = Matrix::Zero(K, J);

    for (int j = 0; j < J; j++) {
        u(0, j) = psi(j * dx);
    }

    double s = a * dt / pow(dx, 2);
    double p = b * 0.5 * dt / dx;

    for (int k = 0; k < K - 1; k++) {
        //заполняем каждый временной слой, кроме граничных точек (0, k) и (L, k):
        for (int j = 1; j < J - 1; j++) {
            u(k + 1, j) = u(k, j - 1) * (s - p) + u(k, j) * (1.0 - 2.0 * s + c * dt) +
                          u(k, j + 1) * (s + p) + dt * f(j * dx, k * dt);
        }

        //граничные условия:
        double z  = 0.5 * pow(dx, 2) / a;
        double tk = (k + 1) * dt;
        u(k + 1, 0) =
            (alpha * (u(k + 1, 1) + z * u(k, 0) / dt + z * f(0.0, tk)) + phi1(tk) * (b * z - dx)) /
            (alpha * (1.0 + z / dt - c * z) + beta * (b * z - dx));
        u(k + 1, J - 1) = (gamma * (u(k + 1, J - 2) + z * u(k, J - 1) / dt + z * f(xMax, tk)) +
                           phi2(tk) * (b * z + dx)) /
                          (gamma * (1.0 + z / dt - c * z) + delta * (b * z + dx));
    }
}

void ParabolicEquation::solveImplicit(double &dx, double &dt, Matrix &u) const
{
    int J, K;
    computeSteps(dx, dt, J, K);

    u = Matrix::Zero(K, J);

    for (int j = 0; j < J; j++) {
        u(0, j) = psi(j * dx);
    }

    double s = a * dt / pow(dx, 2);
    double p = 0.5 * b * dt / dx;

    //векторы для прогонки :
    Matrix lhs = Matrix::Zero(J, J);
    Vector rhs = Vector::Zero(J);

    for (int k = 0; k < K - 1; k++) {
        for (int j = 1; j < J - 1; j++) {
            rhs(j) = -u(k, j) - dt * f(j * dx, (k + 1) * dt);
        }

        double z = 0.5 * pow(dx, 2) / a;

        for (int i = 1; i < J - 1; i++) {
            lhs(i, i - 1) = s - p;
            lhs(i, i)     = -1.0 - 2.0 * s + c * dt;
            lhs(i, i + 1) = s + p;
        }
        lhs(0, 0)         = alpha * (1.0 + z / dt - c * z) + beta * (b * z - dx);
        lhs(0, 1)         = -alpha;
        lhs(J - 1, J - 2) = -gamma;
        lhs(J - 1, J - 1) = gamma * (1.0 + z / dt - c * z) + delta * (b * z + dx);

        rhs(0) =
            alpha * z * (u(k, 0) / dt + f(0.0, (k + 1) * dt)) + phi1((k + 1) * dt) * (b * z - dx);

        rhs(J - 1) = gamma * z * (u(k, J - 1) / dt + f(0.0, (k + 1) * dt)) +
                     phi2((k + 1) * dt) * (b * z + dx);

        u.row(k + 1) = (lhs.colPivHouseholderQr().solve(rhs)).transpose();
    }
}
