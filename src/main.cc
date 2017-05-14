#include "src/mathematical_physics/parabolic_equation.h"

#include <iostream>

static constexpr double a     = 1;
static constexpr double b     = 3;
static constexpr double c     = -4;
static constexpr double alpha = 1;
static constexpr double beta  = 1;
static constexpr double gama  = 1;
static constexpr double delta = 1;


class F : public Function2
{
    double operator()(double /*x*/, double /*t*/) const override
    {
        return 0.0;
    }
};

class Phi1 : public Function
{
    double operator()(double t) const override
    {
        return std::exp((c - a) * t) * (std::cos(b * t) + std::sin(b * t));
    }
};

class Phi2 : public Function
{
    double operator()(double t) const override
    {
        return -1 * std::exp((c - a) * t) * (std::cos(b * t) + std::sin(b * t));
    }
};

class Psi : public Function
{
    double operator()(double x) const override
    {
        return std::sin(x);
    }
};

class Analitic : public Function2
{
    double operator()(double x, double t) const override
    {
        return std::exp((c - a) * t) * std::sin(x + b * t);
    }
};


int main()
{
    double tMax = 1;
    double xMax = Const::PI;

    /*(double _a, double _b, double _c, Function2 *const _f, double _alpha, double _beta,
                      Function *const _phi1, double _gamma, double _delta, Function *const _phi2,
       Function *const _psi,
                      double _xMax, double _tMax);*/
    F        f;
    Phi1     phi1;
    Phi2     phi2;
    Psi      psi;
    Analitic an;

    ParabolicEquation pe(a, b, c, &f, alpha, beta, &phi1, gama, delta, &phi2, &psi, xMax, tMax);

    double dt = 0.1;
    double dx = 0.2;
    double err;
    Matrix u;

    double step = 0.2;
    while (step >= 0.1) {
        dt = dx = step;
        pe.solveExplicit(dx, dt, u);
        pe.computeError(u, &an, dx, dt, err);

        std::cout << "dx    = " << dx << "\n";
        std::cout << "dt    = " << dt << "\n";
        std::cout << "size u: " << u.rows() << " x " << u.cols() << "\n";
        std::cout << "error = " << err << "\n\n";

        step *= 0.5;
    }
    return 0;
}
