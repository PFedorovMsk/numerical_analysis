#ifndef PARABOLIC_EQUATION_H
#define PARABOLIC_EQUATION_H

#include "src/constants.h"
#include "src/matrix.h"
#include "src/nonlinear_equations/function.h"


class ParabolicEquation
{
public:
    ParabolicEquation(double _a, double _b, double _c, Function2 *const _f, double _alpha,
                      double _beta, Function *const _phi1, double _gamma, double _delta,
                      Function *const _phi2, Function *const _psi, double _xMax, double _tMax);

    void solveExplicit(double &h, double &dt, Matrix &u) const;

    void computeError(const Matrix &numeric, Function2 *const analitic, double dx, double dt,
                      double &error);

private:
    void computeSteps(double &dx, double &dt, int &countByX, int &countByT) const;

private:
    double a;
    double b;
    double c;

    double alpha;
    double beta;
    double gamma;
    double delta;

    double xMax;
    double tMax;

    Function &phi1;
    Function &phi2;
    Function &psi;

    Function2 &f;
};


#endif // PARABOLIC_EQUATION_H
