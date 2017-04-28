#ifndef BISECTION_METHOD_H
#define BISECTION_METHOD_H

#include "src/constants.h"
#include "src/math.h"
#include "src/nonlinear_equations/function.h"


class BisectionMethod
{
public:
    BisectionMethod();

    static bool findRoot(const Function &f, double a, double b, double tol, int &iterations, double &root);
};


#endif // BISECTION_METHOD_H
