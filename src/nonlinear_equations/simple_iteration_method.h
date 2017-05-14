#ifndef SIMPLE_ITERATION_METHOD_H
#define SIMPLE_ITERATION_METHOD_H

#include "src/constants.h"
#include "src/math.h"
#include "src/nonlinear_equations/function.h"


class SimpleIterationMethod
{
public:
    SimpleIterationMethod();

    static bool findRoot(const Function &f, double x0, double tol, int &iterations, double &root);
    static bool findRoot(const MultivarFunction &f, const Vector &x0, double tol, int &iterations,
                         Vector &root);
};


#endif // SIMPLE_ITERATION_METHOD_H
