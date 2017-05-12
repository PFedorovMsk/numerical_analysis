#include "newtons_method.h"


NewtonsMethod::NewtonsMethod()
{
}

bool NewtonsMethod::findRoot(const Function &f, const Function &df, double x0, double tol, int &iterations,
                             double &root)
{
    iterations = 0;
    double x   = x0;
    do {
        x0 = x;
        x  = x0 - f(x0) / df(x0);
        ++iterations;
    } while (fabs(x - x0) > tol);

    root = x;
    if (fabs(f(root)) < tol) {
        return true;
    }
    return false;
}
