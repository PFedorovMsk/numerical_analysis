#include "simple_iteration_method.h"


SimpleIterationMethod::SimpleIterationMethod()
{
}

bool SimpleIterationMethod::findRoot(const Function &f, double x0, double tol, int &iterations, double &root)
{
    iterations = 0;
    double x   = x0;
    do {
        x0 = x;
        x  = x0 + f(x0);
        ++iterations;
    } while (fabs(x - x0) > tol);

    root = x;
    if (fabs(f(root)) < tol) {
        return true;
    }
    return false;
}

bool SimpleIterationMethod::findRoot(const MultivarFunction &f, const Vector &initX, double tol, int &iterations,
                                     Vector &root)
{
    iterations = 0;
    Vector x   = initX;
    Vector x0  = initX;
    do {
        x0 = x;
        x  = x0 + f(x0);
        ++iterations;
    } while (fabs((x - x0).norm()) > tol);

    root = x;

    return true;
}
