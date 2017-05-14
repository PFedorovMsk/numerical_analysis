#include "bisection_method.h"


BisectionMethod::BisectionMethod()
{
}

bool BisectionMethod::findRoot(const Function &f, double a, double b, double tol, int &iterations,
                               double &root)
{
    if (fabs(f(a)) < Const::EPS) {
        root = a;
        return true;
    }

    if (fabs(f(b)) < Const::EPS) {
        root = b;
        return true;
    }

    iterations = 0;
    while (b - a > tol) {
        ++iterations;

        double x = a + 0.5 * (b - a);
        if (f(a) * f(x) < 0) {
            b = x;
        } else {
            a = x;
        }
    }

    root = a + 0.5 * (b - a);
    if (fabs(f(root)) < tol) {
        return true;
    }
    return false;
}
