#ifndef NEWTONS_METHOD_H
#define NEWTONS_METHOD_H

#include "src/constants.h"
#include "src/math.h"
#include "src/nonlinear_equations/function.h"


class NewtonsMethod
{
public:
    NewtonsMethod();

    static bool findRoot(const Function &f, const Function &df, double x0, double tol, int &iterations, double &root);
};


#endif // NEWTONS_METHOD_H
