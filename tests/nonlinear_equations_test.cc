#include "nonlinear_equations_test.h"

#include "src/nonlinear_equations/bisection_method.h"
#include "src/nonlinear_equations/newtons_method.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>


class Func1 : public Function
{
    double operator()(double x) const override
    {
        return pow(2.0, x) - pow(x, 2) - 0.5;
    }
};

class DF1 : public Function
{
    double operator()(double x) const override
    {
        return pow(2.0, x) * log(2.0) - 2.0 * x;
    }
};

class Func2 : public Function
{
    double operator()(double x) const override
    {
        return log(x + 2.0) - pow(x, 2.0);
    }
};

class DF2 : public Function
{
    double operator()(double x) const override
    {
        return 1.0 / (x + 2.0) + 2.0 * x;
    }
};

class Func3 : public Function
{
    double operator()(double x) const override
    {
        return sqrt(1 - pow(x, 2)) - pow(Const::E, x) + 0.1;
    }
};

class DF3 : public Function
{
    double operator()(double x) const override
    {
        return -pow(Const::E, x) - x / sqrt(1 - pow(x, 2.0));
    }
};


void testBisectionMethod()
{
    std::ofstream fout("test_bisection.txt");

    fout << "Test bisection method.\n\n";

    Func1 f1;
    Func2 f2;
    Func3 f3;

    std::vector<Function *> funcs(3);
    funcs[0] = &f1;
    funcs[1] = &f2;
    funcs[2] = &f3;

    std::vector<double> a(3), b(3), x(3);

    a[0] = 0.0;
    a[1] = 0.0;
    a[2] = -0.6;

    b[0] = 6.0;
    b[1] = 3.0;
    b[2] = 0.8;

    for (size_t j = 0; j < 3; ++j) {
        fout << "f" << j + 1 << "(x) = 0:\n\n";
        double tol = 0.001;
        int    it;
        for (int i = 0; i < 5; ++i) {
            fout << "tol = " << tol << ":\n";
            if (BisectionMethod::findRoot(*(funcs[j]), a[j], b[j], tol, it, x[j])) {
                fout << "x = " << x[j] << "\n";
                fout << "iterations: " << it << "\n";
            } else {
                fout << "no roots on (" << a[j] << ", " << b[j] << ")\n";
            }
            tol *= 0.1;
            fout << "\n";
        }
        fout << "\n************************************\n\n";
    }


    fout.close();
}

void testNewtonsMethod()
{
    std::ofstream fout("test_newtons.txt");

    fout << "Test Newton's method.\n\n";

    Func1 f1;
    Func2 f2;
    Func3 f3;

    DF1 df1;
    DF2 df2;
    DF3 df3;

    std::vector<Function *> funcs(3);
    funcs[0] = &f1;
    funcs[1] = &f2;
    funcs[2] = &f3;

    std::vector<Function *> df(3);
    df[0] = &df1;
    df[1] = &df2;
    df[2] = &df3;

    std::vector<double> x0(3), x(3);

    x0[0] = 1.0;
    x0[1] = 1.2;
    x0[2] = -0.2;

    for (size_t j = 0; j < 3; ++j) {
        fout << "f" << j + 1 << "(x) = 0:\n\n";
        double tol = 0.001;
        int    it;
        for (int i = 0; i < 5; ++i) {
            fout << "tol = " << tol << ":\n";
            if (NewtonsMethod::findRoot(*(funcs[j]), *(df[j]), x0[j], tol, it, x[j])) {
                fout << "x = " << x[j] << "\n";
                fout << "iterations: " << it << "\n";
            } else {
                fout << "no roots\n";
            }
            tol *= 0.1;
            fout << "\n";
        }
        fout << "\n************************************\n\n";
    }


    fout.close();
}

void testNonlinearEquations()
{
    testBisectionMethod();
    testNewtonsMethod();
}
