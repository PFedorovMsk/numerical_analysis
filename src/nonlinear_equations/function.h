#ifndef FUNCTION_H
#define FUNCTION_H

#include "src/matrix.h"


class Function
{
public:
    Function()
    {
    }
    virtual ~Function()
    {
    }

    virtual double operator()(double) const = 0;
};


class MultivarFunction
{
public:
    MultivarFunction()
    {
    }
    virtual ~MultivarFunction()
    {
    }

    virtual Vector operator()(const Vector &) const = 0;
};


#endif // FUNCTION_H
