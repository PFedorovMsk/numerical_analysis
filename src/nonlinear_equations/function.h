#ifndef FUNCTION_H
#define FUNCTION_H


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


#endif // FUNCTION_H
