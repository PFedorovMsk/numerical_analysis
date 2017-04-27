#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <limits>


class Const
{
public:
    static constexpr double PI    = 3.14159265358979323846264; // Число pi.
    static constexpr double E     = 2.71828182845904523536029; // Число e.
    static constexpr double GAMMA = 0.57721566490153286060651; // Постоянная Эйлера-Маскерони.
    static constexpr double EPS = std::numeric_limits<double>::epsilon(); // Машинное epsilon для типа double.
    static constexpr double MIN =
        std::numeric_limits<double>::min(); // Минимальное положительное значение для типа double.
    static constexpr double MAX =
        std::numeric_limits<double>::max(); // Максимальное возможное значение для типа double.
    static constexpr double INF =
        std::numeric_limits<double>::infinity(); // Любое значение типа double, превышающее MAX.
};


#endif // CONSTANTS_H
