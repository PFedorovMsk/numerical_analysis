#ifndef LU_DECOMPOSITION_H
#define LU_DECOMPOSITION_H

#include "src/matrix.h"


class LUDecomposition
{
public:
    LUDecomposition();

    static Vector solve(const Matrix &lhs, const Vector &rhs, Matrix &lu);
};


#endif // LU_DECOMPOSITION_H
