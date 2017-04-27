#ifndef JACOBI_EIGENVALUES_H
#define JACOBI_EIGENVALUES_H

#include "src/constants.h"
#include "src/matrix.h"


class JacobiEigenvalues
{
public:
    JacobiEigenvalues();

    static void compute(const Matrix &mat, Vector &eigenValues, Matrix &eigenVectors, int &iterations, double &error,
                        double tol = Const::EPS);
};


#endif // JACOBI_EIGENVALUES_H
