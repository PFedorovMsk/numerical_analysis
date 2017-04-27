#ifndef QR_EIGENVALUES_H
#define QR_EIGENVALUES_H

#include "src/constants.h"
#include "src/matrix.h"


class QREigenvalues
{
public:
    QREigenvalues();

    static void householder(Matrix &m, Matrix &Q, Matrix &R);
    static void compute(const Matrix &mat, CVector &eigenValues, int &iterations, double &error,
                        double tol = Const::EPS);
};


#endif // QR_EIGENVALUES_H
