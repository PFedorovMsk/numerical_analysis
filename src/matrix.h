#ifndef MATRIX_H
#define MATRIX_H

#include "third_party/eigen/Dense"


using Matrix    = Eigen::MatrixXd;
using Vector    = Eigen::VectorXd;
using RowVector = Eigen::RowVectorXd;

using CMatrix    = Eigen::MatrixXcd;
using CVector    = Eigen::VectorXcd;
using CRowVector = Eigen::RowVectorXcd;


void MakeBlockVector(const Vector &vec1, const Vector &vec2, Vector &blockVec);
void MakeBlockMatrix(const Matrix &mat1, const Matrix &mat2, Matrix &blockMat,
                     bool vertical = false);


#endif // MATRIX_H
