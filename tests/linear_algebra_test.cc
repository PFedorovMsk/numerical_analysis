#include "linear_algebra_test.h"

#include "src/linear_algebra/gauss_seidel_solver.h"
#include "src/linear_algebra/gauss_solver.h"
#include "src/linear_algebra/jacobi_solver.h"
#include "src/linear_algebra/lu_decomposition.h"
#include "src/linear_algebra/tridiagonal_matrix_solver.h"

#include "src/linear_algebra/jacobi_eigenvalues.h"
#include "src/linear_algebra/qr_eigenvalues.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>


struct ItByTol {
    double tol;
    int    it;
};


void testLUDecomposition()
{
    std::ofstream fout("test_lu.txt");

    fout << "Test LU decomposition.\n\n";

    Matrix a(4, 4);
    a << 1, 2, -2, 6, -3, -5, 14, 13, 1, 2, -2, -2, -2, -4, 5, 10;

    Vector b(4);
    b << 24, 41, 0, 20;

    Matrix lu;

    Vector x = LUDecomposition::solve(a, b, lu);

    fout << "A =\n" << a << "\n\nb =\n" << b << "\n\nx =\n" << x << "\n\n* * * * * * * * * *\n\n";
    fout << "e = A * x - b =\n" << a * x - b << "\n\n* * * * * * * * * *\n\n";
    fout << "LU =\n" << lu << "\n";

    fout.close();
}

void testTridiagonalMatrix()
{
    std::ofstream fout("test_tridiagonal_matrix.txt");

    fout << "Test tridiagonal matrix algorithm.\n\n";

    Matrix a(5, 5);
    a << -11, -9, 0, 0, 0, 5, -15, -2, 0, 0, 0, -8, 11, -3, 0, 0, 0, 6, -15, 4, 0, 0, 0, 3, 6;

    Vector b(5);
    b << -122, -48, -14, -50, 42;

    Vector x;
    bool   ok = TridiagonalMmatrixSolver::solve(a, b, x);

    if (ok) {
        fout << "A =\n"
             << a << "\n\nb =\n"
             << b << "\n\nx =\n"
             << x << "\n\n* * * * * * * * * *\n\n";
        fout << "e = A * x - b =\n" << a * x - b << "\n";
    } else {
        fout << "The condition of a correctness or stability isn't satisfied\n";
    }

    fout.close();
}

void testJacobi()
{
    std::ofstream fout("test_jacobi.txt");

    fout << "Test Jacobi method.\n\n";

    Matrix a(4, 4);
    a << 19, -4, -9, -1, -2, 20, -2, -7, 6, -5, -25, 9, -3, 0, -9, 12;

    Vector b(4);
    b << 100, -5, 34, 69;

    fout << "A =\n" << a << "\n\nb =\n" << b << "\n\n";

    std::vector<ItByTol> itByTol;
    double               tol = 0.01;
    for (int i = 0; i < 7; ++i) {
        int    it = 0;
        Vector x  = JacobiSolver::solve(a, b, it, tol);
        itByTol.push_back({tol, it});

        fout << "\n\n* * * * * * * * * *\n\ntol = " << tol << ", iterations count: " << it;
        fout << "\n\nx =\n" << x << "\n\n";
        fout << "e = A * x - b =\n" << a * x - b;

        tol *= 0.01;
    }

    fout << "\n\n* * * * * * * * * *\n\nTolerance : Iterations\n";
    for (size_t i = 0; i < itByTol.size(); ++i) {
        fout << itByTol[i].tol << " : " << itByTol[i].it << "\n";
    }

    fout.close();
}

void testGaussSeidel()
{
    std::ofstream fout("test_gauss_seidel.txt");

    fout << "Test Gauss-Seidel method.\n\n";

    Matrix a(4, 4);
    a << 19, -4, -9, -1, -2, 20, -2, -7, 6, -5, -25, 9, -3, 0, -9, 12;

    Vector b(4);
    b << 100, -5, 34, 69;

    fout << "A =\n" << a << "\n\nb =\n" << b << "\n\n";

    std::vector<ItByTol> itByTol;
    double               tol = 0.01;
    for (int i = 0; i < 7; ++i) {
        int    it = 0;
        Vector x  = GaussSeidelSolver::solve(a, b, it, tol);
        itByTol.push_back({tol, it});

        fout << "\n\n* * * * * * * * * *\n\ntol = " << tol << ", iterations count: " << it;
        fout << "\n\nx =\n" << x << "\n\n";
        fout << "e = A * x - b =\n" << a * x - b;

        tol *= 0.01;
    }

    fout << "\n\n* * * * * * * * * *\n\nTolerance : Iterations\n";
    for (size_t i = 0; i < itByTol.size(); ++i) {
        fout << itByTol[i].tol << " : " << itByTol[i].it << "\n";
    }

    fout.close();
}

void testGauss()
{
    std::ofstream fout("test_gauss.txt");

    fout << "Test Gauss estimation.\n\n";

    Matrix a(4, 4);
    a << 19, -4, -9, -1, -2, 20, -2, -7, 6, -5, -25, 9, -3, 0, -9, 12;

    Vector b(4);
    b << 100, -5, 34, 69;

    Vector x = GaussSolver::solve(a, b);

    fout << "A =\n" << a << "\n\nb =\n" << b << "\n\nx =\n" << x << "\n\n* * * * * * * * * *\n\n";
    fout << "e = A * x - b =\n" << a * x - b << "\n";

    fout.close();
}


//--------------------------------------------------------------------------//


void testQR()
{
    std::ofstream fout("test_qr.txt");

    fout << "Test QR.\n\n";

    Matrix mat(3, 3);
    mat << 3, -7, -1, -9, -8, 7, 5, 2, 2;

    fout << "A =\n" << mat << "\n\n";

    std::vector<ItByTol> itByTol;
    std::vector<double>  tol = {0.1, 0.075, 0.05, 0.025, 0.01, 0.005};
    for (size_t i = 0; i < tol.size(); ++i) {
        CVector eigenValues;
        int     it;
        double  err;
        QREigenvalues::compute(mat, eigenValues, it, err, tol[i]);
        itByTol.push_back({tol[i], it});

        fout << "\n\n* * * * * * * * * *\n\ntol = " << tol[i] << ", err = " << err
             << ", iterations count: " << it;
        fout << "\n\neigen values:\n" << eigenValues << "\n";
    }

    fout << "\n\n* * * * * * * * * *\n\nTolerance : Iterations\n";
    for (size_t i = 0; i < itByTol.size(); ++i) {
        fout << itByTol[i].tol << " : " << itByTol[i].it << "\n";
    }

    fout.close();
}

void testJacobiEigenvalue()
{
    std::ofstream fout("test_jacobi_eigenvalue.txt");

    fout << "Test Jacobi eigenvalue.\n\n";

    Matrix mat(3, 3);
    mat << -7, 4, 5, 4, -6, -9, 5, -9, -8;

    fout << "A =\n" << mat << "\n\n";

    std::vector<ItByTol> itByTol;
    std::vector<double>  tol = {0.1, 0.075, 0.05, 0.025, 0.01, 0.005};
    for (size_t i = 0; i < tol.size(); ++i) {
        Vector eigenValues;
        Matrix eigenVectors;
        int    it;
        double err;
        JacobiEigenvalues::compute(mat, eigenValues, eigenVectors, it, err, tol[i]);
        itByTol.push_back({tol[i], it});

        fout << "\n\n* * * * * * * * * *\n\ntol = " << tol[i] << ", err = " << err
             << ", iterations count: " << it;
        fout << "\n\neigen values:\n"
             << eigenValues << "\n\neigen vectors:\n"
             << eigenVectors << "\n";
    }

    fout << "\n\n* * * * * * * * * *\n\nTolerance : Iterations\n";
    for (size_t i = 0; i < itByTol.size(); ++i) {
        fout << itByTol[i].tol << " : " << itByTol[i].it << "\n";
    }

    fout.close();
}

//--------------------------------------------------------------------------//

void testLinearAlgebra()
{
    std::cout << "Solve linear systems:\n";

    testLUDecomposition();
    std::cout << " -- LU-decompositions complete\n";

    testTridiagonalMatrix();
    std::cout << " -- Tridiagonal matrix algorithm complete\n";

    testJacobi();
    std::cout << " -- Jacobi method complete\n";

    testGaussSeidel();
    std::cout << " -- Gauss-Seidel method complete\n";

    testGauss();
    std::cout << " -- Gauss estimation complete\n\n";


    std::cout << "Eigen values problem:\n";

    testJacobiEigenvalue();
    std::cout << " -- Jacobi algorithm complete\n";

    testQR();
    std::cout << " -- QR-algorithm complete\n\n";

    std::cout << "End.\n";
}
