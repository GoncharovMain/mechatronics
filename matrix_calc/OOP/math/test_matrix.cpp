#include <iostream>
#include <vector>
#include <iomanip>
#include <string>

#include "declarations.h"

#include "Matrix2D.cpp"
#include "Vector.cpp"
#include "ArrayMatrix.cpp"
#include "operators.cpp"
#include "test_functions.cpp"
#include "Mechatronic.cpp"

using namespace std;
#define Nmax 4

int test_array_matrix();
void test();
int test_kinematic();

int test_kinematic()
{
    int N = 3; // Nmax => 10

    Vector P = Vector({2, 1, 3});
    Vector beta = Vector<int>(3);

    Vector<double> q({0., M_PI / 2, 1.});

    Vector<double> dq({1., 1., 1.});
    Vector<double> ddq({1., 1., 1.});

    ArrayMatrix<Matrix2D<double>> Aii_(4, 4, Nmax);
    ArrayMatrix<Matrix2D<double>> Ai_j(4, 4, Nmax);
    ArrayMatrix<Matrix2D<double>> Aij(4, 4, Nmax);
    ArrayMatrix<Matrix2D<double>> Aok(4, 4, Nmax);
    ArrayMatrix<Matrix2D<double>> D(4, 4, Nmax);
    ArrayMatrix<Matrix2D<double>> DAij(4, 4, Nmax);

    Matrix2D<double> DiAok(4, 4);
    Matrix2D<double> DiiAok(4, 4);
    Matrix2D<double> DijAok(4, 4);

    ArrayMatrix<Matrix2D<double>> M(3, 3, Nmax);
    ArrayMatrix<Vector<double>> L(4, Nmax);

    Vector l_coord = Vector<double>({1.0, 1.0, 1.0, 1.0});

    Vector<double> l_speed(4);
    Vector<double> l_acceler(4);
    Vector<double> a_speed(4);
    Vector<double> a_acceler(4);

    L[2][1] = 1.;
    L[0][2] = 2.;
    L[2][3] = 3.;

    M[0][0][0] = 1.;
    M[1][1][0] = 1.;
    M[2][2][0] = 1.;
    M[0][2][1] = 1.;
    M[1][0][1] = 1.;
    M[2][1][1] = 1.;
    M[0][2][2] = 1.;
    M[1][0][2] = 1.;
    M[2][1][2] = 1.;
    M[0][1][3] = 1.;
    M[1][2][3] = 1.;
    M[2][0][3] = 1.;

    Mechatronic mechatronic;

    cout << setprecision(0) << fixed;

    // L.print("L");
    // M.print("M");
    // Aii_.print("Aii_");

    mechatronic.Aii_Matrix(N, L, M, Aii_);
    cout << "[Aii*]:" << endl;;
    mechatronic.PrintMatrix4N(Aii_);

    mechatronic.BetaVect(N, P, beta);
    beta.print("beta");

    mechatronic.Ai_jMatrix(N, beta, q, Ai_j);
    cout << "[Ai*j]:" << endl;
    mechatronic.PrintMatrix4N(Ai_j);

    mechatronic.AijMatrix(N, Aii_, Ai_j, Aij);
    cout << "[Aij]:" << endl;
    mechatronic.PrintMatrix4N(Aij);

    mechatronic.AoiMatrix(N, Aij, Aok);
    cout << "[Aok]:" << endl;
    mechatronic.PrintMatrix4N(Aok);

    mechatronic.DMatrixAry(N, beta, D);
    cout << "[D]:" << endl;
    mechatronic.PrintMatrix4N(D);

    mechatronic.AijMatrix(N, Aij, D, DAij); // or DAijMatrixAry(N, D, Aij, DAij)
    cout << "[DAij]:" << endl;
    mechatronic.PrintMatrix4N(DAij);


    for (int k = 1; k <= N; k++)
        for (int i = 1; i <= N; i++)
        {
            mechatronic.DiAokMatrix(k, i, Aok, Aij, DAij, DiAok);

            cout << endl << "[dAok/dqi]" << endl << endl;
            cout << "k = " << k << ", i = " << i << endl;

            DiAok.print();
        }

    for (int k = 1; k <= N; k++)
        for (int i = 1; i <= N; i++)
        {
            mechatronic.DiiAokMatrix(k, i, Aok, Aij, DAij, D, DiiAok);

            cout << endl << "[d2Aok/dqidqi]" << endl << endl;
            cout << "k = " << k << ", i = " << i << endl;

            DiiAok.print();
        }

    for (int k = 1; k <= N; k++)
        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
            {
                mechatronic.DijAokMatrix(k, i, j, Aok, Aij, DAij, D, DijAok);

                cout << endl << "[d2Aok/dqidqj]" << endl << endl;
                cout << "k = " << k << ", i = " << i << ", j = " << j << endl;

                DijAok.print();
            }


    mechatronic.lin_speed(3, l_coord, dq, Aok, Aij, DAij, l_speed);
    l_speed.print("l_speed");


    mechatronic.lin_acceler(3, l_coord, dq, ddq, Aok, Aij, D, DAij, l_acceler);
    l_acceler.print("l_acceler");


    mechatronic.angul_acceler(3, beta, dq, ddq, Aok, Aij, DAij, a_acceler);
    a_acceler.print("a_speed");


    return 0;
}


int main()
{

    //test_array_matrix();
    test();
    //test_kinematic();

    return 0;
}