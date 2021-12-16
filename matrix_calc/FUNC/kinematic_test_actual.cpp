//----------------------------------------------
// MS_kinematics_test.cpp
//----------------------------------------------
// Программа формирует геометрическую и
// кинематическую модели манипуляционных систем
// и вычисляет кинематические параметры
// (скорости и ускорения).
// Вывод результатов выполняется на "консоль".
//----------------------------------------------

#pragma hdrstop
#pragma argsused
#include <stdio.h>
#include "geometric_model.h"
#include "kinematic_model.h"
#include <iostream>
#include <iomanip>

using namespace std;

int PrintMatrix3N(int k, double M[][3][Nmax])
{
    int i, j;
    if (k <= 0 || k > Nmax)
    {
        printf("k=%d\n\n", k);
        return 1;
    }
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
            printf("%4.0f", M[i][j][k - 1]);
        printf("\n");
    }
    printf("\n");
    //_getch();
    return 0;
}
// Вывод к-ой матрицы(4х4)на "консоль"
int PrintMatrix4N(int k, double M[][4][Nmax])
{
    int i, j;
    if (k <= 0 || k > Nmax)
    {
        printf("k=%d\n\n", k);
        return 1;
    }
    printf("k=%d\n\n", k);
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
            printf("%4.0f", M[i][j][k - 1]);
        printf("\n");
    }
    printf("\n");
    //getch();
    return 0;
}
// Вывод матрицы(4х4)на "консоль"
int PrintMatrix4(double M[][4])
{
    int i, j;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
            printf("%4.0f", M[i][j]);
        printf("\n");
    }
    printf("\n");
    //getch();
    return 0;
}
// Вывод целочисленного вектора[Nmax] на "консоль"
int PrintVectInt(int N, int Vector[Nmax])
{
    int i;
    for (i = 0; i < N; i++)
        printf("%4d", Vector[i]);
    printf("\n");
    return 0;
}

// Вывод вещественного вектора[Nmax] на
// "консоль"
int PrintVectDou(int N, double Vector[Nmax])
{
    int i;
    for (i = 0; i < N; i++)
    {
        //printf("%4f", Vector[i]);
        cout << Vector[i] << "| ";
    }
    printf("\n");
    return 0;
}
//----------------------------------------------
int main(int argc, char *argv[])
{

    int i, j, k;
    int N = 3;

    int P[3] = {2, 1, 3};
    int beta[3];

    double q[3] = {0., M_PI / 2, 1.};

    double dq[3] = {1., 1., 1.};
    double ddq[3] = {1., 1., 1.};

    double Aii_     [4][4][Nmax];
    double Ai_j     [4][4][Nmax];
    double Aij      [4][4][Nmax];
    double Aok      [4][4][Nmax];
    double D        [4][4][Nmax];
    double DAij     [4][4][Nmax];

    double DiAok    [4][4];
    double DiiAok   [4][4];
    double DijAok   [4][4];

    double M        [3][3][Nmax];
    double L        [4][Nmax];

    double l_coord  [4];
    double l_speed  [4];
    double l_acceler[4];
    double a_speed  [4];
    double a_acceler[4];

#pragma region
    L[0][0] = 0.;
    L[1][0] = 0.;
    L[2][0] = 0.;

    L[0][1] = 0.;
    L[1][1] = 0.;
    L[2][1] = 1.;

    L[0][2] = 2.;
    L[1][2] = 0.;
    L[2][2] = 0.;

    L[0][3] = 0.;
    L[1][3] = 0.;
    L[2][3] = 3.;


    M[0][0][0] = 1.;
    M[0][1][0] = 0.;
    M[0][2][0] = 0.;
    M[1][0][0] = 0.;
    M[1][1][0] = 1.;
    M[1][2][0] = 0.;
    M[2][0][0] = 0.;
    M[2][1][0] = 0.;
    M[2][2][0] = 1.;

    M[0][0][1] = 0.;
    M[0][1][1] = 0.;
    M[0][2][1] = 1.;
    M[1][0][1] = 1.;
    M[1][1][1] = 0.;
    M[1][2][1] = 0.;
    M[2][0][1] = 0.;
    M[2][1][1] = 1.;
    M[2][2][1] = 0.;

    M[0][0][2] = 0.;
    M[0][1][2] = 0.;
    M[0][2][2] = 1.;
    M[1][0][2] = 1.;
    M[1][1][2] = 0.;
    M[1][2][2] = 0.;
    M[2][0][2] = 0.;
    M[2][1][2] = 1.;
    M[2][2][2] = 0.;

    M[0][0][3] = 0.;
    M[0][1][3] = 1.;
    M[0][2][3] = 0.;
    M[1][0][3] = 0.;
    M[1][1][3] = 0.;
    M[1][2][3] = 1.;
    M[2][0][3] = 1.;
    M[2][1][3] = 0.;
    M[2][2][3] = 0.;
#pragma endregion
    cout << setprecision(1) << fixed;

    // cout << "L[4][Nmax]" << endl;
    // for (i = 0; i < 4; i++)
    // {
    //     for (j = 0; j < Nmax; j++)
    //         cout << L[i][j] << " | ";
    //     cout << endl;
    // }

    // cout << "M[3][3][Nmax]" << endl;
    // for (i = 0; i < 3; i++)
    // {
    //     cout << i << endl;
    //     for (j = 0; j < 3; j++)
    //     {
    //         for (k = 0; k < Nmax; k++)
    //             cout << M[i][j][k] << " | ";
    //         cout << endl;
    //     }
    //     cout << endl;
    // }

    // cout << "before Aii_[4][4][Nmax]" << endl;
    // for (i = 0; i < 4; i++)
    // {
    //     cout << i << endl;
    //     for (j = 0; j < 4; j++)
    //     {
    //         for (k = 0; k < Nmax; k++)
    //         {
    //             Aii_[i][j][k] = 0.0;
    //             cout << Aii_[i][j][k]<< " | ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }

    // cout << "after Aii_[4][4][Nmax]" << endl;
    // for (i = 0; i < 4; i++)
    // {
    //     cout << i << endl;
    //     for (j = 0; j < 4; j++)
    //     {
    //         for (k = 0; k < Nmax; k++)
    //         {
    //             Aii_[i][j][k] = 0.0;
    //             cout << Aii_[i][j][k]<< " | ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }


    Aii_Matrix(N, L, M, Aii_);

    printf("\n[Aii*]:\n\n");
    for (k = 0; k <= N; k++)
        PrintMatrix4N(k, Aii_);

    BetaVect(N, P, beta);

    printf("\n{beta}:");
    PrintVectInt(N, beta);

    Ai_jMatrix(N, beta, q, Ai_j);

    printf("\n[Ai*j]:\n\n");
    for (k = 1; k <= N; k++)
        PrintMatrix4N(k, Ai_j);


    AijMatrix(N, Aii_, Ai_j, Aij);

    printf("\n[Aij]:\n\n");
    for (k = 1; k <= N; k++)
        PrintMatrix4N(k, Aij);

    AoiMatrix(N, Aij, Aok);
    //AokMatrix(N, Aij, Aok);

    printf("\n[Aok]:\n\n");
    for (k = 1; k <= N; k++)
        PrintMatrix4N(k, Aok);

    DMatrixAry(N, beta, D);

    printf("\n[D]:\n\n");
    for (k = 1; k <= N; k++)
        PrintMatrix4N(k, D);

    DAijMatrixAry(N, D, Aij, DAij);


    printf("\n[dAij/dqj]:\n\n");
    for (k = 1; k <= N; k++)
        PrintMatrix4N(k, DAij);


    for (k = 1; k <= N; k++)
        for (i = 1; i <= N; i++)
        {
            DiAokMatrix(k, i, Aok, Aij, DAij, DiAok);

            printf("\n[dAok/dqi]:\n\n");
            printf("k=%d", k);
            printf(", i=%d\n\n", i);
            PrintMatrix4(DiAok);
        }
    // //матрицы (4х4) d2Aok/dqidqi:
    for (k = 1; k <= N; k++)
        for (i = 1; i <= N; i++)
        {
            DiiAokMatrix(k, i, Aok, Aij, DAij, D, DiiAok);

            printf("\n[d2Aok/dqidqi]:\n\n");
            printf("k=%d", k);
            printf(", i=%d\n\n", i);
            PrintMatrix4(DiiAok);
        }

    for (k = 1; k <= N; k++)
        for (i = 1; i <= N; i++)
            for (j = 1; j <= N; j++)
            {
                DijAokMatrix(k, i, j, Aok, Aij, DAij, D, DijAok);

                printf("\n[d2Aok/dqidqj]:\n\n");
                printf("k=%d", k);
                printf(", i=%d", i);
                printf(", j=%d\n\n", j);
                PrintMatrix4(DijAok);
            }

    l_coord[0] = 1.;
    l_coord[1] = 1.;
    l_coord[2] = 1.;
    l_coord[3] = 1.;

    lin_speed(k = 3, l_coord, dq, Aok, Aij, DAij, l_speed);

    printf("\n{l_speed}:");
    PrintVectDou(4, l_speed);

    lin_acceler(k = 3, l_coord, dq, ddq, Aok, Aij, D, DAij, l_acceler);

    printf("\n{l_acceler}:");
    PrintVectDou(4, l_acceler);

    //printf("\n{a_speed}:");
    angul_acceler(k = 3, beta, dq, ddq, Aok, Aij, DAij, a_acceler);

    printf("\n{a_acceler}:");
    PrintVectDou(4, a_acceler);


    return 0;
}