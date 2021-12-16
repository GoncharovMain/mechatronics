// MS_kinematis_models.h
//----------------------------------------------
// Модуль формирует кинематическую модель
// манипуляционной системы
//----------------------------------------------
#include <math.h>
#include <iostream>

using namespace std;

// DMatrixAry() - формирует массив матриц
// дифференцирования
// матриц преобразования однородных координат.
// Входные параметры:
// N - число звеньев;
// beta[Nmax] - вектор опред. вид и
// последовательность кинематических пар.
// Выходные параметры:
// D[][4][Nmax] - массив матриц (4х4)
// дифференцирования.
int DMatrixAry(int N,               //_________________in
               int beta[Nmax],      //________in
               double D[][4][Nmax]) //___ou
{
    int i, j, k;
    if (N >= Nmax)
        return 1; // error
    else
        for (k = 0; k < N; k++)
        {
            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                    D[i][j][k] = 0.0;
            if (beta[k] == 1)
            {
                D[0][1][k] = -1.0;
                D[1][0][k] = 1.0;
            }
            else if (beta[k] == 0)
                D[2][3][k] = 1.0;
            else
                return 2; // error
        }                 // end for(k)
    return 0;
} // end DMatrixAry()
//----------------------------------------------
// DAijMatrixAry() - формирует массив
// дифференцированных матриц Aij.
// Входные параметры:
// N - число звеньев;
// D[][4][Nmax] - массив матриц D(4х4)
// дифференцирования,
// Aij[][4][Nmax]- массив матриц Aij(4х4).
// Выходные параметры:
// DAij[][4][Nmax]- массив матриц
// dAij/dqj=AijD(4х4).
int DAijMatrixAry(int N,                  //___________________in
                  double D[][4][Nmax],    //__in
                  double Aij[][4][Nmax],  //__in
                  double DAij[][4][Nmax]) //__ou
{
    int k;
    for (k = 0; k < N; k++)
        MultiMatrix(k, Aij, k, D, k, DAij);
    return 0;
} // end DAijMatrix()
//----------------------------------------------
// MultiMatrix() - выполняет умножение
// матриц(4х4).
// Входные параметры:
// ML[][4] - левая матрица(4х4),
// MR[][4] - правая матрица(4х4).
// Выходные параметры:
// M[][4] - матрица(4х4).
int MultiMatrix(double ML[][4], //____in
    double MR[][4], //____in
    double M [][4]) //____ou
{
    int i, j, l;
    for (i=0; i<4; i++ )
        for (j=0; j<4; j++ )
        {
            M[i][j]= 0.;
            for (l=0; l<4; l++ )
                M[i][j]= M[i][j]+ML[i][l]*MR[l][j];
        }
    return 0;
} //end MultiMatrix()
//-------------------------------------------------------
// DiAokMatrix() - вычисляет матрицу
// dAok/dqi(4х4).
// Входные параметры:
// k - номер матрицы Aok,
// i - номер обобщённой координаты qi
// по которой выполняется дифференцирование
// матрицы Aok,
// Aok[][4][Nmax]- массив матриц (4х4)
// преобразования однородных координат,
// Aij[][4][Nmax]- массив матриц Aij(4х4),
// DAij[][4][Nmax]- массив матриц
// dAij/dqj(4х4).
// Выходные параметры:
// DiAok[][4] - матрица dAok/dqi(4х4).
int DiAokMatrix(int k,                  //_____________________in
                int i,                  //_____________________in
                double Aok[][4][Nmax],  //___in
                double Aij[][4][Nmax],  //____in
                double DAij[][4][Nmax], //_____in
                double DiAok[][4])      //______ou
{
    int j, l1, l2, l3;
    double M[4][4];
    if (i > k)
    {
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                DiAok[l1][l2] = 0.0;
        return 0;
    }
    if (k == 1 && i == k) // и
    {
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                DiAok[l1][l2] = DAij[l1][l2][i - 1];
        return 0;
    }
    if (i == 1)
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                DiAok[l1][l2] = DAij[l1][l2][i - 1];
    else
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
            {
                DiAok[l1][l2] = 0.;
                for (l3 = 0; l3 < 4; l3++)
                    DiAok[l1][l2] =
                        DiAok[l1][l2] + Aok[l1][l3][i - 2] * DAij[l3][l2][i - 1];
            }
    for (j = i; j < k; j++)
    {
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
            {
                M[l1][l2] = 0.;
                for (l3 = 0; l3 < 4; l3++)
                    M[l1][l2] =
                        M[l1][l2] + DiAok[l1][l3] * Aij[l3][l2][j];
            }
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                DiAok[l1][l2] = M[l1][l2];
    }
    return 0;
} // end DiAokMatrix()
//----------------------------------------------
// DiiAokMatrix() - вычисляет матрицу
// d2Aok/dqidqi(4х4).
// Входные параметры:
// k - номер матрицы Aok,
// i - номер обобщённой координаты qi
// по которой выполняется дифференцирование
// матрицы Aok,
// Aok[][4][Nmax]- массив матриц (4х4)
// преобразования однородных координат,
// Aij[][4][Nmax]- массив матриц Aij(4х4),
// DAij[][4][Nmax]- массив матриц
// dAij/dqj(4х4),
// D[][4][Nmax]- массив матриц D(4х4)
// дифференцирования.
// Выходные параметры:
// DiiAok[][4] - матрица d2Aok/dqidqi(4х4).
int DiiAokMatrix(int k,                  //___________________in
                 int i,                  //____________________in
                 double Aok[][4][Nmax],  //__in
                 double Aij[][4][Nmax],  //___in
                 double DAij[][4][Nmax], //____in
                 double D[][4][Nmax],    //_____in
                 double DiiAok[][4])     //______ou
{
    int j, l1, l2, l3;
    double M[4][4];
    if (i > k)
    {
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                DiiAok[l1][l2] = 0.0;
        return 0;
    }
    if (k == 1 && i == k)
    {
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
            {
                DiiAok[l1][l2] = 0.;
                for (l3 = 0; l3 < 4; l3++)
                    DiiAok[l1][l2] = DiiAok[l1][l2] + DAij[l1][l3][i - 1] * D[l3][l2][i - 1];
            }
        return 0;
    }
    for (l1 = 0; l1 < 4; l1++)
        for (l2 = 0; l2 < 4; l2++)
        {
            M[l1][l2] = 0.;
            for (l3 = 0; l3 < 4; l3++)
                M[l1][l2] = M[l1][l2] + Aok[l1][l3][i - 2] * DAij[l3][l2][i - 1];
        }
    for (l1 = 0; l1 < 4; l1++)
        for (l2 = 0; l2 < 4; l2++)
        {
            DiiAok[l1][l2] = 0.;
            for (l3 = 0; l3 < 4; l3++)
                DiiAok[l1][l2] = DiiAok[l1][l2] + M[l1][l3] * D[l3][l2][i - 1];
        }
    for (j = i + 1; j <= k; j++)
    {
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
            {
                M[l1][l2] = 0.;
                for (l3 = 0; l3 < 4; l3++)
                    M[l1][l2] =
                        M[l1][l2] + DiiAok[l1][l3] * Aij[l3][l2][j - 1];
            }
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                DiiAok[l1][l2] = M[l1][l2];
    }
    return 0;
} // end DiiAokMatrix()
//----------------------------------------------
// DijAokMatrix() - вычисляет матрицу
// d2Aok/dqidqj(4х4).
// Входные параметры:
// k - номер матрицы Aok,
// i - номер обобщённой координаты qi,
// j - номер обобщённой координаты qj,
// по которым выполняется дифференцирование
// матрицы Aok,
// Aok[][4][Nmax]- массив матриц (4х4)
// преобразования однородных координат,
// Aij[][4][Nmax]- массив матриц Aij(4х4),
// DAij[][4][Nmax]- массив матриц
// dAij/dqj(4х4),
// D[][4][Nmax]- массив матриц D(4х4)
// дифференцирования.
// Выходные параметры:
// DijAok[][4] - матрица d2Aok/dqidqj(4х4).
int DijAokMatrix(int k,                  //___________________in
                 int i,                  //____________________in
                 int j,                  //_____________________in
                 double Aok[][4][Nmax],  //__in
                 double Aij[][4][Nmax],  //___in
                 double DAij[][4][Nmax], //____in
                 double D[][4][Nmax],    //_____in
                 double DijAok[][4])     //______ou
{
    int l, l1, l2, l3;
    double M[4][4];
    if (i > k || j > k)
    {
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                DijAok[l1][l2] = 0.0;
        return 0;
    }
    if (i == j)
    {
        DiiAokMatrix(k, i, Aok, Aij, DAij, D, DijAok);
        return 0;
    }
    if (i > j)
    {
        l1 = i;
        i = j;
        j = l1;
    }
    for (l1 = 0; l1 < 4; l1++)
        for (l2 = 0; l2 < 4; l2++)
        {
            DijAok[l1][l2] = 0.;
            for (l3 = 0; l3 < 4; l3++)
                DijAok[l1][l2] = DijAok[l1][l2] + Aok[l1][l3][i - 2] * DAij[l3][l2][i - 1];
        }
    for (l = i + 1; l < j; l++)
    {
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
            {
                M[l1][l2] = 0.;
                for (l3 = 0; l3 < 4; l3++)
                    M[l1][l2] =
                        M[l1][l2] + DijAok[l1][l3] * Aij[l3][l2][l - 1];
            }
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                DijAok[l1][l2] = M[l1][l2];
    }
    for (l1 = 0; l1 < 4; l1++)
        for (l2 = 0; l2 < 4; l2++)
        {
            M[l1][l2] = 0.;
            for (l3 = 0; l3 < 4; l3++)
                M[l1][l2] =
                    M[l1][l2] + DijAok[l1][l3] * DAij[l3][l2][j - 1];
        }
    for (l = j + 1; l <= k; l++)
    {
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
            {
                DijAok[l1][l2] = 0.;
                for (l3 = 0; l3 < 4; l3++)
                    DijAok[l1][l2] = DijAok[l1][l2] + M[l1][l3] * Aij[l3][l2][l - 1];
            }
        if (l < k)
            for (l1 = 0; l1 < 4; l1++)
                for (l2 = 0; l2 < 4; l2++)
                    M[l1][l2] = DijAok[l1][l2];
    }
    return 0;
} // end DijAokMatrix()
//----------------------------------------------
// MultiConstMatrix() - выполняет умножение матриц(4х4)
// на константу.
// Входные параметры:
// ML[][4]- матрица(4х4),
// _const - константа,
// Выходные параметры:
// MR[][4]- матрица(4х4).
int MultiConstMatrix(double ML[][4], //______in
                     double _const,  //______in
                     double MR[][4]) //______ou
{
    int l1, l2;
    for (l1 = 0; l1 < 4; l1++)
        for (l2 = 0; l2 < 4; l2++)
            MR[l1][l2] = ML[l1][l2] * _const;
    return 0;
} // end MultiConstMatrix()
//----------------------------------------------
// MultiVectMatrix() - выполняет умножение
// матриц(4х4)на вектор.
// Входные параметры:
// M[][4]- матрица(4х4),
// VectL[4]- вектор(4х1).
// Выходные параметры:
// VectR[4]- вектор(4х1).
int MultiVectMatrix(double M[][4],   //________in
                    double VectL[4], //________in
                    double VectR[4]) //________ou
{
    int l1, l2;
    for (l1 = 0; l1 < 4; l1++)
    {
        VectR[l1] = 0.;
        for (l2 = 0; l2 < 4; l2++)
            VectR[l1] = VectR[l1] + M[l1][l2] * VectL[l2];

    }
    return 0;
} // end MultiVectMatrix()
//----------------------------------------------
// lin_speed()-вычисляет линейную скорость точки
// по её локальным координатам.
// Входные параметры:
// k - номер звена,
// l_coord[4] - локальные координаты
// выбранной точки,
// dq[Nmax] - вектор обобщённых скоростей,
// Aok[][4][Nmax]- массив матриц (4х4)
// преобразования однородных координат,
// Aij[][4][Nmax]- массив матриц Aij(4х4),
// DAij[][4][Nmax]- массив матриц
// dAij/dqj(4х4).
// Выходные параметры:
// l_speed[4] - линейная скорость выбранной
// точки k-го звена.
void print44(double M[][4])
{
    cout << endl;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
            cout << M[i][j] << "  ";
        cout << endl;
    }
    cout << endl;
}
int lin_speed(int k,                  //____________________in
              double l_coord[4],      //_________in
              double dq[Nmax],        //__________in
              double Aok[][4][Nmax],  //___in
              double Aij[][4][Nmax],  //____in
              double DAij[][4][Nmax], //_____in
              double l_speed[4])      //______ou
{
    int i, l1, l2;
    double M1[4][4], M2[4][4], DiAok[4][4];
    if (k > Nmax)
        return 1;
    for (l1 = 0; l1 < 4; l1++)
        for (l2 = 0; l2 < 4; l2++)
            M1[l1][l2] = 0.;
    for (i = 1; i <= k; i++)
    {
        DiAokMatrix(k, i, Aok, Aij, DAij, DiAok);

        MultiConstMatrix(DiAok, dq[i], M2);

        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                M1[l1][l2] = M1[l1][l2] + M2[l1][l2];
    }
    MultiVectMatrix(M1, l_coord, l_speed);
    return 0;
} // end lin_speed()
// lin_acceler() - вычисляет линейное ускорение
// точки
// по её локальным координатам.
// Входные параметры:
// k - номер звена,
// l_coord[4] - локальные координаты
// выбранной точки,
// dq[Nmax] - вектор обобщённых скоростей,
// ddq[Nmax] - вектор обобщённых ускорений,
// Aok[][4][Nmax]- массив матриц (4х4)
// преобразования однородных координат,
// Aij[][4][Nmax]- массив матриц Aij(4х4),
// D[][4][Nmax]- массив матриц D(4х4)
// дифференцирования.
// DAij[][4][Nmax]- массив матриц
// dAij/dqj(4х4).
// Выходные параметры:
// l_acceler[4]- линейная скорость выбранной
// точки k-го звена.
int lin_acceler(int k,                  //____________________in
                double l_coord[4],      //_________in
                double dq[Nmax],        //__________in
                double ddq[Nmax],       //___________in
                double Aok[][4][Nmax],  //__in
                double Aij[][4][Nmax],  //___in
                double D[][4][Nmax],    //____in
                double DAij[][4][Nmax], //_____in
                double l_acceler[4])    //______ou
{
    int i, j, l1, l2;
    double dqidqj;
    double M1[4][4], M2[4][4], M3[4][4];
    if (k > Nmax)
        return 1;
    for (l1 = 0; l1 < 4; l1++)
        for (l2 = 0; l2 < 4; l2++)
            M1[l1][l2] = 0.;

    for (i = 1; i <= k; i++)
    {
        DiAokMatrix(k, i, Aok, Aij, DAij, M2);
        MultiConstMatrix(M2, ddq[i], M3);
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                M1[l1][l2] = M1[l1][l2] + M3[l1][l2];
        for (j = 1; j <= k; j++)
        {
            DijAokMatrix(k, i, j, Aok, Aij, DAij, D, M2);
            dqidqj = dq[i] * dq[j];
            MultiConstMatrix(M2, dqidqj, M3);
            for (l1 = 0; l1 < 4; l1++)
                for (l2 = 0; l2 < 4; l2++)
                    M1[l1][l2] = M1[l1][l2] + M3[l1][l2];
        }
    }
    MultiVectMatrix(M1, l_coord, l_acceler);
    return 0;
} // end lin_acceler()
//----------------------------------------------
// angul_speed() -вычисляет угловую скорость
// звена.
// Входные параметры:
// k - номер звена,
// beta[Nmax] - вектор опред. вид и послед.
// кинем. пар,
// dq[Nmax] - вектор обобщённых скоростей,
// Aok[][4][Nmax]- массив матриц (4х4)
// преобразования однородных координат.
// Выходные параметры:
// a_speed[4] - угловая скорость k-го звена.
int angul_speed(int k,                 //____________________in
                int beta[Nmax],        //__________in
                double dq[Nmax],       //___________in
                double Aok[][4][Nmax], //_____in
                double a_speed[4])     //______ou
{
    int i, l1, l2;
    double w[4], wi[4];
    double M[4][4];
    if (k > Nmax)
        return 1;
    for (l1 = 0; l1 < 4; l1++)
    {
        wi[l1] = 0.;
        a_speed[l1] = 0.;
    }
    for (i = 1; i <= k; i++)
    {
        wi[2] = beta[i] * dq[i];
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                M[l1][l2] = Aok[l1][l2][i];
        MultiVectMatrix(M, wi, w);
        for (l1 = 0; l1 < 4; l1++)
            a_speed[l1] =
                a_speed[l1] + w[l1];
    }
    return 0;
} // end angul_speed()
//----------------------------------------------
// angul_acceler() - вычисляет угловое ускорение
// звена.
// Входные параметры:
// k - номер звена,
// beta[Nmax] - вектор опред. вид и послед.
// кинем. пар,
// dq[Nmax] - вектор обобщённых скоростей,
// ddq[Nmax] - вектор обобщённых ускорений,
// Aok[][4][Nmax]- массив матриц (4х4)
// преобразования однородных координат,
// Aij[][4][Nmax]- массив матриц Aij(4х4),
// DAij[][4][Nmax]- массив матриц
// dAij/dqj(4х4).
// Выходные параметры:
// a_acceler[4] - угловое ускорение k-го
// звена.

int angul_acceler(int k,                  //__________________in
                  int beta[Nmax],         //_______in
                  double dq[Nmax],        //________in
                  double ddq[Nmax],       //_________in
                  double Aok[][4][Nmax],  //_in
                  double Aij[][4][Nmax],  //__in
                  double DAij[][4][Nmax], //___in
                  double a_acceler[4])    //____ou
{
    int i, l, l1, l2;
    double w1[4], w2[4], w3[4], w4[4];
    double M1[4][4], M2[4][4], M3[4][4];
    if (k > Nmax)
        return 1;
    for (l1 = 0; l1 < 4; l1++)
    {
        w1[l1] = 0.;
        w2[l1] = 0.;
        a_acceler[l1] = 0.;
    }
    for (l = 1; l <= k; l++)
    {
        w1[2] = beta[l - 1] * dq[l - 1];
        w2[2] = beta[l - 1] * ddq[l - 1];
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                M1[l1][l2] = 0.;
        for (i = 1; i <= l; i++)
        {
            DiAokMatrix(l, i, Aok, Aij, DAij, M2);
            MultiConstMatrix(M2, dq[i], M3);
            for (l1 = 0; l1 < 4; l1++)
                for (l2 = 0; l2 < 4; l2++)
                    M1[l1][l2] = M1[l1][l2] + M3[l1][l2];
        }
        MultiVectMatrix(M1, w1, w3);
        for (l1 = 0; l1 < 4; l1++)
            for (l2 = 0; l2 < 4; l2++)
                M2[l1][l2] = Aok[l1][l2][l];
        MultiVectMatrix(M2, w2, w4);
        for (l1 = 0; l1 < 4; l1++)
            a_acceler[l1] = a_acceler[l1] + w3[l1] + w4[l1];
    }
    return 0;
} // end angul_acceler()