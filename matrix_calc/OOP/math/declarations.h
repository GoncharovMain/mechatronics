#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <math.h>

#define Nmax 10
#define M_PI 3.14159265358979323846

using namespace std;

template <typename T>
class Matrix2D
{
private:
    int row;
    int column;
    T **matrix;

public:
    Matrix2D(int row, int column);

    Matrix2D(initializer_list<initializer_list<T>> ilist);

    Matrix2D(std::vector<std::vector<T>> vectors);

    template <size_t N, size_t M>
    Matrix2D(T (&array)[N][M]);

    int get_row();
    int get_column();

    void print();
    void print(std::string name);
    Matrix2D copy();

    Matrix2D<T>& operator=(Matrix2D<T> equal_matrix)
    {
        this->row = equal_matrix.row;
        this->column = equal_matrix.column;

        for (int i = 0; i < this->row; i++)
            for (int j = 0; j < this->column; j++)
                this->matrix[i][j] = equal_matrix.matrix[i][j];

        return *this;
    }
    Matrix2D<T>& operator=(initializer_list<initializer_list<T>> ilist)
    {
        return Matrix2D<T>(ilist);
    }

    T *operator[](int x)
    {
        return matrix[x];
    }
};

template <typename T>
class Vector
{
private:
    int row;
    T *vector;
public:
    Vector();
    Vector(int row);

    template <size_t N>
    Vector(T (&array)[N]);

    Vector(initializer_list<T> ilist);

    int get_row();
    void print();
    void print(std::string name);
    Vector copy();

    Vector<T>& operator=(Vector<T> equal_vector)
    {
        this->row = equal_vector.row;

        for (int i = 0; i < this->row; i++)
            this->vector[i] = equal_vector.vector[i];

        return *this;
    }
    Vector<T>& operator=(std::initializer_list<T> ilist)
    {
        return Vector<T>(ilist);
    }

    T &operator[](int x)
    {
        return vector[x];
    }

};

template<class C>
class ArrayMatrix
{
private:
    int count;
    std::vector<C> array;
    ArrayMatrix(int count);
public:
    ArrayMatrix(int count, int row);
    ArrayMatrix(int count, int row, int column);

    C  operator[](int n);
    void print();
    void print(std::string name);
};

class Mechatronic
{
public:
    Mechatronic();

    int PrintMatrix4N(int k, ArrayMatrix<Matrix2D<double>> M);
    int PrintMatrix4N(ArrayMatrix<Matrix2D<double>> M);

    int Aii_Matrix(
            int N,
            ArrayMatrix<Vector<double>> Lvect,
            ArrayMatrix<Matrix2D<double>> Mii,
            ArrayMatrix<Matrix2D<double>> Aii_);

    int BetaVect(
            int N,
            Vector<int> Pvect,
            Vector<int> beta);

    int Ai_jMatrix(
            int N,
            Vector<int> beta,
            Vector<double> qj,
            ArrayMatrix<Matrix2D<double>> Ai_j);

    int AijMatrix(
            int N,
            ArrayMatrix<Matrix2D<double>> Aii_,
            ArrayMatrix<Matrix2D<double>> Ai_j,
            ArrayMatrix<Matrix2D<double>> Aij);

    int MultiMatrix(
            int kL,
            ArrayMatrix<Matrix2D<double>> ML,
            int kR,
            ArrayMatrix<Matrix2D<double>> MR,
            int k,
            ArrayMatrix<Matrix2D<double>> M);

    int AoiMatrix(
            int N,
            ArrayMatrix<Matrix2D<double>> Aij,
            ArrayMatrix<Matrix2D<double>> Aoi);

    int DMatrixAry(
            int N,
            Vector<int> beta,
            ArrayMatrix<Matrix2D<double>> D);

    int DAijMatrixAry(
                    int N,
                    ArrayMatrix<Matrix2D<double>> D,
                    ArrayMatrix<Matrix2D<double>> Aij,
                    ArrayMatrix<Matrix2D<double>> DAij);

    int DiAokMatrix(
            int k,
            int i,
            ArrayMatrix<Matrix2D<double>> Aok,
            ArrayMatrix<Matrix2D<double>> Aij,
            ArrayMatrix<Matrix2D<double>> DAij,
            Matrix2D<double> DiAok);

    int DiiAokMatrix(
            int k,
            int i,
            ArrayMatrix<Matrix2D<double>> Aok,
            ArrayMatrix<Matrix2D<double>> Aij,
            ArrayMatrix<Matrix2D<double>> DAij,
            ArrayMatrix<Matrix2D<double>> D,
            Matrix2D<double> DiiAok);

    int DijAokMatrix(
            int k,
            int i,
            int j,
            ArrayMatrix<Matrix2D<double>> Aok,
            ArrayMatrix<Matrix2D<double>> Aij,
            ArrayMatrix<Matrix2D<double>> DAij,
            ArrayMatrix<Matrix2D<double>> D,
            Matrix2D<double> DijAok);

    int lin_speed(int k,
            Vector<double> l_coord,
            Vector<double> dq,
            ArrayMatrix<Matrix2D<double>> Aok,
            ArrayMatrix<Matrix2D<double>> Aij,
            ArrayMatrix<Matrix2D<double>> DAij,
            Vector<double> l_speed);

    int lin_acceler(int k,
            Vector<double> l_coord,
            Vector<double> dq,
            Vector<double> ddq,
            ArrayMatrix<Matrix2D<double>> Aok,
            ArrayMatrix<Matrix2D<double>> Aij,
            ArrayMatrix<Matrix2D<double>> D,
            ArrayMatrix<Matrix2D<double>> DAij,
            Vector<double> l_acceler);

    int angul_acceler(int k,
                Vector<int> beta,
                Vector<double> dq,
                Vector<double> ddq,
                ArrayMatrix<Matrix2D<double>> Aok,
                ArrayMatrix<Matrix2D<double>> Aij,
                ArrayMatrix<Matrix2D<double>> DAij,
                Vector<double> a_acceler);
};

