template <typename T>
Matrix2D<T> operator*(const double constant, Matrix2D<T> matrix)
{
    Matrix2D result = Matrix2D<T>(matrix.get_row(), matrix.get_column());

    for (int i = 0; i < result.get_row(); i++)
        for (int j = 0; j < result.get_column(); j++)
            result[i][j] = constant * matrix[i][j];

    return result;
}

template <typename T>
Matrix2D<T> operator*(Matrix2D<T> matrix, const double constant)
{
    return constant * matrix;
}

template <typename T>
Matrix2D<T> operator+(Matrix2D<T> L, Matrix2D<T> R)
{
    Matrix2D<T> result = Matrix2D<T>(L.get_row(), L.get_column());

    for (int i = 0; i < L.get_row(); i++)
        for (int j = 0; j < L.get_column(); j++)
            result[i][j] = L[i][j] + R[i][j];

    return result;
}

template <typename T>
Vector<T> operator+(Vector<T> L, Vector<T> R)
{
    Vector<T> result(L.get_row());

    for (int i = 0; i < result.get_row(); i++)
        result[i] = L[i] + R[i];

    return result;
}

template <typename T>
Matrix2D<T> operator*(Matrix2D<T> L, Matrix2D<T> R)
{
    Matrix2D result = Matrix2D<T>(L.get_row(), L.get_column());

    for (int i = 0; i < result.get_row(); i++)
    {
        for (int j = 0; j < result.get_column(); j++)
        {
            for (int k = 0; k < L.get_column(); k++)
                result[i][j] = result[i][j] + L[i][k] * R[k][j];
        }
    }

    return result;
}

template <typename T>
Vector<T> operator*(Vector<T> vect, Matrix2D<T> matrix)
{
    Vector result = Vector<T>(matrix.get_column());

    for(int i = 0; i < matrix.get_row(); i++)
        for(int j = 0; j < matrix.get_column(); j++)
            result[i] += matrix[i][j] * vect[j];

    return result;
}

template <typename T>
Matrix2D<T> operator*(Vector<T> L, Vector<T> R)
{

}
