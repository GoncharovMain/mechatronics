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
    if (L.get_row() != R.get_row() || L.get_column() != R.get_column())
        throw std::string { "Matrix size mismatch." };

    Matrix2D<T> result = Matrix2D<T>(L.get_row(), L.get_column());

    for (int i = 0; i < L.get_row(); i++)
        for (int j = 0; j < L.get_column(); j++)
            result[i][j] = L[i][j] + R[i][j];

    return result;
}

template <typename T>
Vector<T> operator+(Vector<T> L, Vector<T> R)
{
    if (L.get_row() != R.get_row())
        throw std::string{"Vector size mismatch."};

    Vector<T> result(L.get_row());

    for (int i = 0; i < result.get_row(); i++)
        result[i] = L[i] + R[i];

    return result;
}

template <typename T>
Matrix2D<T> operator*(Matrix2D<T> L, Matrix2D<T> R)
{
    if (L.get_column() != R.get_row())
        throw std::string{ "Matrix size mismatch. The condition left_column: " +
                        std::to_string(L.get_column()) + " == right_row: " +
                        std::to_string(R.get_row()) + " is not met." };

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
    if (vect.get_row() != matrix.get_column())
        throw std::string{ "Vect and Matrix size mismatch. The condition left_row: " +
                            std::to_string(vect.get_row()) + " == right_row: " +
                            std::to_string(matrix.get_row()) + " is not met." };
                            
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
