template<typename T>
int Matrix2D<T>::get_row()
{
    return Matrix2D::row;
}

template<typename T>
int Matrix2D<T>::get_column()
{
    return Matrix2D::column;
}

template<typename T>
void Matrix2D<T>::print()
{
    cout << endl;
    for (int i = 0; i < Matrix2D::row; i++)
    {
        cout << "\t";
        for (int j = 0; j < Matrix2D::column; j++)
            cout << Matrix2D::matrix[i][j] << "\t";
        cout << endl;
    }
    cout << endl;
}
template<typename T>
void Matrix2D<T>::print(std::string name)
{
    cout << "[" << name << "]:" << endl;
    Matrix2D::print();
}

template<typename T>
Matrix2D<T> Matrix2D<T>::copy()
{
    std::vector<std::vector<T>> vectors;

    for (int i = 0; i < Matrix2D::row; i++)
    {
        std::vector<T> column;

        for (int j = 0; j < Matrix2D::column; j++)
            column.push_back((T)Matrix2D::matrix[i][j]);

        vectors.push_back(column);
    }

    return Matrix2D<T>(vectors);
}

template<typename T>
Matrix2D<T>::Matrix2D(int row, int column)
{
    Matrix2D::row = row;
    Matrix2D::column = column;

    Matrix2D::matrix = new T *[row];
    for (int i = 0; i < row; i++)
    {
        Matrix2D::matrix[i] = new T[column];
        for (int j = 0; j < column; j++)
            Matrix2D::matrix[i][j] = (T)0.0;
    }
}

template<typename T>
Matrix2D<T>::Matrix2D(initializer_list<initializer_list<T>> ilist)
    :Matrix2D::Matrix2D( row = ilist.size(), column = (*ilist.begin()).size() )
{
    matrix = new T *[this->row];

    for(int i = 0; i < this->row; i++)
    {
        matrix[i] = new T[column];
        auto row_list = (*(ilist.begin() + i));
        std::copy(row_list.begin(), row_list.end(), this->matrix[i]);
    }
}

template<typename T>
Matrix2D<T>::Matrix2D(std::vector<std::vector<T>> vectors)
    : Matrix2D::Matrix2D(row = vectors.size(), column = vectors[0].size())
{
    for (int i = 0; i < Matrix2D::row; i++)
        std::copy(vectors[i].begin(), vectors[i].end(), Matrix2D::matrix[i]);
}

template<typename T>
template <size_t N, size_t M>
Matrix2D<T>::Matrix2D(T (&array)[N][M])
{
    this->row = sizeof(array) / sizeof(array[0]);
    this->column = sizeof(array[0]) / sizeof(array[0][0]);

    matrix = new T *[this->row];

    for (int i = 0; i < this->row; i++)
        matrix[i] = array[i];
}