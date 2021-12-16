template<class C>
ArrayMatrix<C>::ArrayMatrix(int count)
{
    ArrayMatrix::count = count;
}

template<class C>
ArrayMatrix<C>::ArrayMatrix(int count, int row)
    : ArrayMatrix::ArrayMatrix(count)
{
    for (int i = 0; i < count; i++)
        ArrayMatrix::array.push_back(C(row));
}

template<class C>
ArrayMatrix<C>::ArrayMatrix(int count, int row, int column)
    : ArrayMatrix::ArrayMatrix(count)
{
    for (int i = 0; i < count; i++)
        ArrayMatrix::array.push_back(C(row, column));
}

template<class C>
C ArrayMatrix<C>::operator[](int n)
{
    return ArrayMatrix::array[n];
}

template<class C>
void ArrayMatrix<C>::print()
{
    for (C item : ArrayMatrix::array)
        item.print();
}
template<class C>
void ArrayMatrix<C>::print(std::string name)
{
    cout << "[" << name << "]:" << endl;
    ArrayMatrix::print();
}

