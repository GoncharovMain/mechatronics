template<typename T>
int Vector<T>::get_row()
{
    return Vector::row;
}

template<typename T>
void Vector<T>::print()
{
    cout << "\t";
    for (int i = 0; i < Vector::row; i++)
        cout << Vector::vector[i] << "\t";
    cout << endl;
}

template<typename T>
void Vector<T>::print(std::string name)
{
    cout << "[" << name << "]: ";
    Vector::print();
}

template<typename T>
Vector<T> Vector<T>::copy()
{
    Vector<T> copy = Vector(Vector::vector);
    return copy;
}

template<typename T>
Vector<T>::Vector(int row)
{
    Vector::row = row;
    Vector::vector = new T[row];
    for (int i = 0; i < row; i++)
        Vector::vector[i] = (T)0.0;
}

template<typename T>
template <size_t N>
Vector<T>::Vector(T (&array)[N])
    :Vector::Vector(sizeof(array) / sizeof(array[0]))
{
    this->vector = array;
}

template<typename T>
Vector<T>::Vector(initializer_list<T> ilist)
    :Vector::Vector(ilist.size())
{
    std::copy(ilist.begin(), ilist.end(), this->vector);
}