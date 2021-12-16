template <typename T>
Matrix2D<T> get_random(int row, int column)
{
    Matrix2D result = Matrix2D<T>(row, column);

    for (int i = 0; i < row; i++)
        for (int j = 0; j < column; j++)
            result[i][j] = rand() % 9 + 1;

    return result;
}

template <typename T>
Vector<T> get_random(int row)
{
    Vector result = Vector<T>(row);

    for (int i = 0; i < row; i++)
        result[i] = rand() % 9 + 1;

    return result;
}

int PrintMatrix4N(int k, ArrayMatrix<Matrix2D<double>> M)
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

int test_array_matrix()
{
    ArrayMatrix<Vector<double>> vectors(3, 3);
    ArrayMatrix<Matrix2D<double>> matrices(5, 3, 3);

    matrices[0][1][1] = 88;

    for (int i = 0; i < 5; ++i)
    {
        matrices[i].print();
        cout << i << endl;
    }


    return 0;
}

void test()
{

#pragma region_matrix

    Matrix2D left = get_random<int>(3, 3);
    Matrix2D right = get_random<int>(3, 3);

    cout << "left" << endl;
    left.print();

    cout << "right" << endl;
    right.print();

    cout << "result *:" << endl;
    Matrix2D result = left * right;
    result.print();

    cout << "result +:" << endl;
    result = left + right;
    result.print();

    cout << "result 5*M:" << endl;
    result = 5 * left;
    result.print();

    cout << result.get_column() << endl;


    Vector vect = get_random<int>(3);
    cout << "vect " << endl;
    vect.print();
#pragma endregion
    cout << endl << endl;
#pragma region
    int arr[3][3] = {
        {5, 5, 5},
        {6, 6, 6},
        {7, 7, 7}
    };

    Matrix2D matrix_init_list = Matrix2D<int>({
        {5, 5, 5},
        {6, 6, 6},
        {7, 7, 7}
    });

    Matrix2D matrix_init_arr = Matrix2D<int>(arr);

    cout << "matrix init list" << endl;
    matrix_init_list.print();

    cout << "matrix init arr" << endl;
    matrix_init_arr.print();


    Vector vector_init_list = Vector<int>({3, 8, 9});
    Vector vector_init_arr = Vector<int>(arr[0]);

    cout << "vector init list" << endl;
    vector_init_list.print();

    cout << "vector init arr" << endl;
    vector_init_arr.print();

    // cout << "error (0, 0):" << matrix_init_list(0, 0) << endl;
    cout << "[0, 0]:" << matrix_init_list[0][0] << endl;

    Matrix2D<double> arr_matrix[3]{
        get_random<double>(3, 3),
        get_random<double>(3, 3),
        get_random<double>(3, 3)
    };

    Matrix2D A = get_random<int>(3, 3);
    Matrix2D B = A.copy();

    cout << "part B" << endl;


    cout << "\t A" << endl;
    A.print();
    cout << "\t B" << endl;
    B.print();

    A[1][1] = 1000;
    cout << "\t A" << endl;
    A.print();
    cout << "\t B" << endl;
    B.print();
#pragma endregion

}