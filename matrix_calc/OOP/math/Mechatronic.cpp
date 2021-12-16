Mechatronic::Mechatronic(){ }

int Mechatronic::PrintMatrix4N(int k, ArrayMatrix<Matrix2D<double>> M)
{
    cout << " k = " << k << endl << endl;
    if (k <= 0 || k > Nmax)
    {
        return 1;
    }
    for (int i = 0; i < 4; i++)
    {
        cout << "\t";
        for (int j = 0; j < 4; j++)
            cout << M[i][j][k - 1] << "\t";
        cout << endl;
    }
    cout << endl;
    return 0;
}

int Mechatronic::PrintMatrix4N(ArrayMatrix<Matrix2D<double>> M)
{
    for (int k = 1; k <= 3; k++) // k => N = 3 (have in test_matrix.cpp::main())
        Mechatronic::PrintMatrix4N(k, M);

    return 0;
}

int Mechatronic::Aii_Matrix(
                int N,
                ArrayMatrix<Vector<double>> Lvect,
                ArrayMatrix<Matrix2D<double>> Mii,
                ArrayMatrix<Matrix2D<double>> Aii_)
{
    if (N>=Nmax)
        return 1;

    for (int k = 0; k <= N ;k++)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Aii_[i][j][k] = Mii[i][j][k];

        Aii_[0][3][k] = Lvect[0][k];
        Aii_[1][3][k] = Lvect[1][k];
        Aii_[2][3][k] = Lvect[2][k];
        Aii_[3][0][k] = 0.0;
        Aii_[3][1][k] = 0.0;
        Aii_[3][2][k] = 0.0;
        Aii_[3][3][k] = 1.0;
    }
    return 0;
}

int Mechatronic::BetaVect(
            int N,
            Vector<int> Pvect,
            Vector<int> beta)
{
    if (N >= Nmax)
        return 1;

    for (int k = 0; k < N; k++)
        if (Pvect[k] == 1 || Pvect[k] == 2)
            beta[k] = 1;
        else if (Pvect[k] == 3)
            beta[k] = 0;
        else
            return 2;

    return 0;
}

int Mechatronic::Ai_jMatrix(
                int N,
                Vector<int> beta,
                Vector<double> qj,
                ArrayMatrix<Matrix2D<double>> Ai_j)
{
    for (int k = 0; k < N; k++)
    {
        Ai_j[0][0][k] = cos(beta[k] * qj[k]);
        Ai_j[1][0][k] = sin(beta[k] * qj[k]);
        Ai_j[2][3][k] = (1 - beta[k]) * qj[k];
        Ai_j[0][1][k] = -Ai_j[1][0][k];
        Ai_j[1][1][k] = Ai_j[0][0][k];
        Ai_j[0][2][k] = 0.0;
        Ai_j[0][3][k] = 0.0;
        Ai_j[1][2][k] = 0.0;
        Ai_j[1][3][k] = 0.0;
        Ai_j[2][0][k] = 0.0;
        Ai_j[2][1][k] = 0.0;
        Ai_j[2][2][k] = 1.0;
        Ai_j[3][0][k] = 0.0;
        Ai_j[3][1][k] = 0.0;
        Ai_j[3][2][k] = 0.0;
        Ai_j[3][3][k] = 1.0;
    }
    return 0;
}

int Mechatronic::AijMatrix(
                int N,
                ArrayMatrix<Matrix2D<double>> Aii_,
                ArrayMatrix<Matrix2D<double>> Ai_j,
                ArrayMatrix<Matrix2D<double>> Aij)
{
    for (int k = 0; k < N; k++)
        Mechatronic::MultiMatrix(k, Aii_, k, Ai_j, k, Aij);

    return 0;
}

int Mechatronic::MultiMatrix(
                int kL,
                ArrayMatrix<Matrix2D<double>> ML,
                int kR,
                ArrayMatrix<Matrix2D<double>> MR,
                int k,
                ArrayMatrix<Matrix2D<double>> M)
{
    for (int i = 0; i < 4; i++ )
        for (int j = 0; j < 4; j++ )
            for (int l = 0; l < 4; l++ )
                M[i][j][k] = M[i][j][k] + ML[i][l][kL] * MR[l][j][kR];

    return 0;
}

int Mechatronic::AoiMatrix(
                int N,
                ArrayMatrix<Matrix2D<double>> Aij,
                ArrayMatrix<Matrix2D<double>> Aoi)
{

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            Aoi[i][j][0] = Aij[i][j][0];

    for (int k = 1; k < N; k++)
        Mechatronic::MultiMatrix(k - 1, Aoi, k, Aij, k, Aoi);

    return 0;
}

int Mechatronic::DMatrixAry(
                int N,
                Vector<int> beta,
                ArrayMatrix<Matrix2D<double>> D)
{
    if (N >= Nmax)
        return 1;

    for (int k = 0; k < N; k++)
    {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                D[i][j][k] = 0.0;
        if (beta[k] == 1)
        {
            D[0][1][k] = -1.0;
            D[1][0][k] = 1.0;
        }
        else if (beta[k] == 0)
            D[2][3][k] = 1.0;
        else
            return 2;
    }

    return 0;
}

int Mechatronic::DAijMatrixAry(
                    int N,
                    ArrayMatrix<Matrix2D<double>> D,
                    ArrayMatrix<Matrix2D<double>> Aij,
                    ArrayMatrix<Matrix2D<double>> DAij)
{
    for (int k = 0; k < N; k++)
        MultiMatrix(k, Aij, k, D, k, DAij);

    return 0;
}

int Mechatronic::DiAokMatrix(
                int k,
                int i,
                ArrayMatrix<Matrix2D<double>> Aok,
                ArrayMatrix<Matrix2D<double>> Aij,
                ArrayMatrix<Matrix2D<double>> DAij,
                Matrix2D<double> DiAok)
{
    Matrix2D<double> M(4, 4);
    if (i > k)
        return 0;

    if (k == 1 && i == k) // и
    {
        for (int l1 = 0; l1 < 4; l1++)
            for (int l2 = 0; l2 < 4; l2++)
                DiAok[l1][l2] = DAij[l1][l2][i - 1];
        return 0;
    }
    if (i == 1)
        for (int l1 = 0; l1 < 4; l1++)
            for (int l2 = 0; l2 < 4; l2++)
                DiAok[l1][l2] = DAij[l1][l2][i - 1];
    else
        for (int l1 = 0; l1 < 4; l1++)
            for (int l2 = 0; l2 < 4; l2++)
            {
                DiAok[l1][l2] = 0.;
                for (int l3 = 0; l3 < 4; l3++)
                    DiAok[l1][l2] =
                        DiAok[l1][l2] + Aok[l1][l3][i - 2] * DAij[l3][l2][i - 1];
            }

    for (int j = i; j < k; j++)
    {
        for (int l1 = 0; l1 < 4; l1++)
            for (int l2 = 0; l2 < 4; l2++)
            {
                M[l1][l2] = 0.;
                for (int l3 = 0; l3 < 4; l3++)
                    M[l1][l2] =
                        M[l1][l2] + DiAok[l1][l3] * Aij[l3][l2][j];
            }
        DiAok = M;
    }
    return 0;
}

int Mechatronic::DiiAokMatrix(
                int k,
                int i,
                ArrayMatrix<Matrix2D<double>> Aok,
                ArrayMatrix<Matrix2D<double>> Aij,
                ArrayMatrix<Matrix2D<double>> DAij,
                ArrayMatrix<Matrix2D<double>> D,
                Matrix2D<double> DiiAok)
{

    Matrix2D<double> M(4, 4);
    if (i > k)
    {
        DiiAok = Matrix2D<double>(4, 4);
        return 0;
    }

    if (k == 1 && i == k)
    {
        DiiAok = Matrix2D<double>(4, 4);

        for (int l1 = 0; l1 < 4; l1++)
            for (int l2 = 0; l2 < 4; l2++)
                for (int l3 = 0; l3 < 4; l3++)
                    DiiAok[l1][l2] += DAij[l1][l3][i - 1] * D[l3][l2][i - 1];

        return 0;
    }

    for (int l1 = 0; l1 < 4; l1++)
        for (int l2 = 0; l2 < 4; l2++)
            for (int l3 = 0; l3 < 4; l3++)
                M[l1][l2] += Aok[l1][l3][i - 2] * DAij[l3][l2][i - 1];

    DiiAok = Matrix2D<double>(4, 4);

    for (int l1 = 0; l1 < 4; l1++)
        for (int l2 = 0; l2 < 4; l2++)
            for (int l3 = 0; l3 < 4; l3++)
                DiiAok[l1][l2] += M[l1][l3] * D[l3][l2][i - 1];

    M = Matrix2D<double>(4, 4);

    for (int j = i + 1; j <= k; j++)
    {
        for (int l1 = 0; l1 < 4; l1++)
            for (int l2 = 0; l2 < 4; l2++)
                for (int l3 = 0; l3 < 4; l3++)
                    M[l1][l2] += DiiAok[l1][l3] * Aij[l3][l2][j - 1];

        DiiAok = M;
    }

    return 0;
}

int Mechatronic::DijAokMatrix(int k,
                 int i,
                 int j,
                 ArrayMatrix<Matrix2D<double>> Aok,
                 ArrayMatrix<Matrix2D<double>> Aij,
                 ArrayMatrix<Matrix2D<double>> DAij,
                 ArrayMatrix<Matrix2D<double>> D,
                 Matrix2D<double> DijAok)
{
    Matrix2D<double> M(4, 4);

    if (i > k || j > k)
    {
        DijAok = Matrix2D<double>(4, 4);
        return 0;
    }

    if (i == j)
    {
        DiiAokMatrix(k, i, Aok, Aij, DAij, D, DijAok);
        return 0;
    }

    if (i > j)
    {
        int temp = i;
        i = j;
        j = temp;
    }

    DijAok = Matrix2D<double>(4, 4);

    for (int l1 = 0; l1 < 4; l1++)
        for (int l2 = 0; l2 < 4; l2++)
            for (int l3 = 0; l3 < 4; l3++)
                DijAok[l1][l2] += Aok[l1][l3][i - 2] * DAij[l3][l2][i - 1];

    for (int l = i + 1; l < j; l++)
    {
        for (int l1 = 0; l1 < 4; l1++)
            for (int l2 = 0; l2 < 4; l2++)
                for (int l3 = 0; l3 < 4; l3++)
                    M[l1][l2] += DijAok[l1][l3] * Aij[l3][l2][l - 1];
        DijAok = M;
    }

    M = Matrix2D<double>(4, 4);

    for (int l1 = 0; l1 < 4; l1++)
        for (int l2 = 0; l2 < 4; l2++)
            for (int l3 = 0; l3 < 4; l3++)
                M[l1][l2] += DijAok[l1][l3] * DAij[l3][l2][j - 1];

    DijAok = Matrix2D<double>(4, 4);

    for (int l = j + 1; l <= k; l++)
    {
        for (int l1 = 0; l1 < 4; l1++)
            for (int l2 = 0; l2 < 4; l2++)
                for (int l3 = 0; l3 < 4; l3++)
                    DijAok[l1][l2] = DijAok[l1][l2] + M[l1][l3] * Aij[l3][l2][l - 1];

        if (l < k)
            M = DijAok;
    }
    return 0;
}

int Mechatronic::lin_speed(
                int k,
                Vector<double> l_coord,
                Vector<double> dq,
                ArrayMatrix<Matrix2D<double>> Aok,
                ArrayMatrix<Matrix2D<double>> Aij,
                ArrayMatrix<Matrix2D<double>> DAij,
                Vector<double> l_speed)
{

    Matrix2D<double> M1(4, 4), M2(4, 4), DiAok(4, 4);

    if (k > Nmax)
        return 1;

    for (int i = 1; i <= k; i++)
    {
        Mechatronic::DiAokMatrix(k, i, Aok, Aij, DAij, DiAok);

        M2 = dq[i] * DiAok;
        M1 = M1 + M2;
    }

    l_speed = l_coord * M1;

    return 0;
}

int Mechatronic::lin_acceler(int k,
                Vector<double> l_coord,
                Vector<double> dq,
                Vector<double> ddq,
                ArrayMatrix<Matrix2D<double>> Aok,
                ArrayMatrix<Matrix2D<double>> Aij,
                ArrayMatrix<Matrix2D<double>> D,
                ArrayMatrix<Matrix2D<double>> DAij,
                Vector<double> l_acceler)
{
    int i, j, l1, l2;
    double dqidqj;

    Matrix2D<double> M1(4, 4), M2(4, 4), M3(4, 4);

    if (k > Nmax)
        return 1;

    for (i = 1; i <= k; i++)
    {
        DiAokMatrix(k, i, Aok, Aij, DAij, M2);

        M3 = ddq[i] * M2;
        M1 = M1 + M3;

        for (j = 1; j <= k; j++)
        {
            DijAokMatrix(k, i, j, Aok, Aij, DAij, D, M2);
            dqidqj = dq[i] * dq[j];

            M3 = dqidqj * M2;
            M1 = M1 + M2;
        }
    }
    l_acceler = l_coord * M1;
    return 0;
}

int Mechatronic::angul_acceler(int k,
                  Vector<int> beta,
                  Vector<double> dq,
                  Vector<double> ddq,
                  ArrayMatrix<Matrix2D<double>> Aok,
                  ArrayMatrix<Matrix2D<double>> Aij,
                  ArrayMatrix<Matrix2D<double>> DAij,
                  Vector<double> a_acceler)
{
    Vector<double> w1(4), w2(4), w3(4), w4(4);
    Matrix2D<double> M1(4, 4), M2(4, 4), M3(4, 4);

    if (k > Nmax)
        return 1;

    for (int l = 1; l <= k; l++)
    {
        w1[2] = beta[l - 1] * dq[l - 1];
        w2[2] = beta[l - 1] * ddq[l - 1];

        M1 = Matrix2D<double>(4, 4);
        for (int i = 1; i <= l; i++)
        {
            DiAokMatrix(l, i, Aok, Aij, DAij, M2);
            M3 = dq[i] * M2;
            M1 = M1 + M3;
        }

        w3 = w1 * M1;

        for (int l1 = 0; l1 < 4; l1++)
            for (int l2 = 0; l2 < 4; l2++)
                M2[l1][l2] = Aok[l1][l2][l];

        w4 = w2 * M2;

        a_acceler = a_acceler + w3 + w4;
    }

    return 0;
}