#include "Matrix33.h"
#include <iostream>
#include <cstdlib>
#include <ctime>

Matrix33 identity_matrix()
{
    Matrix33 m;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (i == j)
                m(i, j) = 1.0;
            else
                m(i, j) = 0.0;
        }
    }

    return m;
}

Matrix33 random_matrix()
{
    Matrix33 m;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            m(i, j) = 2.0 * double(std::rand()) / double(RAND_MAX) - 1.0;
        }
    }

    return m;
}

int main()
{
    std::srand(std::time(NULL));

    Matrix33 A = random_matrix();
    //Matrix33 A = identity_matrix();
    Matrix33 Ainv = A.Inverse();

    Matrix33 I = A * Ainv;

    std::cout << A << std::endl;
    std::cout << Ainv << std::endl;
    std::cout << I << std::endl;
}

