// $Id: TestStructuredExpr.cpp 122 2011-08-18 16:13:30Z kato $

#include "Structured.h"
#include <iostream>

int main()
{
    IndexRange range(0, 0, 0, 1, 1, 0);
    Structured<double> A(1, range), R(1, range);

    A = 1.0;
    std::cout << "A = " << A << std::endl;

    R = 2.0 + A;
    std::cout << "2.0 + A = " << R << std::endl;

    R = 1.2 + 2.0 + A;
    std::cout << "1.2 + 2.0 + A = " << R << std::endl;

    R = A + 2.0;
    std::cout << "A + 2.0 = " << R << std::endl;

    R = 1.2 + A + 2.0;
    std::cout << "1.2d + A + 2.0 = " << R << std::endl;

    R = 2.0 * A;
    std::cout << "2.0 * A = " << R << std::endl;

    R = A * 2.0;
    std::cout << "A * 2.0 = " << R << std::endl;

    R = 2.0 * A + 1.2;
    std::cout << "2.0 * A + 1.2 = " << R << std::endl;

    R = A / 0.3;
    std::cout << "A / 0.3 = " << R << std::endl;

    R = 2.0 * A / 0.3;
    std::cout << "2.0 * A / 0.3 = " << R << std::endl;
}

