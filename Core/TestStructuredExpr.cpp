/*
    gus.mb, an open source flow solver.
    Copyright (C) 2016 Hiromasa Kato <hiromasa at gmail.com>

    This file is part of gus.mb.

    gus.mb is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    gus.mb is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
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

