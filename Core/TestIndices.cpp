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
// $Id: TestIndices.cpp 73 2010-09-19 15:24:30Z kato $

#include "Connectivity1to1.h"
#include "IndexUtils.h"
#include <iostream>

int main()
{
    int t[3] = {-1, 3, -2};
    IndexTransform transform(t, IndexIJK(0, 0, 0), IndexIJK(0, 0, 0));
    std::cout << transform << std::endl;
    std::cout << transform.CellIndexTransform() << std::endl;

    std::cout << std::endl;
    IndexRange meshRange(18, 0, 0, 24, 10, 0);
    t[0] = 3; t[1] = 2; t[2] = -1;
    IndexTransform transform2(t, meshRange.Start, IndexIJK(0, 0, 0));
    std::cout << transform2 << std::endl;
    std::cout << transform2.CellIndexTransform() << std::endl;

    IndexRange cellRangeToRecv = IndexUtils::FromMeshRangeToCellRange(meshRange, IndexUtils::Opposite(K), 2);

    IndexTransform transform3(transform2.CellIndexTransform());
    IndexRange donorCellRangeToRecv = transform3.Apply(cellRangeToRecv).Canonical();
    std::cout << cellRangeToRecv << ", " << donorCellRangeToRecv << std::endl;

    IndexIterator itor(IndexRange(0, 0, 0, 5, 5, 2));
    while (!itor.IsEnd())
    {
        IndexIJK i = itor.Index();
        std::cout << i << ' ';
        itor.Advance();
    }
    std::cout << std::endl;

    return 0;
}
