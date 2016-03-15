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
// $Id: Connectivity.h 44 2010-07-08 16:15:52Z kato $
#ifndef INCLUDED_CONNECTIVITY_H__
#define INCLUDED_CONNECTIVITY_H__

#include "BCPlanar.h"

class Connectivity : public BCPlanar
{
public:
    Connectivity(const IndexRange& meshRange, Direction direction) : BCPlanar(meshRange, direction) {}
    virtual ~Connectivity() {}

protected:

private:
};

#endif // INCLUDED_CONNECTIVITY_H__

