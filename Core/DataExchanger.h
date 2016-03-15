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
// $Id: DataExchanger.h 21 2010-05-17 13:50:02Z kato $
#ifndef INCLUDED_DATA_EXCHANGER_H__
#define INCLUDED_DATA_EXCHANGER_H__

class DataExchanger
{
public:
    virtual ~DataExchanger() {}

    virtual void Start() = 0;
    virtual void Finish() = 0;

protected:

private:
};

#endif // INCLUDED_DATA_EXCHANGER_H__

