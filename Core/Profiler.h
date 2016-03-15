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
// $Id: Profiler.h 189 2012-01-19 07:40:08Z kato $
#ifndef INCLUDED_PROFILER_H__
#define INCLUDED_PROFILER_H__

class Profiler
{
public:
    static Profiler* GetInstance();

    virtual ~Profiler();

    double CheckPoint(); // returns the elapsed time since the last checkpoint.

protected:
    Profiler();

private:
    static Profiler* mInstance;

    double mTimeStamp;
};

#endif // INCLUDED_PROFILER_H__

