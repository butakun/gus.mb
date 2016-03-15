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
// $Id: Communicator.h 225 2012-05-01 23:18:07Z kato $
#ifndef INCLUDED_COMMUNICATOR_H__
#define INCLUDED_COMMUNICATOR_H__

#include <mpi.h>
#include <iostream>

class Communicator
{
public:
    static void Initialize(int* argc, char*** argv);

    static Communicator* GetInstance();

    virtual ~Communicator() {}

    virtual int Size() const = 0;
    virtual int MyRank() const = 0;
    virtual void Barrier() const = 0;
    virtual void Finalize() const = 0;

    virtual void Broadcast(int* buffer, int count, int root) const = 0;
    //virtual void AllGather() const = 0;
    virtual bool Any(bool flag) const = 0; // AllReduce and LogicalOr operation on boolean flags from all ranks.

    virtual int ReduceSum(int value, int root) const = 0;
    virtual int AllReduceSum(int value) const = 0;

    virtual std::ostream& Console() = 0;

    virtual int ReserveTag(void* key) = 0;

protected:

private:
    static Communicator* Instance;
};

#endif // INCLUDED_COMMUNICATOR_H__

