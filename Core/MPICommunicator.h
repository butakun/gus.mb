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
// $Id: MPICommunicator.h 225 2012-05-01 23:18:07Z kato $
#ifndef INCLUDED_MPI_COMMUNICATOR_H__
#define INCLUDED_MPI_COMMUNICATOR_H__

#include "Communicator.h"
#include <fstream>
#include <map>

class MPICommunicator : public Communicator
{
public:
    MPICommunicator(int* argc, char*** argv);
    virtual ~MPICommunicator();

    virtual int Size() const;
    virtual int MyRank() const;
    virtual void Barrier() const;
    virtual void Finalize() const;
    virtual std::ostream& Console();
    virtual int ReserveTag(void* key);

    virtual void Broadcast(int* buffer, int count, int root) const;
    virtual bool Any(bool flag) const;

    virtual int ReduceSum(int value, int root) const;
    virtual int AllReduceSum(int value) const;

protected:

private:
    typedef std::map<void*, int> TagMap;

    std::ofstream mConsole;
    int mNextTag;
    TagMap mTagMap;
};

#endif // INCLUDED_MPI_COMMUNICATOR_H__

