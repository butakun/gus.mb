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
// $Id: DummyCommunicator.h 89 2011-01-28 14:16:26Z kato $
#ifndef INCLUDED_DUMMY_COMMUNICATOR_H__
#define INCLUDED_DUMMY_COMMUNICATOR_H__

#include "Communicator.h"

class DummyCommunicator : public Communicator
{
public:
    DummyCommunicator(int* argc, char*** argv) {}
    virtual ~DummyCommunicator() {}

    virtual int Size() const { return 1; }
    virtual int MyRank() const { return 0; }
    virtual void Barrier() const {}
    virtual void Finalize() const {}
    virtual std::ostream& Console() { return std::cout; }

    virtual void Broadcast(int* buffer, int count, int root) const {}

protected:

private:
};

#endif // INCLUDED_DUMMY_COMMUNICATOR_H__

