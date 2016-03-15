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
// $Id: Communicator.cpp 27 2010-05-23 15:51:53Z kato $

#include "Communicator.h"
#include "DummyCommunicator.h"
#include "MPICommunicator.h"
#include <cstdlib>
#include <cassert>

Communicator* Communicator::Instance = NULL;

void Communicator::Initialize(int* argc, char*** argv)
{
    assert(Communicator::Instance == NULL);
    //Communicator::Instance = new DummyCommunicator(argc, argv);
    Communicator::Instance = new MPICommunicator(argc, argv);
}

Communicator* Communicator::GetInstance()
{
    assert(Communicator::Instance != NULL);
    return Communicator::Instance;
}

