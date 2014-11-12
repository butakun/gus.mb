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

