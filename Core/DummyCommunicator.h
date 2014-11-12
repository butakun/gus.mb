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

