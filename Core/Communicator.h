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

    virtual std::ostream& Console() = 0;

    virtual int ReserveTag(void* key) = 0;

protected:

private:
    static Communicator* Instance;
};

#endif // INCLUDED_COMMUNICATOR_H__

