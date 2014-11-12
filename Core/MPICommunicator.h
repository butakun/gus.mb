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

protected:

private:
    typedef std::map<void*, int> TagMap;

    std::ofstream mConsole;
    int mNextTag;
    TagMap mTagMap;
};

#endif // INCLUDED_MPI_COMMUNICATOR_H__

