// $Id: InterfaceDataExchanger.cpp 315 2014-01-31 10:01:26Z kato $

#include "Communicator.h"
#include "InterfaceDataExchanger.h"
#include "AbuttingInterface.h"
#include "Model.h"
#include "Block.h"
#include "Roster.h"
#include <string>
#include <mpi.h>

class InterfaceDataAdaptor : public InterfaceDataAdaptorBase
{
public:
    InterfaceDataAdaptor(InterfaceDataExchanger* idx, const char* dataName)
    : mIDX(idx), mDataName(dataName)
    {}

    virtual Structured<double>& GetBlockData(int blockID) const
    {
        return Roster::GetInstance()->GetBlockData(blockID, mDataName.c_str());
    }    

    virtual Structured<double>& GetBlockPatchData(int blockPatchUniqueID) const
    {
        //Structured<double>* data = mIDX->mRecvStore.find(blockPatchUniqueID)->second;
        InterfaceDataExchanger::UniqueIDToDataMap::const_iterator i = mIDX->mRecvStore.find(blockPatchUniqueID);
        assert(i != mIDX->mRecvStore.end());
        Structured<double>* data = i->second;
        return *data;
    }

protected:

private:
    InterfaceDataExchanger* mIDX;
    std::string mDataName;
};

InterfaceDataExchanger::InterfaceDataExchanger(const AbuttingInterface& interface, Model* model, const char* dataName)
:   mInterface(interface), mModel(model), mDataName(dataName)
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    size_t numSelfPatches = mInterface.SelfBlockPatches().size();
    size_t numDonorPatches = mInterface.DonorBlockPatches().size();

    for (size_t i = 0; i < numSelfPatches; ++i)
    {
        const BlockPatch& bp = mInterface.SelfBlockPatches()[i];
        if (Roster::GetInstance()->GetRankOf(bp.BlockID()) == Communicator::GetInstance()->MyRank())
        {
            // Made sure this patch is indeed in the local blocks.
            const Structured<double>& data = Roster::GetInstance()->GetBlockData(bp.BlockID(), mDataName.c_str());
            int uniqueID = bp.UniqueID();
            assert(mSendStore.find(uniqueID) == mSendStore.end());
            mSendStore[uniqueID] = new Structured<double>(data.DOF(), bp.CellRange());
        }
    }

    size_t dof = mSendStore.begin()->second->DOF();
    for (size_t i = 0; i < numDonorPatches; ++i)
    {
        const BlockPatch& bp = mInterface.DonorBlockPatches()[i];
        assert(mRecvStore.find(bp.UniqueID()) == mRecvStore.end());
        mRecvStore[bp.UniqueID()] = new Structured<double>(dof, bp.CellRange());
    }
}

InterfaceDataExchanger::~InterfaceDataExchanger()
{
    for (std::map<int, Structured<double>*>::iterator i = mSendStore.begin();
        i != mSendStore.end(); ++i)
    {
        Structured<double>* s = i->second;
        delete[] s->Data;
        delete s;
    }

    for (std::map<int, Structured<double>*>::iterator i = mRecvStore.begin();
        i != mRecvStore.end(); ++i)
    {
        Structured<double>* s = i->second;
        delete[] s->Data;
        delete s;
    }
}

void
InterfaceDataExchanger::Exchange()
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    assert(mSendReqs.size() == 0);
    assert(mRecvReqs.size() == 0);

    // Send
    for (AbuttingInterface::BlockPatches::const_iterator i = mInterface.SelfBlockPatches().begin();
        i != mInterface.SelfBlockPatches().end(); ++i)
    {
        const BlockPatch& self = *i;
        Block* block = dynamic_cast<Block*>(Roster::GetInstance()->GetBlock(self.BlockID()));
        assert(block != NULL);
        const Structured<double>& sendData = Roster::GetInstance()->GetBlockData(*block, mDataName.c_str());
        //LOG << "IDF::Exchange: preparing to send named data " << mDataName << " from Block " << block->ID() << " buf = " << sendData.Data << std::endl;
        for (AbuttingInterface::BlockPatches::const_iterator j = mInterface.DonorBlockPatches().begin();
            j != mInterface.DonorBlockPatches().end(); ++j)
        {
            const BlockPatch& donor = *j;
            int donorRank = Roster::GetInstance()->GetRankOf(donor.BlockID());
            if (donorRank == Communicator::GetInstance()->MyRank())
            {
                SendLocal(self, donor, sendData);
            }
            else
            {
                SendRemote(self, donor, sendData);
            }
        }
    }

    // Receive
    for (AbuttingInterface::BlockPatches::const_iterator i = mInterface.DonorBlockPatches().begin();
        i != mInterface.DonorBlockPatches().end(); ++i)
    {
        const BlockPatch& donor = *i;
        //LOG << "IDF::Exchange: receiving data from donor patch " << donor << std::endl;
        VirtualBlock* block = Roster::GetInstance()->GetBlock(donor.BlockID());
        assert(block != NULL);
        int donorRank = Roster::GetInstance()->GetRankOf(donor.BlockID());
        if (donorRank == Communicator::GetInstance()->MyRank())
        {
            RecvLocal(donor);
        }
        else
        {
            RecvRemote(donor);
        }
    }

    // Wait for the completion of all the pending send/recvs.
    Finish();

    // Then map the solution on the interface. FIXME: does this have to be here?
    InterfaceDataAdaptor ida(this, mDataName.c_str());
    mInterface.MapData(*mModel, &ida);
}

void
InterfaceDataExchanger::PrepareDataForSend(const BlockPatch& self, Structured<double>& sendBuf, const Structured<double>& data)
{
    // Copy the data into a consectutive memory and transform them to the inertial frame of reference
    Block* block = dynamic_cast<Block*>(Roster::GetInstance()->GetBlock(self.BlockID()));
    assert(block != NULL);

    const IndexRange& cr = self.CellRange();
    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* src = data(i, j, k);
                double* dst = sendBuf(i, j, k);
                mModel->FromLocalToGlobal(dst, src, *block, IndexIJK(i, j, k));
            }
        }
    }
}

void
InterfaceDataExchanger::SendLocal(const BlockPatch& self, const BlockPatch& donor, const Structured<double>& data)
{
}

void
InterfaceDataExchanger::SendRemote(const BlockPatch& self, const BlockPatch& donor, const Structured<double>& data)
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    Structured<double>& sendBuf = *mSendStore[self.UniqueID()];

    PrepareDataForSend(self, sendBuf, data);

    int donorRank = Roster::GetInstance()->GetRankOf(donor.BlockID());
    int tag = self.UniqueID();
    MPI_Request req;
    int err;

    //LOG << "IDF::SendRemote: sending (async) to rank " << donorRank << std::endl;
    err = MPI_Isend(
        sendBuf.Data,
        sendBuf.GetRange().Count() * sendBuf.DOF(),
        MPI_DOUBLE,
        donorRank,
        tag,
        MPI_COMM_WORLD,
        &req
        );
    assert(err == MPI_SUCCESS);

    //LOG << "IDF::SendRemote: async send request = " << req << std::endl;
    mSendReqs[tag] = req;
}

void
InterfaceDataExchanger::RecvLocal(const BlockPatch& donor)
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    const Block* block = dynamic_cast<const Block*>(Roster::GetInstance()->GetBlock(donor.BlockID()));
    assert(block != NULL);

    const Structured<double>& srcData = Roster::GetInstance()->GetBlockData(block->ID(), mDataName.c_str());

    Structured<double>& recvBuf = *mRecvStore[donor.UniqueID()];
    //LOG << "IDX::RecvLocal: donor = " << donor << ", CellRange = " << donor.CellRange() << std::endl;
    //LOG << "IDX::RecvLocal: recvBuf = " << recvBuf.GetRange() << ", srcData.CellRange = " << srcData.GetRange() << std::endl;

    int dof = recvBuf.DOF();
    const IndexRange& cr = donor.CellRange();
    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* src = srcData(i, j, k);
                double* dst = recvBuf(i, j, k);
                mModel->FromLocalToGlobal(dst, src, *block, IndexIJK(i, j, k));
            }
        }
    }
}

void
InterfaceDataExchanger::RecvRemote(const BlockPatch& donor)
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    Structured<double>& recvBuf = *mRecvStore[donor.UniqueID()];
    //LOG << "IDF::RecvRemote: receiving data into temporary storage " << recvBuf.Data << std::endl;

    int donorRank = Roster::GetInstance()->GetRankOf(donor.BlockID());
    int tag = donor.UniqueID();
    MPI_Request req;
    int err;

    //LOG << "IDF::RecvRemote: receiving (async) from rank " << donorRank << std::endl;
    err = MPI_Irecv(
        recvBuf.Data,
        recvBuf.GetRange().Count() * recvBuf.DOF(),
        MPI_DOUBLE,
        donorRank,
        tag,
        MPI_COMM_WORLD,
        &req
        );
    assert(err == MPI_SUCCESS);

    //LOG << "IDF::RecvRemote: async recv request = " << req << std::endl;
    mRecvReqs[tag] = req;
}

void
InterfaceDataExchanger::Finish()
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    int count = mSendReqs.size() + mRecvReqs.size();
    MPI_Request* reqs = new MPI_Request[count];
    MPI_Status* statuses = new MPI_Status[count];

    MPI_Request* r = &reqs[0];
    for (ReqMap::const_iterator i = mSendReqs.begin(); i != mSendReqs.end(); ++i)
    {
        *r++ = i->second;
    }
    for (ReqMap::const_iterator i = mRecvReqs.begin(); i != mRecvReqs.end(); ++i)
    {
        *r++ = i->second;
    }

    //LOG << "IDF::Finish: waiting for the pending send/recv requests" << std::endl;
    int err;
    err = MPI_Waitall(count, reqs, statuses);
    assert(err == MPI_SUCCESS);
    //LOG << "IDF::Finish: all send/recv requests completed" << std::endl;

    delete[] reqs;
    delete[] statuses;

    mSendReqs.clear();
    mRecvReqs.clear();

    FinalizeDataAfterRecv();
}

void
InterfaceDataExchanger::FinalizeDataAfterRecv()
{
}

