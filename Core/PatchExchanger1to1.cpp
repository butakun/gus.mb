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

#include "Communicator.h"
#include "PatchExchanger1to1.h"
#include "Connectivity1to1.h"
#include "Roster.h"
#include "StructuredDataExchanger.h"

PatchExchanger1to1::PatchExchanger1to1(StructuredDataExchanger* ex, const Connectivity1to1& conn, Structured<double>& data)
:   mEx(ex), mConn(conn), mData(data),
    mDataToSend(data.DOF(), conn.CellRangeToSend()),
    mDataToRecv(data.DOF(), conn.DonorCellRangeToRecv().Canonical())
{
    assert(mConn.CellRangeToSend().IsCanonical());
}

PatchExchanger1to1::~PatchExchanger1to1()
{
    delete[] mDataToSend.Data;
    delete[] mDataToRecv.Data;
}

void
PatchExchanger1to1::Start()
{
    int myRank = Communicator::GetInstance()->MyRank();
    int donorRank = Roster::GetInstance()->GetRankOf(mConn.DonorBlockID());
    if (donorRank == myRank)
    {
        StartLocal();
    }
    else
    {
        StartRemote(myRank, donorRank);
    }
}

void
PatchExchanger1to1::Finish()
{
    int myRank = Communicator::GetInstance()->MyRank();
    int donorRank = Roster::GetInstance()->GetRankOf(mConn.DonorBlockID());
    if (donorRank == myRank)
    {
        FinishLocal();
    }
    else
    {
        FinishRemote(myRank, donorRank);
    }
}

void
PatchExchanger1to1::StartLocal()
{
}

void
PatchExchanger1to1::StartRemote(int myRank, int donorRank)
{
    IndexRange crToSend = mConn.CellRangeToSend();

#if 0
    IndexRange crToRecv = mConn.CellRangeToRecv();
    std::ostream& LOG = Communicator::GetInstance()->Console();
    LOG << "PatchExchanger1to1::StartRemote: Block " << mConn.BlockID() << " <-> DonorBlock " << mConn.DonorBlockID() << ", crToSend = " << crToSend << ", crToRecv = " << crToRecv << std::endl;
#endif

    for (int k = crToSend.Start.K; k <= crToSend.End.K; ++k)
    {
        for (int j = crToSend.Start.J; j <= crToSend.End.J; ++j)
        {
            for (int i = crToSend.Start.I; i <= crToSend.End.I; ++i)
            {
                double* self = mData(i, j, k);
                double* sendbuf = mDataToSend(i, j, k);
                for (int l = 0; l < mData.DOF(); ++l)
                {
                    sendbuf[l] = self[l];
                }
            }
        }
    }

    int err;
    err = MPI_Isend(
        mDataToSend.Data,
        mDataToSend.GetRange().Count() * mDataToSend.DOF(),
        MPI_DOUBLE,
        donorRank, // dest
        mConn.Tag(), // tag
        MPI_COMM_WORLD,
        &mSendReq
        );
    assert(err == MPI_SUCCESS);
    //LOG << "PatchExchanger1to1::StartRemote: Block " << mConn.BlockID() << ", Send request = " << mSendReq << std::endl;

    err = MPI_Irecv(
        mDataToRecv.Data,
        mDataToRecv.GetRange().Count() * mDataToRecv.DOF(),
        MPI_DOUBLE,
        donorRank, // src
        mConn.Tag(),
        MPI_COMM_WORLD,
        &mRecvReq
        );
    assert(err == MPI_SUCCESS);
    //LOG << "PatchExchanger1to1::StartRemote: Block " << mConn.BlockID() << ", Recv request = " << mRecvReq << std::endl;
}

void
PatchExchanger1to1::FinishLocal()
{
    Structured<double>& dataDonor = *mEx->BlockData(mConn.DonorBlockID());

    IndexRange crToRecv = mConn.CellRangeToRecv();

#if 0
    IndexRange crToSend = mConn.CellRangeToSend();
    std::ostream& LOG = Communicator::GetInstance()->Console();
    LOG << "PatchExchanger1to1::FinishLocal:" << mConn.Tag() << " Block " << mConn.BlockID() << " <-> DonorBlock " << mConn.DonorBlockID() << ", crToSend = " << crToSend << ", crToRecv = " << crToRecv << std::endl;
#endif

    assert(crToRecv.IsCanonical());

    for (int k = crToRecv.Start.K; k <= crToRecv.End.K; ++k)
    {
        for (int j = crToRecv.Start.J; j <= crToRecv.End.J; ++j)
        {
            for (int i = crToRecv.Start.I; i <= crToRecv.End.I; ++i)
            {
                IndexIJK iSelf(i, j, k);
                IndexIJK iDonor = mConn.DonorCellIndex(iSelf);
                double* self = mData(iSelf);
                double* donor = dataDonor(iDonor);
                for (int l = 0; l < mData.DOF(); ++l)
                {
                    self[l] = donor[l];
                }
            }
        }
    }
}

std::string mpi_error_string(int errcode)
{
    char buf[MPI_MAX_ERROR_STRING];
    int len;
    int err;
    err = MPI_Error_string(errcode, buf, &len);
    assert(err == MPI_SUCCESS);
    return std::string(buf);
}

void
PatchExchanger1to1::FinishRemote(int myRank, int donorRank)
{
    int err;
    MPI_Status status;

    err = MPI_Wait(&mSendReq, &status);
    assert(err == MPI_SUCCESS);

    err = MPI_Wait(&mRecvReq, &status);
    assert(err == MPI_SUCCESS);

    IndexRange crToRecv = mConn.CellRangeToRecv();
    std::ostream& LOG = Communicator::GetInstance()->Console();
    //LOG << "PatchExchanger1to1::FinishRemote:" << mConn.Tag() << " Block " << mConn.BlockID() << " selfCellRange " << crToRecv << std::endl;
    for (int k = crToRecv.Start.K; k <= crToRecv.End.K; ++k)
    {
        for (int j = crToRecv.Start.J; j <= crToRecv.End.J; ++j)
        {
            for (int i = crToRecv.Start.I; i <= crToRecv.End.I; ++i)
            {
                IndexIJK selfIndex(i, j, k);
                IndexIJK donorIndex = mConn.DonorCellIndex(selfIndex);
                double* self = mData(selfIndex);
                double* recvbuf = mDataToRecv(donorIndex);
                for (int l = 0; l < mData.DOF(); ++l)
                {
                    self[l] = recvbuf[l];
                }
            }
        }
    }
}

