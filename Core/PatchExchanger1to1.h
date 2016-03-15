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
// $Id: PatchExchanger1to1.h 39 2010-06-21 02:42:48Z kato $
#ifndef INCLUDED_PATCH_EXCHANGER_1to1_H__
#define INCLUDED_PATCH_EXCHANGER_1to1_H__

#include "PatchExchanger.h"
#include "Structured.h"
#include <mpi.h> // FIXME: MPI-specific types should not be here.

class Connectivity1to1;
class StructuredDataExchanger;

class PatchExchanger1to1 : public PatchExchanger
{
public:
    PatchExchanger1to1(StructuredDataExchanger* ex, const Connectivity1to1& conn, Structured<double>& data);
    virtual ~PatchExchanger1to1();

    virtual void Start();
    virtual void Finish();

protected:
    void StartLocal();
    void StartRemote(int myRank, int donorRank);
    void FinishLocal();
    void FinishRemote(int myRank, int donorRank);

private:
    StructuredDataExchanger* mEx;
    const Connectivity1to1& mConn;
    Structured<double>& mData;
    Structured<double> mDataToSend;
    Structured<double> mDataToRecv;

    MPI_Request mSendReq;
    MPI_Request mRecvReq;
};

#endif // INCLUDED_PATCH_EXCHANGER_1to1_H__

