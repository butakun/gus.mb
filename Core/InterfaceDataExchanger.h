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
// $Id: InterfaceDataExchanger.h 245 2012-06-07 08:25:14Z kato $
#ifndef INCLUDED_INTERFACE_DATA_EXCHANGER_H__
#define INCLUDED_INTERFACE_DATA_EXCHANGER_H__

#include "Structured.h"
#include "AbuttingInterface.h"
#include "BlockPatch.h"
#include <map>
#include <string>
#include <mpi.h>

class Model;

class InterfaceDataExchanger
{
public:
    InterfaceDataExchanger(const AbuttingInterface& interface, Model* model, const char* dataName);
    virtual ~InterfaceDataExchanger();

    void Exchange();

protected:
    void SendLocal(const BlockPatch& self, const BlockPatch& donor, const Structured<double>& Data);
    void SendRemote(const BlockPatch& self, const BlockPatch& donor, const Structured<double>& Data);
    void RecvLocal(const BlockPatch& donor);
    void RecvRemote(const BlockPatch& donor);
    void Finish(); // Wait until all pending send/recv requests are completed.

    void PrepareDataForSend(const BlockPatch& self, Structured<double>& sendBuf, const Structured<double>& data);
    void FinalizeDataAfterRecv();

private:
    friend class InterfaceDataAdaptor;

    const AbuttingInterface& mInterface;

    Model* mModel;

    std::string mDataName;

    // the data are indexed by an id unique to the patch (BlockPatch::UniqueID())
    typedef std::map<int, Structured<double>*> UniqueIDToDataMap;
    UniqueIDToDataMap mSendStore;
    UniqueIDToDataMap mRecvStore;

    typedef std::map<int, MPI_Request> ReqMap;
    ReqMap mSendReqs;
    ReqMap mRecvReqs;
};

#endif // INCLUDED_INTERFACE_DATA_EXCHANGER_H__

