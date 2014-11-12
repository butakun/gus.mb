// $Id: DomainInfo.h 58 2010-08-20 14:57:19Z kato $
#ifndef INCLUDED_DOMAIN_INFO_H__
#define INCLUDED_DOMAIN_INFO_H__

#include "Structured.h"
#include "Connectivity1to1.h"
#include "Vector3.h"
#include <vector>
#include <functional>
#include <algorithm>
#include <string>

class Connectivity1to1Info
{
public:
    Connectivity1to1Info() {}
    Connectivity1to1Info(const char* name, int zone, const IndexRange& meshRange,
        int donorZone, const IndexRange& donorMeshRange, int* transform)
    :   mName(name), mZone(zone), mMeshRange(meshRange),
        mDonorZone(donorZone), mDonorMeshRange(donorMeshRange),
        mTransform(transform, mMeshRange.Start, mDonorMeshRange.Start),
        mPeriodicity(Connectivity1to1::NONE)
    {
    }

    Connectivity1to1Info(const Connectivity1to1Info& conn)
    :   mName(conn.mName), mZone(conn.mZone), mMeshRange(conn.mMeshRange),
        mDonorZone(conn.mDonorZone), mDonorMeshRange(conn.mDonorMeshRange),
        mTransform(conn.mTransform),
        mPeriodicity(conn.mPeriodicity), mRotCenter(conn.mRotCenter), mRotAngle(conn.mRotAngle)
    {
    }

    const IndexRange& MeshRange() const { return mMeshRange; }
    const IndexTransform& Transform() const { return mTransform; }
    int Zone() const { return mZone; }
    int DonorZone() const { return mDonorZone; }

    void SetPeriodicity(Connectivity1to1::Periodicity periodicity, const Vector3& rotCenter, const Vector3& rotAngle)
    {
        mPeriodicity = periodicity;
        mRotCenter = rotCenter;
        mRotAngle = rotAngle;
    }
    Connectivity1to1::Periodicity GetPeriodicity() const { return mPeriodicity; }
    const Vector3& GetRotationCenter() const { return mRotCenter; }
    const Vector3& GetRotationAngle() const { return mRotAngle; }

    std::ostream& Dump(std::ostream& o) const
    {
        o << "Name: " << mName << std::endl
        << "Zone: " << mZone << std::endl
        << "MeshRange: " << mMeshRange << std::endl
        << "DonorZone: " << mDonorZone << std::endl
        << "DonorMeshRange: " << mDonorMeshRange << std::endl
        << "Transform: " << mTransform;
        if (mPeriodicity != Connectivity1to1::NONE)
        {
            o << std::endl << "Periodicity: center = " << mRotCenter << ", angle = " << mRotAngle;
        }
        return o;
    }

protected:

private:
    std::string mName;
    int mZone;
    IndexRange mMeshRange;
    int mDonorZone;
    IndexRange mDonorMeshRange;
    IndexTransform mTransform;

    Connectivity1to1::Periodicity mPeriodicity;
    Vector3 mRotCenter, mRotAngle;
};

inline
std::ostream& operator << (std::ostream& o, const Connectivity1to1Info& conn)
{
    return conn.Dump(o);
}

typedef std::vector<Connectivity1to1Info> Conn1to1s;

class BCInfo
{
public:
    BCInfo() {}
    BCInfo(const char* name, const char* type, const IndexRange& meshRange)
    : mName(name), mType(type), mMeshRange(meshRange)
    {
    }

    const char* Name() const { return mName.c_str(); }
    const char* Type() const { return mType.c_str(); }
    const IndexRange& MeshRange() const { return mMeshRange; }

    bool IsName(const char* name) const { return mName == name; }

    std::ostream& Dump(std::ostream& o) const
    {
        o << "BC: " << Name() << ", Type " << Type() << ", " << MeshRange();
        return o;
    }

protected:

private:
    std::string mName;
    std::string mType;
    IndexRange mMeshRange;
};

inline
std::ostream& operator << (std::ostream& o, const BCInfo& bci)
{
    return bci.Dump(o);
}

typedef std::vector<BCInfo> BCInfos;

class BlockInfo
{
public:
    BlockInfo(const char* name, int zone, const IndexRange& meshRange)
    :   mName(name), mZone(zone), mMeshRange(meshRange)
    {}

    const char* Name() const { return mName.c_str(); }
    int Zone() const { return mZone; }
    IndexRange MeshRange() const { return mMeshRange; }
    size_t NumBCs() const { return mBCInfos.size(); }
    size_t NumConn1to1s() const { return mConn1to1s.size(); }
    size_t NumConnGens() const { return 0; } // FIXME
    const Vector3& AngularVelocity() const { return mAngularVelocity; }
    Vector3& AngularVelocity() { return mAngularVelocity; }

    const BCInfos& GetBCInfos() const { return mBCInfos; }
    BCInfos& GetBCInfos() { return mBCInfos; }
    const Conn1to1s& GetConn1to1s() const { return mConn1to1s; }
    Conn1to1s& GetConn1to1s() { return mConn1to1s; }

    bool operator == (const char* name) const { return mName == name; }
    bool operator == (int zone) const { return mZone == zone; }

    const BCInfo& FindBCInfoByName(const char* name) const
    {
        BCInfos::const_iterator i;
        i = std::find_if(mBCInfos.begin(), mBCInfos.end(), std::bind2nd(std::mem_fun_ref(&BCInfo::IsName), name));
        assert(i != mBCInfos.end());
        return *i;
    }

    std::ostream& Dump(std::ostream& o) const
    {
        o << "Name: " << Name() << std::endl
        << "Zone: " << Zone() << std::endl
        << "MeshRange: " << MeshRange() << std::endl
        << "# of BCs: " << NumBCs() << std::endl
        << "# of 1-to-1 connectivities: " << NumConn1to1s() << std::endl
        << "# of general connectivities: " << NumConnGens() << std::endl;
        for (size_t i = 0; i < NumBCs(); ++i)
            o << "  " << i + 1 << ":" << mBCInfos[i] << std::endl;
        for (size_t i = 0; i < NumConn1to1s(); ++i)
            o << mConn1to1s[i] << std::endl;
        return o;
    }

protected:

private:
    std::string mName;
    int mZone;
    IndexRange mMeshRange;
    Conn1to1s mConn1to1s;
    BCInfos mBCInfos;
    Vector3 mAngularVelocity;
};

inline
std::ostream& operator << (std::ostream& o, const BlockInfo& bi)
{
    return bi.Dump(o);
}

typedef std::vector<BlockInfo> BlockInfos;

class DomainInfo
{
public:
    DomainInfo() {}

    void AddBlockInfo(const BlockInfo& bi) { mBlocks.push_back(bi); }
    const BlockInfo& FindBlockInfo(const char* name) const
    {
        BlockInfos::const_iterator i = std::find(mBlocks.begin(), mBlocks.end(), name);
        return *i;
    }
    BlockInfo& FindBlockInfo(const char* name)
    {
        BlockInfos::iterator i = std::find(mBlocks.begin(), mBlocks.end(), name);
        return *i;
    }
    const BlockInfo& FindBlockInfo(int zone) const
    {
        BlockInfos::const_iterator i = std::find(mBlocks.begin(), mBlocks.end(), zone);
        return *i;
    }
    BlockInfo& FindBlockInfo(int zone)
    {
        BlockInfos::iterator i = std::find(mBlocks.begin(), mBlocks.end(), zone);
        return *i;
    }

    const BlockInfos& GetBlockInfos() const
    {
        return mBlocks;
    }

protected:

private:
    BlockInfos mBlocks;
};

#endif // INCLUDED_DOMAIN_INFO_H__

