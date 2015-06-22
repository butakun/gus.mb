// $Id: Block.h 255 2012-12-05 07:47:00Z kato $
#ifndef INCLUDED_HIRO_BLOCK_H__
#define INCLUDED_HIRO_BLOCK_H__

#include "VirtualBlock.h"
#include "Structured.h"
#include "Vector3.h"
#include <vector>

class BC;

typedef std::vector<BC*> BCs;

class Block : public VirtualBlock
{
public:
    static Block* New(int id, const IndexRange& meshRange, bool unsteady = false);

    virtual ~Block();

    int GhostLayers() const { return mNGhostLayers; }
    const IndexRange& Dim() const { return mDim; }
    const IndexRange& CellRangeWithGhosts() const { return mDim; }

    Structured<double>& XYZ() { return mXYZ; }
    Structured<double>& U() { return mU; }
    Structured<double>& UT() { return mUT; }
    Structured<double>& MuK() { return mMuK; }
    Structured<double>& TurMuK() { return mTurMuK; }
    Structured<double>& Sxi() { return mSxi; }
    Structured<double>& Seta() { return mSeta; }
    Structured<double>& Szeta() { return mSzeta; }
    Structured<double>& Vol() { return mVol; }
    Structured<double>& Radius() { return mRadius; }
    const Structured<double>& XYZ() const { return mXYZ; }
    const Structured<double>& U() const { return mU; }
    const Structured<double>& UT() const { return mUT; }
    const Structured<double>& MuK() const { return mMuK; }
    const Structured<double>& TurMuK() const { return mTurMuK; }
    const Structured<double>& Sxi() const { return mSxi; }
    const Structured<double>& Seta() const { return mSeta; }
    const Structured<double>& Szeta() const { return mSzeta; }
    const Structured<double>& Vol() const { return mVol; }
    const Structured<double>& Radius() const { return mRadius; }
    const Structured<int>& Mask() const { return mMask; }

    void ComputeMetrics();

    void ComputeTransportProperties();
    void ComputeTransportProperties(Structured<double>& MuK, Structured<double>& U, const Structured<double>& Radius) const;
    void ComputeTurbulentTransportProperties();

    bool CheckNegatives() const { std::vector<IndexIJK> indices; return CheckNegatives(mU, indices); }
    bool CheckNegatives(std::vector<IndexIJK>& indices) const { return CheckNegatives(mU, indices); }
    bool CheckNegatives(const Structured<double>& U, std::vector<IndexIJK>& indices) const;

    BCs& GetBCs() { return mBCs; }
    const BCs& GetBCs() const { return mBCs; }
    void RegisterBC(BC* bcfunctor);
    void ApplyBCs();
    void ApplyBCs(Structured<double> U);
    void ApplyTurbBCs();
    void ApplyViscousWallBC();
    void FinalizeConnectivities(Structured<double>& U);

    bool IsZone(int id) const { return ID() == id; }
    bool operator == (int id) const { return ID() == id; }

    void FillCornerGhosts();
    void FillCornerGhosts(Structured<double>& U) const;

    Structured<double>& UStorage(int i);
    Structured<double>& UTStorage(int i);
    void ShiftTime();

protected:
    Block(int id, const IndexRange& meshRange, int temporalstore = 1, int nGhosts = 2);

private:
    int mNGhostLayers;      // # of ghost cell layers. (rind cells)
    IndexRange mDim;        // actual dimensions of the matrices.

    Structured<double> mXYZ;
    Structured<double> mU, mUT;
    Structured<double> mU2, mU3, mUT2, mUT3;
    Structured<double> mSxi, mSeta, mSzeta, mVol, mRadius;
    Structured<double> mMuK;
    Structured<double> mTurMuK;
    Structured<double> mDistance;
    Structured<int> mMask;

    BCs mBCs;
    Structured<BC*> mBCIMin, mBCIMax, mBCJMin, mBCJMax, mBCKMin, mBCKMax;

    //bool mRotating;
    //Vector3 mAxisOrigin, mAngularVelocity;
};

class Blocks : public std::vector<Block*>
{
public:
    Blocks() {}

protected:

private:
};

#endif // INCLUDED_HIRO_BLOCK_H__

