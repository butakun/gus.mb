#undef NDEBUG
#include "Communicator.h"
#include "BitmapMeshInterpolator.h"
#include "Block.h"
#include "BlockPatch.h"
#include "Structured.h"
#include "VTKWriter.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>

void CreateCylinderMesh(Structured<double>& XYZ, const IndexRange& meshRange, double z1, double z2, double r1, double r2, double theta1, double theta2)
{
    std::cout << XYZ.GetRange() << " == " << meshRange << std::endl;
    assert(XYZ.GetRange() == meshRange);

    // I -> axial direction, J -> radial direction, K -> circumferential direction
    // Thetas are measure from positive Y axis
    IndexIJK shape = meshRange.Shape();
    for (IndexIterator itor(meshRange); !itor.IsEnd(); itor.Advance())
    {
        IndexIJK ijk = itor.Index();
        double xi, eta, zeta;
        xi   = double(ijk.I - meshRange.Start.I) / double(meshRange.End.I - meshRange.Start.I);
        eta  = double(ijk.J - meshRange.Start.J) / double(meshRange.End.J - meshRange.Start.J);
        zeta = double(ijk.K - meshRange.Start.K) / double(meshRange.End.K - meshRange.Start.K);
        double x, y, z, r, theta;
        x = (1.0 - xi) * z1 + xi * z2;
        r = (1.0 - eta) * r1 + eta * r2;
        theta = (1.0 - zeta) * theta1 + zeta * theta2;
        y = r * std::cos(theta);
        z = r * std::sin(theta);
        XYZ(ijk)[0] = x;
        XYZ(ijk)[1] = y;
        XYZ(ijk)[2] = z;
    }
}

void WriteVTKPoints(const char* filename, const std::vector<Vector3>& pixels)
{
    std::ofstream f(filename);

    f << "# vtk DataFile Version 2.0" << std::endl
      << "pixels" << std::endl
      << "ASCII" << std::endl
      << "DATASET POLYDATA" << std::endl
      << "POINTS " << pixels.size() << " float" << std::endl;
    for (std::vector<Vector3>::const_iterator i = pixels.begin();
        i != pixels.end(); ++i)
    {
        const Vector3& p = *i;
        f << p.X() << ' ' << p.Y() << ' ' << p.Z() << std::endl;
    }
}

int main(int argc, char** argv)
{
    Communicator::Initialize(&argc, &argv);

#if 1
    const int B1_IEND = 5, B1_JEND = 11, B1_KEND = 21;
    const int B2_IEND = 5, B2_JEND = 17, B2_KEND = 33;
    const double DTHETADEG = 60.0;
#else
    const int B1_IEND = 1, B1_JEND = 1, B1_KEND = 2;
    const int B2_IEND = 1, B2_JEND = 1, B2_KEND = 2;
    const double DTHETADEG = 10.0;
#endif

    IndexRange mr1(IndexIJK(0, 0, 0), IndexIJK(B1_IEND, B1_JEND, B1_KEND));
    IndexRange mr2(IndexIJK(0, 0, 0), IndexIJK(B2_IEND, B2_JEND, B2_KEND));

    Block* block1 = Block::New(1, mr1);
    Block* block2 = Block::New(2, mr2);

    CreateCylinderMesh(
        block1->XYZ(), mr1,
        0.0, 1.0, 1.0, 1.1, 0.0, DTHETADEG * M_PI / 180.0
        );

    CreateCylinderMesh(
        block2->XYZ(), mr2,
        1.0, 2.0, 1.0, 1.1, 0.0, DTHETADEG * M_PI / 180.0
        );

    VTKWriter* writer;
    writer = new VTKWriter("cylinder1.vtk");
    writer->AddMesh(block1->XYZ());
    writer->Write();
    delete writer;
    writer = new VTKWriter("cylinder2.vtk");
    writer->AddMesh(block2->XYZ());
    writer->Write();
    delete writer;

    std::cout << "Block 1: XYZ.Data = " << block1->XYZ().Data << std::endl;
    std::cout << "Block 2: XYZ.Data = " << block2->XYZ().Data << std::endl;

    BlockPatch bp1 = BlockPatch::New(1, IndexRange(B1_IEND, 0, 0, B1_IEND, B1_JEND, B1_KEND), 1);
    BlockPatch bp2 = BlockPatch::New(2, IndexRange(0, 0, 0, 0, B2_JEND, B2_KEND), 2);

    BitmapMeshInterpolator ipr; 
    ipr.AddPatch(0, bp1, block1->XYZ());
    ipr.AddPatch(1, bp2, block2->XYZ());

    ipr.GenerateBitmapPixels();

    WriteVTKPoints("pixels.vtk", ipr.GetPixels());

    const std::vector<BitmapMeshInterpolator::PatchIndex>& pixelMapping = ipr.GetPixelMapping();
    for (size_t i = 0; i < pixelMapping.size(); ++i)
    {
        const BitmapMeshInterpolator::PatchIndex& pi = pixelMapping[i];
        std::cout << "Pixel " << i << ": " << pi.PatchID[0] << ' ' << pi.PatchIJK[0] << " => " << pi.PatchID[1] << ' ' << pi.PatchIJK[1] << std::endl;
    }

    ipr.ComputeCellFaceWeights();
    std::map<int, Structured<double> >& cfws = ipr.GetCellFaceWeights();
    for (int ibp = 0; ibp < 2; ++ibp)
    {
        const BlockPatch& bp = ibp == 0 ? bp1 : bp2;
        const Structured<double>& cfw = cfws[bp.UniqueID()];
        const IndexRange& cfr = cfw.GetRange();
        std::cout << "BlockPatch: " << bp << std::endl;
        for (IndexIterator itor(cfr); !itor.IsEnd(); itor.Advance())
        {
            const IndexIJK& ijk = itor.Index();
            std::cout << " " << ijk << " = " << cfw(ijk)[0] << std::endl;
        }
    }
}

