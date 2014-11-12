// $Id: TestAbuttingInterface.cpp 257 2012-12-05 07:48:51Z kato $

#include "Communicator.h"
#include "Block.h"
#include "AbuttingInterface.h"
#include "PlanarMapping.h"
#include "Structured.h"
#include "VTKWriter.h"
#include <cmath>

void BaumKuchen(Structured<double>& XYZ, double rHub, double rShroud, double xUp, double xDown, double theta0, double dTheta)
{
    IndexRange range = XYZ.GetRange();
    for (IndexIterator i(range); !i.IsEnd(); i.Advance())
    {
        IndexIJK ijk = i.Index();
        double xi = double(ijk.I - range.Start.I) / double(range.End.I - range.Start.I);
        double eta = double(ijk.J - range.Start.J) / double(range.End.J - range.Start.J);
        double zeta = double(ijk.K - range.Start.K) / double(range.End.K - range.Start.K);
        double theta = (1.0 - zeta) * theta0 + zeta * (theta0 + dTheta);
        double r = (1.0 - eta) * rHub + eta * rShroud;
        double x = (1.0 - xi) * xUp + xi * xDown;
        double* xyz = XYZ(ijk);
        xyz[0] = x;
        xyz[1] = r * std::cos(theta);
        xyz[2] = r * std::sin(theta);
    }
}

void Rotate(Structured<double>& xyz, const Vector3& rot)
{
    double R[3][3];
    Vector3::RotationMatrix(R, rot);
    for (IndexIterator i(xyz.GetRange()); !i.IsEnd(); i.Advance())
    {
        double* p = xyz(i.Index());
        Vector3 p1(p);
        p1.Apply(R);
        p[0] = p1.X();
        p[1] = p1.Y();
        p[2] = p1.Z();
    }
}

int main(int argc, char** argv)
{
    Communicator::Initialize(&argc, &argv);

    IndexRange mr1(0, 0, 0, 4, 19, 19);
    IndexRange mr2(0, 0, 0, 4, 19, 19);

    Block* block1 = Block::New(1, mr1);
    Block* block2 = Block::New(2, mr2);

    Structured<double>& XYZ1 = block1->XYZ();
    Structured<double>& XYZ2 = block2->XYZ();
    double rHub = 1.0, rShroud = 2.0;
    double x1 = 0.0, x2 = 1.0, x3 = 2.0;
    double theta1 = 0.0, dTheta1 = 10.0 * M_PI / 180.0;
    double theta2 = 5.0 * M_PI / 180.0, dTheta2 = 10.0 * M_PI / 180.0;
    //double theta2 = 0.0 * M_PI / 180.0, dTheta2 = 10.0 * M_PI / 180.0;

    block1->SetPeriodicity(Vector3(dTheta1, 0.0, 0.0));
    block2->SetPeriodicity(Vector3(dTheta2, 0.0, 0.0));

    BaumKuchen(XYZ1, rHub, rShroud, x1, x2, theta1, dTheta1);
    BaumKuchen(XYZ2, rHub, rShroud, x2, x3, theta2, dTheta2);
    //Rotate(XYZ2, Vector3(dTheta2, 0.0, 0.0));

    VTKWriter writer1("baumkuchen_1.vtk");
    writer1.AddMesh(XYZ1);
    writer1.Write();
    VTKWriter writer2("baumkuchen_2.vtk");
    writer2.AddMesh(XYZ2);
    writer2.Write();

    BlockPatch patch1 = BlockPatch::New(1, IndexRange(mr1.End.I, mr1.Start.J, mr1.Start.K, mr1.End.I, mr1.End.J, mr1.End.K), 1);
    BlockPatch patch2 = BlockPatch::New(2, IndexRange(mr2.Start.I, mr2.Start.J, mr2.Start.K, mr2.Start.I, mr2.End.J, mr2.End.K), 2);
//    BlockPatch patch1(1, IndexRange(mr1.End.I, mr1.Start.J, mr1.Start.K, mr1.End.I, mr1.End.J, mr1.End.K));
//    BlockPatch patch2(2, IndexRange(mr2.Start.I, mr2.Start.J, mr2.Start.K, mr2.Start.I, mr2.End.J, mr2.End.K));
    AbuttingInterface::BlockPatches patches1, patches2;
    patches1.push_back(patch1);
    patches2.push_back(patch2);
    AbuttingInterface* interface = AbuttingInterface::New(patches1, patches2);
    interface->SetPatchMesh(1, XYZ1);
    interface->SetPatchMesh(2, XYZ2);

    IterationContext iteration;
    interface->MapMesh(iteration);

    const PlanarMapping& mapper = interface->GetMapper(patch1.UniqueID());
    const Structured<int> cells = mapper.DonorCells();
    const Structured<double> dists = mapper.Distances();
    for (IndexIterator i(cells.GetRange()); !i.IsEnd(); i.Advance())
    {
        IndexIJK ijk = i.Index();
        int* cell = cells(ijk);
        int donorBlockID = cell[0];
        IndexIJK ijkd(cell[1], cell[2], cell[3]);
        double* dist = dists(ijk);
        std::cout << ijk << " = " << donorBlockID << ":" << ijkd << ", d = " << dist[0] << std::endl;
    }
}

