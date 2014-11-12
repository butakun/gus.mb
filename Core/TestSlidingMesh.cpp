// $Id: TestSlidingMesh.cpp 211 2012-03-29 22:26:18Z kato $

#include "Clock.h"
#include "Communicator.h"
#include "Roster.h"
#include "CylinderMesh.h"
#include "Block.h"
#include "RigidBodyMotion.h"
#include "VTKWriter.h"
#include <cstdlib>

void Rotate(Structured<double>& xyz, const RigidBodyMotion& motion)
{
    IndexRange range = xyz.GetRange();
    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                double* v = xyz(i, j, k);
                Vector3 p(v);
                Vector3 p2 = motion.GetPosition(p);
                v[0] = p2.X();
                v[1] = p2.Y();
                v[2] = p2.Z();
            }
        }
    }
}

int
main(int argc, char** argv)
{
    Clock::Initialize(0.0);
    Communicator::Initialize(&argc, &argv);

    IndexRange meshRange1(0, 0, 0, 10, 50, 50);
    IndexRange meshRange2(0, 0, 0, 10, 50, 50);
#if 0
    Block block1(1, meshRange1);
    Block block2(2, meshRange2);
    Roster::GetInstance()->RegisterBlock(0, 1);
    Roster::GetInstance()->RegisterBlock(0, 2);
#else
    Block* block1 = Block::New(1, meshRange1);
    Block* block2 = Block::New(2, meshRange2);
#endif

    CylinderMesh mesher;

    mesher.GenerateMesh(block1->XYZ(), 0.0, 1.0, 10.0, 11.0, -5.0, 10.0);
    mesher.GenerateMesh(block2->XYZ(), 1.0, 1.0, 10.0, 11.0, 0.0, 10.0);

    double omega = 10.0 * M_PI / 180.0; // 10 deg/sec
    RotationalMotion motion(Vector3(0.0, 0.0, 0.0), Vector3(-omega, 0.0, 0.0));

    VTKWriter writer1("cylinder1.vtk");
    writer1.AddMesh(block1->XYZ());
    writer1.Write();
    VTKWriter writer2("cylinder2.vtk");
    writer2.AddMesh(block2->XYZ());
    writer2.Write();

    Clock::GetInstance()->Advance(0.2);
    Rotate(block2->XYZ(), motion);
    VTKWriter writer3("cylinder2_2.vtk");
    writer3.AddMesh(block2->XYZ());
    writer3.Write();

    return EXIT_SUCCESS;
}

