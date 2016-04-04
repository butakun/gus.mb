#    gus.mb, an open source flow solver.
#    Copyright (C) 2016 Hiromasa Kato <hiromasa at gmail.com>
#
#    This file is part of gus.mb.
#
#    gus.mb is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    gus.mb is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
# $Id: CMakeLists.txt 291 2013-07-25 10:25:11Z kato $

include_directories($(CMAKE_CURRENT_SOURCE_DIR) ${Albastru_SOURCE_DIR}/IO)

find_package(CGNS REQUIRED)
include_directories(${CGNS_INCLUDE_DIR})

#find_package(ANN REQUIRED)
#include_directories(${ANN_INCLUDE_DIR})

add_library(Core
	Structured.cpp
	Physics.cpp
	Clock.cpp
	VirtualBlock.cpp
	Block.cpp
	Roster.cpp
	BlockPatch.cpp
	Profiler.cpp
	Communicator.cpp
	MPICommunicator.cpp
	StructuredDataExchanger.cpp
	PatchExchanger1to1.cpp
	IndexUtils.cpp
	SolverFunctions.cpp
	Residual.cpp
	ViscousFluxEvaluator.cpp
	SolverUtils.cpp
	RectMesh.cpp
	BLMesh.cpp
	RampMesh.cpp
	CylinderMesh.cpp
	BC.cpp
	BCPlanar.cpp
	BCPlanarLocal.cpp
	BCExtrapolate.cpp
	BCSymmetry.cpp
	BCViscousWall.cpp
	BCRotatingViscousWall.cpp
	BCOutletStaticPressure.cpp
	BCInletTotal.cpp
	BCCopy.cpp
	BCFix.cpp
	BCEnforce2D.cpp
	BCForcePull.cpp
	Connectivity1to1.cpp
	ConnectivityAbutting.cpp
	Integrator.cpp
	FlowModel.cpp
	TurbulenceSpec.cpp
	KOmega1988.cpp
	AbuttingInterface.cpp
	SimplePlanarAbuttingInterface.cpp
	PlanarMapping.cpp
	GeometricMapping.cpp
	InterfaceDataExchanger.cpp
	FieldInitializer.cpp
	IterationContext.cpp
	BitmapMeshInterpolator.cpp
	)

#add_executable(Test Test.cpp)
#target_link_libraries(Test Core)

add_executable(TestMatrix33 TestMatrix33.cpp)
target_link_libraries(TestMatrix33 Core IO ${CGNS_LIBRARIES})

add_executable(TestIndices TestIndices.cpp)
target_link_libraries(TestIndices Core IO ${CGNS_LIBRARIES})

add_executable(TestStructuredExpr TestStructuredExpr.cpp)
target_link_libraries(TestStructuredExpr Core IO ${CGNS_LIBRARIES})

add_executable(ShockTube ShockTube.cpp)
target_link_libraries(ShockTube Core IO ${CGNS_LIBRARIES})

add_executable(TestFlatPlate TestFlatPlate.cpp)
target_link_libraries(TestFlatPlate Core IO ${CGNS_LIBRARIES})

add_executable(TestRamp TestRamp.cpp)
target_link_libraries(TestRamp Core IO ${CGNS_LIBRARIES})

add_executable(TestChannel TestChannel.cpp)
target_link_libraries(TestChannel Core IO ${CGNS_LIBRARIES})

add_executable(TestSlidingMesh TestSlidingMesh.cpp)
target_link_libraries(TestSlidingMesh Core IO ${CGNS_LIBRARIES})

add_executable(TestGMRES TestGMRES.cpp)
target_link_libraries(TestGMRES Core IO ${CGNS_LIBRARIES})

add_executable(TestBitmapMeshInterpolator TestBitmapMeshInterpolator.cpp)
target_link_libraries(TestBitmapMeshInterpolator Core IO ${CGNS_LIBRARIES})

