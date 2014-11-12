# $Id: CMakeLists.txt 291 2013-07-25 10:25:11Z kato $

include_directories($(CMAKE_CURRENT_SOURCE_DIR) ${Albastru_SOURCE_DIR}/IO)

find_package(CGNS REQUIRED)
include_directories(${CGNS_INCLUDE_DIR})

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
	BCInviscidWall.cpp
	BCViscousWall.cpp
	BCOutletStaticPressure.cpp
	BCInletTotal.cpp
	BCCopy.cpp
	BCFix.cpp
	BCEnforce2D.cpp
	BCForcePull.cpp
	Connectivity1to1.cpp
	Integrator.cpp
	FlowModel.cpp
	TurbulenceSpec.cpp
	KOmega1988.cpp
	AbuttingInterface.cpp
	PlanarMapping.cpp
	InterfaceDataExchanger.cpp
	FieldInitializer.cpp
	IterationContext.cpp
	)

#add_executable(Test Test.cpp)
#target_link_libraries(Test Core)

add_executable(TestIndices TestIndices.cpp)
target_link_libraries(TestIndices Core)

add_executable(TestStructuredExpr TestStructuredExpr.cpp)
target_link_libraries(TestStructuredExpr Core)

add_executable(ShockTube ShockTube.cpp)
target_link_libraries(ShockTube Core)

add_executable(TestFlatPlate TestFlatPlate.cpp)
target_link_libraries(TestFlatPlate IO Core ${CGNS_LIBRARIES})

add_executable(TestRamp TestRamp.cpp)
target_link_libraries(TestRamp IO Core ${CGNS_LIBRARIES})

add_executable(TestChannel TestChannel.cpp)
target_link_libraries(TestChannel IO Core ${CGNS_LIBRARIES})

add_executable(TestSlidingMesh TestSlidingMesh.cpp)
target_link_libraries(TestSlidingMesh IO Core)

add_executable(TestGMRES TestGMRES.cpp)
target_link_libraries(TestGMRES Core)

