# $Id: CMakeLists.txt 307 2013-10-02 07:03:39Z kato $

find_package(CGNS REQUIRED)
include_directories(${CGNS_INCLUDE_DIR})

include_directories(
	${Albastru_SOURCE_DIR}/Core
	${Albastru_SOURCE_DIR}/IO
	#${Albastru_SOURCE_DIR}/Comm
	${CGNS_INCLUDE_DIR}
	)

add_executable(TestConnectivity1to1 TestConnectivity1to1.cpp)
target_link_libraries(TestConnectivity1to1 Core)

add_executable(TestAbuttingInterface TestAbuttingInterface.cpp)
target_link_libraries(TestAbuttingInterface Core IO)

#add_executable(SimpleSolver SimpleSolver.cpp)
#target_link_libraries(SimpleSolver IO Core ${CGNS_LIBRARIES})

add_executable(Solver Solver.cpp)
target_link_libraries(Solver IO Core ${CGNS_LIBRARIES})

install(TARGETS Solver RUNTIME DESTINATION .)
#install(FILES ${CGNS_LIBRARIES} DESTINATION .)

