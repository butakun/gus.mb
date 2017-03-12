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

find_package(CGNS REQUIRED)
include_directories(${CGNS_INCLUDE_DIR})

find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

find_package(ANN REQUIRED)

include_directories(
	${Albastru_SOURCE_DIR}/Core
	${Albastru_SOURCE_DIR}/IO
	#${Albastru_SOURCE_DIR}/Comm
	${CGNS_INCLUDE_DIR}
	)

add_executable(TestConnectivity1to1 TestConnectivity1to1.cpp)
target_link_libraries(TestConnectivity1to1 Core IO ${CGNS_LIBRARIES})

add_executable(TestConnectivity1to1Periodicity TestConnectivity1to1Periodicity.cpp)
target_link_libraries(TestConnectivity1to1Periodicity Core IO ${CGNS_LIBRARIES})

add_executable(TestAbuttingInterface TestAbuttingInterface.cpp)
target_link_libraries(TestAbuttingInterface Core IO ${CGNS_LIBRARIES})

#add_executable(SimpleSolver SimpleSolver.cpp)
#target_link_libraries(SimpleSolver IO Core ${CGNS_LIBRARIES})

add_executable(Solver Solver.cpp)
target_link_libraries(Solver IO Core ${ANN_LIBRARIES} ${CGNS_LIBRARIES} ${CGAL_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES})

install(TARGETS Solver RUNTIME DESTINATION .)
#install(FILES ${CGNS_LIBRARIES} DESTINATION .)

