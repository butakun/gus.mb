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
# $Id: CMakeLists.txt 314 2013-12-24 15:47:30Z kato $

find_package(CGNS REQUIRED)
include_directories(${CGNS_INCLUDE_DIR})

include_directories(${Albastru_SOURCE_DIR}/Core)

add_library(IO
	PLOT3DMeshReader.cpp
	PLOT3DMeshWriter.cpp
	PLOT3DSolutionWriter.cpp
	CGNSMeshReader.cpp
	CGNSWriter.cpp
	CGNSReader.cpp
	VTKWriter.cpp
	)

#add_executable(TestCGNSMeshReader TestCGNSMeshReader.cpp)
#target_link_libraries(TestCGNSMeshReader Core IO ${CGNS_LIBRARIES})

add_executable(ConvertCGNSToPLOT3D ConvertCGNSToPLOT3D.cpp)
target_link_libraries(ConvertCGNSToPLOT3D Core IO ${CGNS_LIBRARIES})

#add_library(CGNSFuncs CGNSFuncs.cpp)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
find_package(PythonInterp)
find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})
set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(CGNS.i PROPERTIES CPLUSPLUS ON)
swig_add_module(CGNS python CGNS.i CGNSFuncs.cpp)
swig_link_libraries(CGNS ${CGNS_LIBRARIES})

install(FILES _CGNS.so CGNS.py CGNSFile.py PointRange.py LoadBalance.py GenerateBCs.py SolverSetup.py Periodicity.py PeriodicityUtils.py
	DumpCGNS.py IGGBCS.py ApplyIGGBCs.py ExportToVTK.py
	Pre.py Roster.py SceneObject.py PatchSelectionPanel.py wxVTKRenderWindowInteractor.py
	DESTINATION .)
install(FILES P Launch Tools Pre DESTINATION . PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

