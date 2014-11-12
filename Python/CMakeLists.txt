# $Id: CMakeLists.txt 191 2012-01-20 02:38:43Z kato $

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
find_package(PythonInterp)
find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})
set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(MPI.i PROPERTIES CPLUSPLUS ON)
swig_add_module(MPI python MPI.i MPIFuncs.cpp)

install(FILES _MPI.so MPI.py DESTINATION .)
#install(FILES P Launch Tools DESTINATION . PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

