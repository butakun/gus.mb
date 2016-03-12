# $Id: FindCGNS.cmake 318 2014-03-28 05:25:10Z kato $
# Find CGNSLib
# This module defines
# CGNS_INCLUDE_DIR
# CGNS_LIBRARIES
# CGNS_FOUND

find_path(CGNS_INCLUDE_DIR cgnslib.h)
find_library(CGNS_LIBRARY cgns)

if(CGNS_INCLUDE_DIR AND CGNS_LIBRARY)
  set(CGNS_FOUND TRUE)
  #find_package(HDF5)
endif(CGNS_INCLUDE_DIR AND CGNS_LIBRARY)

if(CGNS_FOUND)
  set(CGNS_LIBRARIES ${CGNS_LIBRARY} ${HDF5_LIBRARIES})
  message(STATUS "Found CGNS ${CGNS_LIBRARIES}")
  if(HDF5_FOUND)
    message(STATUS "Found HDF5 ${HDF5_LIBRARIES}")
  endif(HDF5_FOUND)
endif(CGNS_FOUND)

