Prerequisite
============

gus.mb requires the following programs for compilation and execution:

* CMake 2.8 or later, preferrably 3.0.
* SWIG 2.0 or later.
* CGNSLib 3.2.1 or later.

How To Configure & Build CGNSLib For gus.mb
===========================================

1. Download cgnslib_3.2.1.tar.gz and extract the tarball.

2. In the same directory as you found this very file,

        mkdir Libs

3. Go to the cgnslib source directory extracted in 1, and run

        ccmake .

4. Enable the following switches

        CGNS_BUILD_SHARED
        CGNS_ENABLE_64BIT
        CGNS_ENABLE_SCOPING <--- important
        CGNS_USE_SHARED

   and set

        CGNS_INSTALL_PREFIX /path/to/gus.mb/Libs

5. Configure (type c) and generate (type g) Makefile.

6. Build and install

        make
        make install

How To Build gus.mb
===================

0. We assume executables, shared libraries, and Python/shell scripts
   will all be installed in "path/to/gus.mb/build" directory.

1. Configure paths to CGNSLib, SWIG, and MPI environment.

        ccmake .

2. Set the following CMake variables.

        CGNS_INCLUDE_DIR     path/to/gus.mb/Libs/include
        CGNS_LIBRARY         path/to/gus.mb/Libs/lib/libcgns.so
        CMAKE_CXX_COMPILER   path/to/your/MPI/C++/compiler  (Note 1)
        CMAKE_C_COMPILER     path/to/your/MPI/C/compiler    (Note 1)
        CMAKE_INSTALL_PREFIX path/to/gus.mb/build
        SWIG_DIR             /usr/share/swig/2.0.12         (Note 2)
        SWIG_EXECUTABLE      /usr/bin/swig                  (Note 2)

   Note 1. They should be automatically set if MPI environment
           is properly configured.
   Note 2. They should be automatically set too.

3. Build and install.

        make
        make install

4. You might need to manually copy Libs/lib/libcgns.so to the build directory.

Parallel Execution
==================

mpirun -n NUMBER_OF_SLOTS /path/to/gus.mb/build/Launch YourFlowCase.in

