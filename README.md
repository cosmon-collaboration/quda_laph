# quda_laph
Lattice QCD computations with stochastic Laph and distillation using the QUDA library.

Quda_laph carries out lattice QCD computations involving Laplacian-Heaviside (LapH) smeared
Wilson-type quarks.  It utilizes the QUDA library in order to carry out its calculations on GPUs.
It uses no other libraries, except HDF5 (in the future) and those needed by QUDA.
An installed library of QUDA is needed, compiled using the CUDA or ROCm toolkit.
The installed QUDA library does NOT need to be compiled with any of the USQCD libraries, 
such as QIO, QMP, QDP, or Chroma.  QUDA must be compiled with QUDA_INTERFACE_QDP=ON since QDP
ordering of the lattice sites is used on the hosts, but all other USQCD interfaces can be 
safely turned OFF. Either a serial or an MPI-based parallel implementation with or without 
OpenMP threading can be compiled.

IMPORTANT:  currently, you must use the feature/trlm_3d branch!!

There are few BLAS routines used on the CPU, so either openblas or gslcblas should be available
for the quda_laph compile set via the ccmake.

## Building

Here we recommend use of ccmake, create a build directory, cd into it and run ccmake <path to source directory>
this will open a dialog where you can configure options such as:

     PARALLEL_BUILD, OPENMP, HOSTCBLAS, BLAS_LIB, BLAS_INC, GPUTOOLKIT_DIR, QUDA_DIR, CMAKE_BUILD_TYPE

There is usually some simple information so that you can choose which options suit you. This will then perform the
includes, set the macros in Quda_Laph and do the appropriate linking.

## Running

Once a quda_laph executable is created, it can be run as follows:

- (parallel)  mpirun -n 4 --host host1,host2,host3,host4  quda_laph  -npartitions 1 1 2 2 -i input.xml >& output.log
- (serial)    quda_laph -i input.xml >& output.log

The tasks to be performed by quda_laph must be specified in an input XML file.  The header
files provide information on the necessary form of the XML file.  Sample input XML files
are provided.

## Remarks:

- Currently, only "CERN" gauge configurations can be read.  Support for SZIN format will come soon. As an intermediate step you can use https://github.com/RJHudspith/GLU to convert NERSC, Scidac, ILDG, HiRep to CERN format

- A NamedObjMap is available for persistent memory between different tasks.

- An IOMap is used for input/output.  When specifying a path+filename, three different behaviors can occur:

      <FileName>NOM_gauge_field</FileName>            --  read/write from NamedObjMap memory
      <FileName>/path/filename</FileName>             --  read/write to file using fstreams(serial) MPI-IO (parallel)
      <FileName>/path/filename[grouppath]</FileName>  --  read/write to an HDF5 file/group (Not yet supported)

- Tasks currently available are 
      - smearing the gauge field, 
      - computing the LapH eigenvectors, 
      - computing quark perambulators,
      - computing stochastic LapH quark sinks,
