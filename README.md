# quda_laph
Lattice QCD computations with stochastic Laph and distillation using the QUDA library.

Date: January 26, 2025
Version 0.9  (not yet 1.0)

Quda_laph carries out lattice QCD computations involving Laplacian-Heaviside (LapH) smeared
Wilson-type quarks.  It utilizes the QUDA library in order to carry out its calculations on GPUs.
It uses no other libraries, except HDF5 (in the future) and those needed by QUDA.
An installed library of QUDA is needed, compiled using the CUDA or ROCm toolkit.
The installed QUDA library does NOT need to be compiled with any of the USQCD helper libraries, 
such as QIO, QMP, QDP, or Chroma.  QUDA must be compiled with QUDA_INTERFACE_QDP=ON since QDP
ordering of the lattice sites is used on the hosts, but all other USQCD interfaces can be 
safely turned OFF. Either a serial or an MPI-based parallel implementation with or without 
OpenMP threads can be compiled.

IMPORTANT:  currently, you must use the develop branch of QUDA

Current remarks:
  -- This latest version uses a HostGlobal singleton instead of a NamedObjMap for
     persistent memory.  The user no longer needs to use  NOM_id for saving into
     memory.  
  -- Currently, there is a BUG in the split-grid communicator of QUDA that causes
     issues with peer-to-peer communications.  The bug can easily be avoided by using
     an MPI rank distribution that does not keep separate time slices on the same node.

There are few BLAS routines used on the CPU, so either gsl or openblas should be available
for the quda_laph compile.  There are ruby build scripts which provide cmake compilation
information and set needed macros.

Compilation macros:

- define either ARCH_PARALLEL or ARCH_SERIAL for parallel (MPI) or serial architecture
- define either USE_GSL_CBLAS or USE_OPENBLAS to use BLAS on the CPU

Once a quda_laph executable is created, it can be run as follows:

- (parallel)  mpirun -n 4 --host host1,host2,host3,host4  quda_laph  -npartitions 1 1 2 2 -i input.xml >& output.log
- (serial)    quda_laph -i input.xml >& output.log

The tasks to be performed by quda_laph must be specified in an input XML file.  The header
files provide information on the necessary form of the XML file.  Sample input XML files
are provided.

Remarks:

- Currently, only CLS gauge configurations can be read.  Support for SZIN format will come soon.

- A HostGlobal singleton is available for persistent memory between different tasks.

- An IOMap is used for input/output.  When specifying a path+filename, three different behaviors can occur:

      <FileName>/path/filename</FileName>             --  read/write to file using fstreams(serial) MPI-IO (parallel)
      <FileName>/path/filename[grouppath]</FileName>  --  read/write to an HDF5 file/group (Not yet supported)

- Tasks currently available are 
      - smearing the gauge field
      - computing the LapH eigenvectors 
      - computing quark perambulators
      - computing stochastic LapH quark sinks
