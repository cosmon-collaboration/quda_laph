cmake_minimum_required(VERSION 3.18)

project(QudaLaph
  VERSION 1.0
  DESCRIPTION "Lattice QCD applications using QUDA and the stochastic Laph method"
  LANGUAGES CXX C)

set(LiveLaphLibs "")

set(CMAKE_CXX_STANDARD 17)
include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

option(PARALLEL_BUILD "Compile with mpicc/mpic++" ON )

option(LAPH_DOMAIN_WALL "Compile the domain wall version" OFF )
if(LAPH_DOMAIN_WALL)
  add_compile_definitions(LAPH_DOMAIN_WALL)
endif()

set(VALID_BLAS_LIBS
  OPENBLAS
  GSLCBLAS)
set(HOSTCBLAS "OPENBLAS" CACHE STRING "Choose the host BLAS library ${VALID_BLAS_LIBS}")
set_property(CACHE HOSTCBLAS PROPERTY STRINGS ${VALID_BLAS_LIBS})

set(BLAS_INC "/usr/include/openblas/" CACHE STRING "Include directory for BLAS library")
set(BLAS_LIB "/usr/lib/openblas/" CACHE STRING "Lib directory for blas library")
set(GPUTOOLKIT_DIR "/usr/local/cuda-12.4/" CACHE STRING "Specify where GPU toolkit is")
set(QUDA_DIR "" CACHE STRING "Tell us where the Quda install is")

option(OPENMP "Compile/link openmp" ON )

set(VALID_BUILD_TYPES
  DEVEL
  RELEASE
  DEBUG)
set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Choose the type of build, options are: ${VALID_BUILD_TYPES}")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${VALID_BUILD_TYPES})

set(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(CXX)

set(CMAKE_INSTALL_PREFIX "." )

## default compile options
add_compile_options(-march=native -Wall)

if(OPENMP)
  add_compile_options(-fopenmp)
  add_compile_definitions(OPENMP)
endif()

if(PARALLEL_BUILD)
  message(STATUS "parallel build")
  set(CMAKE_C_COMPILER mpicc)
  set(CMAKE_CXX_COMPILER mpic++)
  add_compile_definitions(MPI_COMMS ARCH_PARALLEL MULTI_GPU)
else()
  message(STATUS "serial build")
  set(CMAKE_C_COMPILER gcc)
  set(CMAKE_CXX_COMPILER g++)
  add_compile_definitions(ARCH_SERIAL)
endif()

set(QUDA_LIBDIR ${QUDA_DIR}/lib)
set(QUDA_INCDIR ${QUDA_DIR}/include)

#  Include the GPUTOOLKIT library
set(GPUTOOLKIT_INCDIR ${GPUTOOLKIT_DIR}/include)
set(GPUTOOLKIT_LIBDIR ${GPUTOOLKIT_DIR}/lib64)
add_library(gputoolkitlib INTERFACE IMPORTED )
set_target_properties( gputoolkitlib PROPERTIES IMPORTED_LOCATION ${GPUTOOLKIT_LIBDIR} )

include_directories(${QUDA_INCDIR})
include_directories(${GPUTOOLKIT_INCDIR})
link_directories(${QUDA_LIBDIR})
link_directories(${GPUTOOLKIT_LIBDIR})

include_directories(include)
add_subdirectory(source)

add_executable(quda_laph quda_laph.cc)

## link openmp gomp library
if(OPENMP)
  string(APPEND LiveLaphLibs "-lgomp")
endif()

## BLAS studd can go down here
include_directories(${BLAS_INC})
link_directories(${BLAS_LIB})
if(HOSTCBLAS STREQUAL "GSLCBLAS")
   message(STATUS "using GSL cblas on hosts")
   add_compile_definitions(USE_GSL_CBLAS)
   string(APPEND LiveLaphLibs " -lcblas -lgsl -lgslcblas")
elseif(HOSTCBLAS STREQUAL "OPENBLAS")
   message(STATUS "using OPENblas on hosts")
   add_compile_definitions(USE_OPENBLAS)
   string(APPEND LiveLaphLibs " -lopenblas")
else()
   message( FATAL_ERROR "Must define HOSTCBLAS to be ${VALID_BLAS_LIBS}")
endif()

## compile and link
target_link_libraries(quda_laph PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

## option for building the tests
option(BUILD_TESTS "Build the test binaries" ON )
if(BUILD_TESTS)
  add_subdirectory(tests)
endif()
