cmake_minimum_required(VERSION 3.18)

project(QudaLaph
    VERSION 1.0
    DESCRIPTION "Lattice QCD applications using QUDA and the stochastic Laph method"
    LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(ARCH "PARALLEL" CACHE STRING "Must call with either -DARCH=PARALLEL or -DARCH=SERIAL")
set(HOSTCBLAS "OPENBLAS" CACHE STRING "Must define HOSTCBLAS")
set(BLAS_INC "/usr/include/openblas/" CACHE STRING "Includ directory for BLAS library")
set(BLAS_LIB "/usr/lib/openblas/" CACHE STRING "Lib directory for blas library")
set(GPUTOOLKIT_DIR "/usr/local/cuda-11.4/" CACHE STRING "Specify where GPYU toolkit is")
set(QUDA_DIR "" CACHE STRING "Tell us where the Quda install is")
option(OPENMP "Compile in openmp?" ON )

set(CMAKE_BUILD_TYPE "Release" )
set(CMAKE_VERBOSE_MAKEFILE ON)

## default compile options
add_compile_options(-march=native -Wall)

if(OPENMP)
  add_compile_options(-fopenmp)
endif()

if(ARCH STREQUAL "PARALLEL")
  message(STATUS "parallel build")
  set(CMAKE_C_COMPILER mpicc)
  set(CMAKE_CXX_COMPILER mpic++)
  add_compile_definitions(ARCH_PARALLEL)
  add_compile_definitions(MPI_COMMS)
  add_compile_definitions(MULTI_GPU)
elseif(ARCH STREQUAL "SERIAL")
  message(STATUS "serial build")
  set(CMAKE_C_COMPILER gcc)
  set(CMAKE_CXX_COMPILER g++)
  add_compile_definitions(ARCH_SERIAL)
else()
   message( FATAL_ERROR "Must define ARCH to be PARALLEL or SERIAL" )
endif()

if(NOT GPUTOOLKIT_DIR)
   message(FATAL_ERROR "GPUTOOLKIT_DIR must be set")
endif()

if(NOT QUDA_DIR)
   message(FATAL_ERROR "QUDA_DIR must be set")
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

add_executable(quda_laph)

target_sources(quda_laph PRIVATE
    quda_laph.cc
    byte_handler.cc
    dilution_handler.cc
    dilution_scheme_info.cc
    field_ops.cc
    field_smearing_info.cc
    filelist_info.cc
    gauge_configuration_handler.cc
    gauge_configuration_info.cc
    gluon_smearing_handler.cc
    inverter_info.cc
    io_fm_handler.cc
    laph_eigen_info.cc
    laph_noise.cc
    laph_noise_info.cc
    laph_stdio.cc
    latt_field.cc
    layout_info.cc
    perambulator_handler.cc
    quark_action_info.cc
    quark_handler.cc
    quark_smearing_handler.cc
    quda_info.cc
    read_gauge_field.cc
    stop_watch.cc
    task_quark_line_ends.cc
    task_quark_perambulators.cc
    task_smear_gauge_field.cc
    task_smear_quark_field.cc
    utils.cc
    xml_handler.cc
    )

target_link_libraries(quda_laph PUBLIC cblas)
target_link_libraries(quda_laph PUBLIC gputoolkitlib)
target_link_libraries(quda_laph PUBLIC quda)

## BLAS studd can go down here
include_directories(${BLAS_INC})
link_directories(${BLAS_LIB})
if(HOSTCBLAS STREQUAL "GSL")
   message(STATUS "using GSL cblas on hosts")
   add_compile_definitions(USE_GSL_CBLAS)
   target_link_libraries(quda_laph PUBLIC -lgsl -lgslcblas)
elseif(HOSTCBLAS STREQUAL "OPENBLAS")
   message(STATUS "using OPENblas on hosts")
   add_compile_definitions(USE_OPENBLAS)
   target_link_libraries(quda_laph PUBLIC -lopenblas)
else()
   message( FATAL_ERROR "Must define HOSTCBLAS to be GSL or OPENBLAS" )
endif()
