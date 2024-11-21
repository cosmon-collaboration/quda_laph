
include_directories(${QudaLaph}/include/)

add_library(QudaLaph STATIC init_quda_laph.cc)

add_subdirectory(Dilution)
add_subdirectory(Field)
add_subdirectory(IO)
add_subdirectory(Quarks)
add_subdirectory(Quda_Interface)
add_subdirectory(tasks)
add_subdirectory(utils)

## need to set these
if(PARALLEL_BUILD)
  message(STATUS "parallel build")
  set(CMAKE_C_COMPILER mpicc)
  set(CMAKE_CXX_COMPILER mpic++)
  add_compile_definitions(ARCH_PARALLEL MPI_COMMS MULTI_GPU)
else()
  message(STATUS "serial build")
  set(CMAKE_C_COMPILER gcc)
  set(CMAKE_CXX_COMPILER g++)
  add_compile_definitions(ARCH_SERIAL)
endif()

## Need to set these here
include_directories(${BLAS_INC})
if(HOSTCBLAS STREQUAL "GSLCBLAS")
  add_compile_definitions(USE_GSL_CBLAS)
elseif(HOSTCBLAS STREQUAL "OPENBLAS")
  add_compile_definitions(USE_OPENBLAS)
else()
  message( FATAL_ERROR "Must define HOSTCBLAS to be ${VALID_BLAS_LIBS}")
endif()
