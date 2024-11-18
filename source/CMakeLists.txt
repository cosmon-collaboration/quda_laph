
include_directories(${QudaLaph}/include/)

add_library(QudaLaph STATIC
  byte_handler.cc
  dilution_handler.cc
  dilution_scheme_info.cc
  field_ops.cc
  field_smearing_info.cc
  filelist_info.cc
  gauge_configuration_handler.cc
  gauge_configuration_info.cc
  gluon_smearing_handler.cc
  init_quda_laph.cc
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

## need to set these
if(PARALLEL_BUILD)
  message(STATUS "parallel build")
  set(CMAKE_C_COMPILER mpicc)
  set(CMAKE_CXX_COMPILER mpic++)
  add_compile_definitions(ARCH_PARALLEL)
  add_compile_definitions(MPI_COMMS)
  add_compile_definitions(MULTI_GPU)
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
