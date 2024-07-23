cmake_minimum_required(VERSION 3.6)

project(QudaLaph
    VERSION 1.0
    DESCRIPTION "Lattice QCD applications using QUDA and the stochastic Laph method"
    LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

option(ARCH "Must call with either -DARCH=PARALLEL or -DARCH=SERIAL" "PARALLEL")

if(ARCH STREQUAL "PARALLEL")
   message(STATUS "parallel build")
   add_compile_definitions(ARCH_PARALLEL)
   add_compile_definitions(MPI_COMMS)
   add_compile_definitions(MULTI_GPU)
elseif(ARCH STREQUAL "SERIAL")
   message(STATUS "serial build")
   add_compile_definitions(ARCH_SERIAL)
else()
   message( SEND_ERROR "Must define ARCH to be PARALLEL or SERIAL" )
endif()

if(HOSTCBLAS STREQUAL "GSL")
   message(STATUS "using GSL cblas on hosts")
   add_compile_definitions(USE_GSL_CBLAS)
elseif(HOSTCBLAS STREQUAL "OPENBLAS")
   message(STATUS "using OPENblas on hosts")
   add_compile_definitions(USE_OPENBLAS)
else()
   message( SEND_ERROR "Must define HOSTCBLAS to be GSL or OPENBLAS" )
endif()

if(NOT CUDA_DIR)
   message(FATAL_ERROR "CUDA_DIR must be set")
endif()
if(NOT QUDA_DIR)
   message(FATAL_ERROR "QUDA_DIR must be set")
endif()

set(QUDA_LIBDIR ${QUDA_DIR}/lib)
set(QUDA_INCDIR ${QUDA_DIR}/include)

#  Include the CUDA library
set(CUDA_INCDIR ${CUDA_DIR}/include)
set(CUDA_LIBDIR ${CUDA_DIR}/lib64)
add_library(cudalib INTERFACE IMPORTED )
set_target_properties( cudalib PROPERTIES IMPORTED_LOCATION ${CUDA_LIBDIR} )

include_directories(${QUDA_INCDIR})
include_directories(${CUDA_INCDIR})
link_directories(${QUDA_LIBDIR})
link_directories(${CUDA_LIBDIR})

#  openmp library
#set(OMP_LIB /usr/lib64/libgomp.so.1)
#add_library(omplib SHARED IMPORTED)
#set_target_properties( omplib PROPERTIES IMPORTED_LOCATION ${OMP_LIB} )

#  GSL library
include_directories(${GSL_INCDIR})

set(QUDA_DIR ${QUDA_DIR}/lib/cmake/QUDA)
find_package(QUDA REQUIRED)
#find_package(MPI REQUIRED)

if(BUILD_TESTS)
   add_compile_definitions(TESTING)
   add_compile_definitions(SAFETY_FLAG)
endif()

add_executable(quda_laph)

target_sources(quda_laph PRIVATE
    quda_laph.cc
    byte_handler.cc
    dilution_handler.cc
    dilution_scheme_info.cc
    field_smearing_info.cc
    filelist_info.cc
    gauge_configuration_handler.cc
    gauge_configuration_info.cc
    gluon_smearing_handler.cc
    inverter_info.cc
    io_handler_fm.cc
    laph_eigen_info.cc
    laph_noise.cc
    laph_noise_info.cc
    laph_stdio.cc
    latt_field.cc
    layout_info.cc
    quark_action_info.cc
    quark_handler.cc
    quark_smearing_handler.cc
    quda_info.cc
    read_gauge_field.cc
    stop_watch.cc
    task_quark_line_ends.cc
    task_smear_gauge_field.cc
    task_smear_quark_field.cc
    utils.cc
    xml_handler.cc
    )

if(ARCH STREQUAL "SERIAL")

add_executable(laph_data_cmp)

target_sources(laph_data_cmp PRIVATE
    laph_data_cmp.cc
    byte_handler.cc
    field_smearing_info.cc
    io_handler_fm.cc
    laph_stdio.cc
    latt_field.cc
    layout_info.cc
    xml_handler.cc
    )

endif()

if(BUILD_TESTS)
   target_sources(quda_laph PRIVATE
       task_tests.cc
       test_array.cc
       test_spin_conversions.cc
       test_filelistinfo.cc
       test_iohandler.cc
       test_latt_field.cc
       test_layout.cc
       test_nom.cc
       test_read_gauge.cc
       test_xmlhandler.cc
       )
endif()

target_link_libraries(quda_laph PUBLIC gslcblas)
target_link_libraries(quda_laph PUBLIC cudalib)
target_link_libraries(quda_laph PUBLIC "quda")
