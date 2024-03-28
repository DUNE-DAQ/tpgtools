cmake_minimum_required(VERSION 3.12)
project(tpgtools VERSION 1.1.0)

find_package(daq-cmake REQUIRED)

daq_setup_environment()

find_package(ers REQUIRED)
find_package(logging REQUIRED)
find_package(readoutlibs REQUIRED)
find_package(fdreadoutlibs REQUIRED)
find_package(daqdataformats REQUIRED)
find_package(detchannelmaps REQUIRED)
find_package(detdataformats REQUIRED)
find_package(fddetdataformats REQUIRED)
find_package(trigger REQUIRED)
find_package(hdf5libs REQUIRED)


find_package(Boost COMPONENTS unit_test_framework REQUIRED)

##############################################################################
# Dependency sets
set(TPGTOOLS_DEPENDENCIES
#tools
  ers::ers
  logging::logging
#dunedaq
  readoutlibs::readoutlibs
  fdreadoutlibs::fdreadoutlibs
  daqdataformats::daqdataformats
  detdataformats::detdataformats
  fddetdataformats::fddetdataformats
  detchannelmaps::detchannelmaps
  hdf5libs::hdf5libs
)


##############################################################################
# Extra options and tweaks
set(TPGTOOLS_USE_INTRINSICS ON)

if(${TPGTOOLS_USE_INTRINSICS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
endif()

##############################################################################
# Main library
daq_add_library(
	  *.cpp
    LINK_LIBRARIES ${TPGTOOLS_DEPENDENCIES}
)

##############################################################################
# Applications

daq_add_application(wib2_tpg_algorithms_emulator wib2_tpg_algorithms_emulator.cxx LINK_LIBRARIES ${PROJECT_NAME} CLI11::CLI11)
daq_add_application(wibeth_tpg_validation wibeth_tpg_validation.cxx LINK_LIBRARIES tpgtools  CLI11::CLI11)
daq_add_application(wibeth_tpg_workload_emulator wibeth_tpg_workload_emulator.cxx LINK_LIBRARIES tpgtools  CLI11::CLI11)
daq_add_application(wibeth_tpg_algorithms_emulator wibeth_tpg_algorithms_emulator.cxx LINK_LIBRARIES ${PROJECT_NAME}  CLI11::CLI11)
daq_add_application(wibeth_tpg_pattern_generator wibeth_tpg_pattern_generator.cxx LINK_LIBRARIES ${PROJECT_NAME}  CLI11::CLI11)

daq_add_application(wib2_binary_frame_modifier wib2_binary_frame_modifier.cxx LINK_LIBRARIES ${PROJECT_NAME} )
daq_add_application(wibeth_binary_frame_modifier wibeth_binary_frame_modifier.cxx LINK_LIBRARIES ${PROJECT_NAME} )

daq_add_application(wib2_binary_frame_reader wib2_binary_frame_reader.cxx LINK_LIBRARIES ${PROJECT_NAME} CLI11::CLI11)
daq_add_application(wibeth_binary_frame_reader wibeth_binary_frame_reader.cxx LINK_LIBRARIES ${PROJECT_NAME} CLI11::CLI11)

daq_add_application(streamed_TPs_to_text streamed_TPs_to_text.cxx LINK_LIBRARIES  ${PROJECT_NAME} CLI11::CLI11)
daq_add_application(check_fragment_TPs check_fragment_TPs.cxx LINK_LIBRARIES  ${PROJECT_NAME} CLI11::CLI11)


##############################################################################
# Unit test

# See https://dune-daq-sw.readthedocs.io/en/latest/packages/daq-cmake/#daq_add_unit_test
daq_add_unit_test(Placeholder_test LINK_LIBRARIES)  # Placeholder_test should be replaced with real unit tests


##############################################################################
# Installation
daq_install()
