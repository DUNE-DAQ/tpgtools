cmake_minimum_required(VERSION 3.12)
project(tpgtools VERSION 0.0.0)

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

##############################################################################
# Dependency sets
set(TPGTOOLS_DEPENDENCIES
#tools
  Folly::folly
  ers::ers
  logging::logging
#dunedaq
  appfwk::appfwk
  readoutlibs::readoutlibs
  fdreadoutlibs::fdreadoutlibs
  opmonlib::opmonlib
  daqdataformats::daqdataformats
  detdataformats::detdataformats
  fddetdataformats::fddetdataformats
  trigger::trigger
  triggeralgs::triggeralgs
  detchannelmaps::detchannelmaps
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

daq_add_application(wib2_tpg_algorithms_emulator wib2_tpg_algorithms_emulator.cxx LINK_LIBRARIES tpgtools hdf5libs::hdf5libs CLI11::CLI11)
#daq_add_application(wibeth_simple_tpg_emulator wibeth_simple_tpg_emulator.cxx LINK_LIBRARIES tpgtools hdf5libs::hdf5libs CLI11::CLI11)
daq_add_application(wibeth_tpg_algorithms_emulator wibeth_tpg_algorithms_emulator.cxx LINK_LIBRARIES tpgtools hdf5libs::hdf5libs CLI11::CLI11)

daq_add_application(wib2_binary_frame_modifier wib2_binary_frame_modifier.cxx LINK_LIBRARIES tpgtools hdf5libs::hdf5libs)
daq_add_application(wibeth_binary_frame_modifier wibeth_binary_frame_modifier.cxx LINK_LIBRARIES tpgtools hdf5libs::hdf5libs)

daq_add_application(wib2_binary_frame_reader wib2_binary_frame_reader.cxx LINK_LIBRARIES tpgtools hdf5libs::hdf5libs)
daq_add_application(wibeth_binary_frame_reader wibeth_binary_frame_reader.cxx LINK_LIBRARIES tpgtools hdf5libs::hdf5libs)


##############################################################################
# Installation
#daq_install()
