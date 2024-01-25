/**
 * @file TPGPatternGenerator.hpp 
 * Pattern Generator classes for TPG applications
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef TPGTOOLS_INCLUDE_TPGTOOLS_TPG_PATTERN_ALGORITHMS_HPP_
#define TPGTOOLS_INCLUDE_TPGTOOLS_TPG_PATTERN_ALGORITHMS_HPP_

#include "tpgtools/TPGPatternGenerator.hpp"

/*
#include "fdreadoutlibs/wibeth/WIBEthFrameProcessor.hpp"
#include "hdf5libs/HDF5FileLayout.hpp"
#include "fddetdataformats/WIBEthFrame.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "logging/Logging.hpp"

// system
#include "CLI/CLI.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>


#include <cstring>
#include <immintrin.h>
#include <cstdio>
#include <array>
#include <chrono>
#include <stdio.h>
#include <stdint.h>
#include <memory>

#include <map>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
*/

// =================================================================
//                     PATTERN GENERATION ALGORITHMS
// =================================================================
struct PattgenAlgs 
{
  PattgenAlgs(){};
void
pattgen_function_golden(PattgenInfo& info) {
  std::cout << "********** GENERATED PATTERN: PATT_GOLDEN " << std::endl;
  int patt_time[9]{0};
  int patt_adc[9]{500, 502, 504, 505, 506, 505, 504, 502, 500};
  patt_time[0] = info.time_tick_offset;
  std::cout << "DBG pattgen info " << info.time_tick_offset << std::endl;
  int npatt = sizeof(patt_time) / sizeof(*patt_time);
  for (int ipatt=1; ipatt<npatt; ipatt++) {
    patt_time[ipatt] = info.time_tick_offset+ipatt < 64 ? info.time_tick_offset+ipatt : info.time_tick_offset-64+ipatt;
  }

  for (int itime=0; itime<64; ++itime) {
    for (int ch=0; ch<64; ++ch) {
      info.output_frame->set_adc(ch, itime, 0); // set pedestal value, make it configurable
    }
    //if (i==0 && itime < npatt && std::find(patt_time, patt_time + npatt, itime) == patt_time + npatt ) {
    if (info.iframe==0 && itime < patt_time[0]) {
      std::cout << "Nothing to do for first frame" << std::endl;
    } else {
      for (long unsigned int ipatt=0; ipatt<sizeof(patt_time)/sizeof(patt_time[0]); ipatt++) {
        if (itime == patt_time[ipatt]) info.output_frame->set_adc(info.input_ch, itime, patt_adc[ipatt]);
      }
    }
    uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
    //std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << info.input_ch << " \t\tTimeSample: " << itime <<  std::endl;
  }
}

void
pattgen_function_pulse(PattgenInfo& info) {
    std::cout << "********** GENERATED PATTERN: PULSE " << std::endl;
    for (int itime=0; itime<64; ++itime) {
        for (int ch=0; ch<64; ++ch) {
          info.output_frame->set_adc(ch, itime, 0);
        }
        if (itime == 0 && info.iframe==0) {
          std::cout << "Nothing to do for first frame" << std::endl;
        } else {
          info.output_frame->set_adc(info.input_ch, itime, 666);
        }
        uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
        //std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
    }
}

void
pattgen_function_square(PattgenInfo& info) {
    std::cout << "********** GENERATED PATTERN: EDGE_SQUARE " << std::endl;
    for (int itime=0; itime<64; ++itime) {
      for (int ch=0; ch<64; ++ch) {
        info.output_frame->set_adc(ch, itime, 0);
      }
      if (itime >= 0 && itime<=62 && info.iframe==0) {
        std::cout << "Nothing to do for first frame" << std::endl;
      } else {
        if (itime == 0) info.output_frame->set_adc(info.input_ch, itime, 500);
        if (itime == 63) info.output_frame->set_adc(info.input_ch, itime, 500);
      }
      uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
      //std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
    }
}
}


#endif // TPGTOOLS_INCLUDE_TPGTOOLS_TPG_PATTERN_ALGORITHMS_HPP_
