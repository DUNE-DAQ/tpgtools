/**
 * @file TPGPatternGenerator.hpp 
 * Pattern Generator classes for TPG applications
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef TPGTOOLS_INCLUDE_TPGTOOLS_TPG_PATTERN_GENERATOR_HPP_
#define TPGTOOLS_INCLUDE_TPGTOOLS_TPG_PATTERN_GENERATOR_HPP_

#include "fdreadoutlibs/wibeth/WIBEthFrameProcessor.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessNaive.hpp"

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

// =================================================================
//                     FRAME BINARY FILE
// =================================================================
/*
 * Helper class to handle input/output frame binary files (*.bin)
*/
class FrameFile
{
public:

    FrameFile(const char* filename)
        : m_file(filename, std::ifstream::binary),
          m_buffer(new char[sizeof(dunedaq::fddetdataformats::WIBEthFrame)])
    {
        if(m_file.bad() || m_file.fail() || !m_file.is_open()){
            throw std::runtime_error(std::string("Bad file ")+std::string(filename));
        }
        // Calculate the length of the file
        m_file.seekg(0, m_file.end);
        m_length = m_file.tellg();
        m_file.seekg(0, m_file.beg);
        if(m_length==0){
            throw std::runtime_error("Empty file");
        }
        //if(m_length%sizeof(dunedaq::fddetdataformats::WIBEthFrame)!=0){
        //    throw std::runtime_error("File does not contain an integer number of frames");
        //}
        m_n_frames=m_length/sizeof(dunedaq::fddetdataformats::WIBEthFrame);

        // Reinterpret the frame as WIBEthFrame

        m_file.read(m_buffer, m_file.eof());
        m_wibeth_frame = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(m_buffer);
    }

    ~FrameFile()
    {
        m_file.close();
        delete[] m_buffer;
    }

    // Length of the file in bytes
    size_t length() const {return m_length;}
    // Number of frames in the file
    size_t num_frames() const { return m_n_frames; }

    dunedaq::fddetdataformats::WIBEthFrame* get_wibeth_frame() const { return m_wibeth_frame; }

    dunedaq::fddetdataformats::WIBEthFrame* frame(size_t i)
    {
        if(i>=num_frames()) return nullptr;
        // Seek to the right place in the file
        m_file.seekg(i*sizeof(dunedaq::fddetdataformats::WIBEthFrame));
        // Check we didn't go past the end
        if(m_file.bad() || m_file.eof()) return nullptr;
        // Actually read the fragment into the buffer
        m_file.read(m_buffer,sizeof(dunedaq::fddetdataformats::WIBEthFrame));
        if(m_file.bad() || m_file.eof()) return nullptr;
        return reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(m_buffer);
    }

protected:
    std::ifstream m_file;
    char* m_buffer;
    dunedaq::fddetdataformats::WIBEthFrame* m_wibeth_frame = nullptr;
    size_t m_length;
    size_t m_n_frames;
};

// =================================================================
//                     PATTERN GENERATION INFO
// =================================================================
/*
 * Contain all relevant parameters and pass that info to algorithms
*/
struct PattgenInfo
{
  PattgenInfo(dunedaq::fddetdataformats::WIBEthFrame* output_frame_
		, int num_frames_
		, std::string pattern_name_
                , uint64_t first_timestamp_
                , int clock_tick_offset_
		, std::string out_prefix_
		, std::string path_output_
                , int input_ch_
		, int tpg_threshold_
                , int iframe_
		, bool overwrite_wibeth_header_
		, bool verbose_
                )
  : output_frame(output_frame_)
  , num_frames(num_frames_)
  , pattern_name(pattern_name_)
  , first_timestamp(first_timestamp_)
  , clock_tick_offset(clock_tick_offset_)
  , out_prefix(out_prefix_)
  , path_output(path_output_)
  , input_ch(input_ch_)
  , tpg_threshold(tpg_threshold_)
  , iframe(iframe_)
  , overwrite_wibeth_header(overwrite_wibeth_header_)
  , verbose(verbose_)
  {}
  dunedaq::fddetdataformats::WIBEthFrame* output_frame;
  int num_frames;
  std::string pattern_name;
  uint64_t first_timestamp;
  int clock_tick_offset;
  std::string out_prefix;
  std::string path_output;
  int input_ch;
  int tpg_threshold;
  int iframe;
  bool overwrite_wibeth_header;
  bool verbose;
};
class PattgenHandler {
  public:
  PattgenHandler(){};
  ~PattgenHandler(){};
  std::unique_ptr<PattgenInfo> m_pattgen_info;

  void initialize() {
    m_pattgen_info = std::make_unique<PattgenInfo>(nullptr
		    , 0
		    , "patt_golden"
		    , 0
                    , 0
		    , ""
		    , ""
                    , 0
		    , 0
                    , 0
		    , false
		    , false
                    );
  }
};

// =================================================================
//                     PATTERN GENERATION ALGORITHMS
// =================================================================
/*
 * Extendable class containing the existing algorithms
 * New algorithms can be added here 
*/
struct PattgenAlgs 
{
  PattgenAlgs(std::unique_ptr<PattgenInfo>& info_)
  : info(info_)
  {
    init();
  }
  std::unique_ptr<PattgenInfo>& info;
  typedef void (*fnp)(PattgenInfo& info);
  typedef std::map<std::string,fnp> tpgpg_t;
  tpgpg_t tpgpg_map;

  static void pattgen_function_golden(PattgenInfo& info);
  static void pattgen_function_pulse(PattgenInfo& info);
  static void pattgen_function_square(PattgenInfo& info);
  static void pattgen_function_square_left(PattgenInfo& info);
  static void pattgen_function_square_right(PattgenInfo& info);

  void check_header(PattgenInfo& info);
   
  void init() {
    tpgpg_map["patt_golden"] = &pattgen_function_golden;
    tpgpg_map["patt_pulse"] = &pattgen_function_pulse;
    tpgpg_map["patt_square"] = &pattgen_function_square;
    tpgpg_map["patt_square_left"] = &pattgen_function_square_left;
    tpgpg_map["patt_square_right"] = &pattgen_function_square_right;
  } 

  void run_algorithm(PattgenInfo& info) {
    tpgpg_map[info.pattern_name](info);
  }
};

// header validation 
void
PattgenAlgs::check_header(PattgenInfo& info)
{
  uint64_t timestamp = info.output_frame->get_timestamp();
  uint64_t expected_timestamp = info.first_timestamp + info.iframe*2048;
  if (info.iframe>0 && timestamp != expected_timestamp) {
    if (info.verbose) {
      TLOG() << " |______ TIMESTAMP WILL BE OVERWRITTEN! OLD: " << timestamp << ", NEW: " << expected_timestamp;
    }
    info.output_frame->set_timestamp(expected_timestamp);
  }
  // for offline channel map, avoid error e.g. Invalid stream number 176
  if (info.overwrite_wibeth_header) {
    if (info.verbose) {
       TLOG() << "CRATE " << info.output_frame->daq_header.crate_id; 
       TLOG() << "SLOT " << info.output_frame->daq_header.slot_id; 
       TLOG() << "STREAM " << info.output_frame->daq_header.stream_id; 
    }
    //if crate number looks wrong 
    if (info.output_frame->daq_header.crate_id != 6) {
      if (info.verbose) TLOG() << " |______ CRATE-ID WILL BE OVERWRITTEN! OLD: " << info.output_frame->daq_header.crate_id << ", NEW: " << 6;
      info.output_frame->daq_header.crate_id = 6;
    }
    //if slot number looks wrong
    if (info.output_frame->daq_header.slot_id != 0) {
      if (info.verbose) TLOG() << " |______ SLOT-ID WILL BE OVERWRITTEN! OLD: " << info.output_frame->daq_header.slot_id << ", NEW: " << 0;
      info.output_frame->daq_header.slot_id = 0;
    }
    //if stream number looks wrong (not 0,1,2,3 or 64,65,66,67)
    if (info.output_frame->daq_header.stream_id != 0) {
      if (info.verbose) TLOG() << " |______ STREAM-ID WILL BE OVERWRITTEN! OLD: " << info.output_frame->daq_header.stream_id << ", NEW: " << 0;
      info.output_frame->daq_header.stream_id = 0;
    }
  }
}

// one global and two smaller peaks
void
PattgenAlgs::pattgen_function_golden(PattgenInfo& info) 
{
  TLOG() << "********** GENERATED PATTERN: PATT_GOLDEN ";

  int patt_time[9]{0};
  int patt_adc[9]{500, 502, 504, 505, 506, 505, 504, 502, 500};
  patt_time[0] = info.clock_tick_offset;
  TLOG() << "pattgen info, time tick offset: " << info.clock_tick_offset;

  int npatt = sizeof(patt_time) / sizeof(*patt_time);
  for (int ipatt=1; ipatt<npatt; ipatt++) {
    patt_time[ipatt] = info.clock_tick_offset+ipatt < 64 ? info.clock_tick_offset+ipatt : info.clock_tick_offset-64+ipatt;
  }

  for (int itime=0; itime<64; ++itime) {
    for (int ch=0; ch<64; ++ch) {
      info.output_frame->set_adc(ch, itime, 0); // NB set pedestal value, make it configurable
    }
    if (info.iframe==0 && itime < patt_time[0]) {
      TLOG() << "Nothing to do for first frame";
    } else {
      for (long unsigned int ipatt=0; ipatt<sizeof(patt_time)/sizeof(patt_time[0]); ipatt++) {
        if (itime == patt_time[ipatt]) info.output_frame->set_adc(info.input_ch, itime, patt_adc[ipatt]);
      }
    }
    if (info.verbose) {
      uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
      TLOG() << "Output ADC value: " << adc_val << "\t\t\tFrame: " << itime << " \t\tChannel: " << info.input_ch << " \t\tTimeSample: " << itime;
    }
  }
}

// single peak on single channel at single clock tick
void
PattgenAlgs::pattgen_function_pulse(PattgenInfo& info) 
{
    TLOG() << "********** GENERATED PATTERN: PULSE ";
    for (int itime=0; itime<64; ++itime) {
        for (int ch=0; ch<64; ++ch) {
          info.output_frame->set_adc(ch, itime, 0);
        }
        if (itime == 0 && info.iframe==0) {
          TLOG() << "Nothing to do for first frame";
        } else {
          info.output_frame->set_adc(info.input_ch, itime, 666);
        }
	if (info.verbose) {
          uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
          TLOG() << "Output ADC value: " << adc_val << "\t\t\tFrame: " << itime << " \t\tChannel: " << info.input_ch << " \t\tTimeSample: " << itime;
        }
    }
}

// square shape fully contained within 64 clock ticks
void
PattgenAlgs::pattgen_function_square(PattgenInfo& info) 
{
    TLOG() << "********** GENERATED PATTERN: SQUARE ";
    for (int itime=0; itime<64; ++itime) {
      for (int ch=0; ch<64; ++ch) {
        info.output_frame->set_adc(ch, itime, 0);
      }
      if (itime >= 0 && itime<=62 && info.iframe==0) {
        TLOG() << "Nothing to do for first frame";
      } else {
        if (itime == 0) info.output_frame->set_adc(info.input_ch, itime, 500);
        if (itime == 63) info.output_frame->set_adc(info.input_ch, itime, 500);
      }
      if (info.verbose) {
        uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
        TLOG() << "Output ADC value: " << adc_val << "\t\t\tFrame: " << itime << " \t\tChannel: " << info.input_ch << " \t\tTimeSample: " << itime;
      }
    }
}

// simple square shape across 2 WIBEth frames, peak is in first frame
void
PattgenAlgs::pattgen_function_square_left(PattgenInfo& info) {
    TLOG() << "********** GENERATED PATTERN: EDGE_LEFT ";
    for (int itime=0; itime<64; ++itime) {
      for (int ch=0; ch<64; ++ch) {
        info.output_frame->set_adc(ch, itime, 0);
      }
      if (itime >= 0 && itime<=62 && info.iframe==0) {
        TLOG() << "Nothing to do for first frame";
      } else {
        if (itime == 0) info.output_frame->set_adc(info.input_ch, itime, 500);
        if (itime == 63) info.output_frame->set_adc(info.input_ch, itime, 501);
      }
      if (info.verbose) {
        uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
        TLOG() << "Output ADC value: " << adc_val << "\t\t\tFrame: " << itime << " \t\tChannel: " << info.input_ch << " \t\tTimeSample: " << itime;
      }
    }
}

// simple square shape across 2 WIBEth frames, peak is in second frame
void
PattgenAlgs::pattgen_function_square_right(PattgenInfo& info) {
    TLOG() << "********** GENERATED PATTERN: EDGE_RIGHT ";
    for (int itime=0; itime<64; ++itime) {
      for (int ch=0; ch<64; ++ch) {
        info.output_frame->set_adc(ch, itime, 0);
      }
      if (itime >= 0 && itime<=62 && info.iframe==0) {
        TLOG() << "Nothing to do for first frame";
      } else {
        if (itime == 0) info.output_frame->set_adc(info.input_ch, itime, 501);
        if (itime == 63) info.output_frame->set_adc(info.input_ch, itime, 500);
      }
      if (info.verbose) {
        uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
        TLOG() << "Output ADC value: " << adc_val << "\t\t\tFrame: " << itime << " \t\tChannel: " << info.input_ch << " \t\tTimeSample: " << itime;
      }
    }
}


// =================================================================
//                       PATTERN GENERATOR
// =================================================================
/*
 * Invokes the algorithms using function pointers 
*/
void 
execute_tpgpg(PattgenInfo& info, PattgenAlgs& pa) 
{
  /*
  // Parse the WIBEth frames
  uint64_t timestamp = info.output_frame->get_timestamp();
  uint64_t expected_timestamp = info.first_timestamp + info.iframe*2048;
  if (info.iframe>0 && timestamp != expected_timestamp) {
    if (info.verbose) {
      TLOG() << " |______ TIMESTAMP WILL BE OVERWRITTEN! OLD: " << timestamp << ", NEW: " << expected_timestamp;
    }
    info.output_frame->set_timestamp(expected_timestamp);
  }
  // for offline channel map, avoid error e.g. Invalid stream number 176
  if (info.overwrite_wibeth_header) {
    if (info.verbose) {
       TLOG() << "CRATE " << info.output_frame->daq_header.crate_id; 
       TLOG() << "SLOT " << info.output_frame->daq_header.slot_id; 
       TLOG() << "STREAM " << info.output_frame->daq_header.stream_id; 
    }
    //if crate number looks wrong 
    if (info.output_frame->daq_header.crate_id != 6) {
      if (info.verbose) TLOG() << " |______ CRATE-ID WILL BE OVERWRITTEN! OLD: " << info.output_frame->daq_header.crate_id << ", NEW: " << 6;
      info.output_frame->daq_header.crate_id = 6;
    }
    //if slot number looks wrong
    if (info.output_frame->daq_header.slot_id != 0) {
      if (info.verbose) TLOG() << " |______ SLOT-ID WILL BE OVERWRITTEN! OLD: " << info.output_frame->daq_header.slot_id << ", NEW: " << 0;
      info.output_frame->daq_header.slot_id = 0;
    }
    //if stream number looks wrong (not 0,1,2,3 or 64,65,66,67)
    if (info.output_frame->daq_header.stream_id != 0) {
      if (info.verbose) TLOG() << " |______ STREAM-ID WILL BE OVERWRITTEN! OLD: " << info.output_frame->daq_header.stream_id << ", NEW: " << 0;
      info.output_frame->daq_header.stream_id = 0;
    }
  }
  */
  pa.check_header(info);
  pa.run_algorithm(info);
}

// =================================================================
//                       PATTERN GENERATOR VALIDATION
// =================================================================
/*
 * Standalone algorithm for ADC waveform analysis and signal identification, i.e. hit finding 
 * Extracts hits (or trigger primitives) in an independent way using continuous waveform as input
 * Hits are stored in a text file for reference and validation purposes 
*/
void waveform_analyser(std::vector<uint16_t>& adcs, std::vector<std::vector<int>>& out, const int& channel, 
		const int& threshold) {

  int m_tov_min = 2;

  std::vector<int> tover; // ADCs above threshold 
  std::vector<uint16_t>::iterator it_adcs = adcs.begin();
  while ((it_adcs = find_if(it_adcs, adcs.end(), [&threshold](int x){return x > threshold; })) != adcs.end())
  {
     tover.push_back(distance(adcs.begin(), it_adcs));
     it_adcs++;
  }

  if((int)tover.size() < m_tov_min){
    return void();
  }

  std::vector<int> tover_diff;
  adjacent_difference (tover.begin(), tover.end(), back_inserter(tover_diff));
  tover_diff.erase(tover_diff.begin());

  // find start and end of hits
  std::vector<int> istart;
  std::vector<int> iend;
  istart.push_back(0);
  std::vector<int>::iterator it_tover = tover_diff.begin();
  while ((it_tover = find_if(it_tover, tover_diff.end(), [ ](int x){return x != 1; })) != tover_diff.end())
  {
    istart.push_back(distance(tover_diff.begin(), it_tover)+1);
    iend.push_back(distance(tover_diff.begin(), it_tover));
    it_tover++;
  }
  iend.push_back(tover.size()-1);

  std::vector<int> start;
  std::vector<int> end;
  std::vector<int> hitcontinue;
  for(long unsigned int i = 0; i < istart.size(); i++){
    start.push_back(tover[istart[i]]);
    end.push_back(tover[iend[i]]);
    if(end[i] == 63){
      hitcontinue.push_back(1);
    }else{
      hitcontinue.push_back(0);
    }
  }

  // find hit sums
  std::vector<int> sums;
  for(long unsigned int j = 0; j < start.size(); j++){
      std::vector<uint16_t>::iterator it_start = adcs.begin()+start[j];
      std::vector<uint16_t>::iterator it_end = adcs.begin()+end[j]+1;
      int sum = accumulate(it_start, it_end, 0);
      sums.push_back(sum);
  }

  // find peak adcs and times
  std::vector<int> peak_adcs;
  std::vector<int> peak_times;
  for(long unsigned int j = 0; j < start.size(); j++){
    std::vector<uint16_t>::iterator it_start = adcs.begin()+start[j];
    std::vector<uint16_t>::iterator it_end = adcs.begin()+end[j]+1;
    std::vector<uint16_t>::iterator max = max_element(it_start, it_end);
    int peak = *max;
    int time = distance(adcs.begin(), max);
    peak_adcs.push_back(peak);
    peak_times.push_back(time);
  }

  // check output hits fullfil the m_tov_min condition
  for(long unsigned int k = 0; k < start.size(); k++){
    if (peak_adcs[k] > 16384) {
      continue; // 2**14
    }
    if (end[k]-start[k] >= m_tov_min-1){
      std::vector<int> aux = {start[k], end[k], peak_times[k], channel, sums[k], peak_adcs[k], hitcontinue[k]};
      out.push_back(aux);
    }
  }
}
	
void
execute_tpgpg_validation(PattgenInfo& info)
{
  // 1. pedestal subtraction
  dunedaq::fddetdataformats::WIBEthFrame* output_frame_pedsub; 
  FrameFile input_file_fake = FrameFile((info.path_output+"/"+info.out_prefix+"_wibeth_output.bin").c_str()); 

  std::fstream output_file_pedsub;
  output_file_pedsub.open(info.path_output+"/"+info.out_prefix+"_wibeth_output_pedsub.bin", std::ios::app | std::ios::binary);
  int16_t median_ssr[64] = { 0 };
  int16_t accum_ssr[64] = { 0 };
  int16_t median = 0;
  int16_t accum = 0;

  for (int i=0; i<info.num_frames; i++) {
    TLOG() << "========== FRAME_NUM " << i;
    output_frame_pedsub = input_file_fake.frame(i);
    for (int ch=0; ch<64; ++ch) {
      if (i==0) {
        median_ssr[ch] = output_frame_pedsub->get_adc(ch, 0);
      }	
      median = median_ssr[ch];
      accum = accum_ssr[ch];
      for (int itime=0; itime<64; ++itime) {
        int16_t sample = output_frame_pedsub->get_adc(ch, itime);
	if (ch == info.input_ch) {
          swtpg_wibeth::frugal_accum_update(median, sample, accum, 10);
          median_ssr[ch] = median;
          accum_ssr[ch] = accum;
	  sample -= median;
	}
	// IH TDAQ ERROR when ADC value out of range 
	if (sample > 0 && sample < INT16_MAX) {
	  output_frame_pedsub->set_adc(ch, itime, sample);
	} else {
          output_frame_pedsub->set_adc(ch, itime, 0);
        }
      }
    }
    output_file_pedsub.write(reinterpret_cast<char*>(output_frame_pedsub), sizeof(dunedaq::fddetdataformats::WIBEthFrame) );  
  }
  output_file_pedsub.close();

  // 2. standalone hit finding
  dunedaq::fddetdataformats::WIBEthFrame* input_frame_pedsub; 
  std::string input_file_name = info.path_output+"/"+info.out_prefix+"_wibeth_output.bin";

  std::string file_name_hits = info.path_output+"/"+info.out_prefix+"_wibeth_output_hits.txt";
  input_file_name = info.path_output+"/"+info.out_prefix+"_wibeth_output_pedsub.bin";
  file_name_hits = info.path_output+"/"+info.out_prefix+"_wibeth_output_pedsub_hits.txt";

  FrameFile input_file_pedsub = FrameFile(input_file_name.c_str()); 

  std::ofstream output_file_pedsub_hits;
  output_file_pedsub_hits.open(file_name_hits.c_str(), std::ofstream::app);

  for (int ch=0; ch<64; ++ch) {
    if (ch != info.input_ch) continue;
    std::vector<uint16_t> adcs;
    std::vector<std::vector<int>> tmp_out;
    uint64_t timestamp = 0;
    for (int i=0; i<info.num_frames; i++) {
      TLOG() << "========== FRAME_NUM " << i << " max " << std::vector<int>().max_size();
      input_frame_pedsub = input_file_pedsub.frame(i);
      if (i==0) {
	timestamp = input_frame_pedsub->get_timestamp();
      }
      for (int itime=0; itime<64; ++itime) {
	uint16_t adc = input_frame_pedsub->get_adc(ch, itime);
	adcs.push_back(adc);
      }
    }
    waveform_analyser(adcs, tmp_out, ch, info.tpg_threshold);
    for (auto& it : tmp_out) {
      // start, end, peak_times, channel, sums, peak_adcs, hitcontinue
      uint64_t ts_tov = (it[1]-it[0])*32;
      uint64_t ts_start = it[0] * 32 + timestamp;
      uint64_t ts_peak = it[2] * 32 + timestamp;
      output_file_pedsub_hits << it[3] << "," << ts_start << "," << ts_tov << "," << ts_peak << "," << it[4] << "," << it[5] << "\n";
    }
  }
  output_file_pedsub_hits.close();
}

#endif // TPGTOOLS_INCLUDE_TPGTOOLS_TPG_PATTERN_GENERATOR_HPP_
