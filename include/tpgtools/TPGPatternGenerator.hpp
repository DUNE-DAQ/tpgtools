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
struct PattgenInfo
{
  PattgenInfo(dunedaq::fddetdataformats::WIBEthFrame* output_frame_
		, int num_frames_
		, std::string pattern_name_
                , uint64_t first_timestamp_
                , int time_tick_offset_
		, std::string out_prefix_
                , int input_ch_
		, int tpg_threshold_
                , int iframe_
                )
  : output_frame(output_frame_)
  , num_frames(num_frames_)
  , pattern_name(pattern_name_)
  , first_timestamp(first_timestamp_)
  , time_tick_offset(time_tick_offset_)
  , out_prefix(out_prefix_)
  , input_ch(input_ch_)
  , tpg_threshold(tpg_threshold_)
  , iframe(iframe_)
  {}
  dunedaq::fddetdataformats::WIBEthFrame* output_frame;
  int num_frames;
  std::string pattern_name;
  uint64_t first_timestamp;
  int time_tick_offset;
  std::string out_prefix;
  int input_ch;
  int tpg_threshold;
  int iframe;
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
                    , 0
		    , 0
                    , 0
                    );
  }
};

// =================================================================
//                     PATTERN GENERATION ALGORITHMS
// =================================================================
void
pattgen_function_golden(PattgenInfo& info) 
{
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
      info.output_frame->set_adc(ch, itime, 0); // NB set pedestal value, make it configurable
    }
    if (info.iframe==0 && itime < patt_time[0]) {
      std::cout << "Nothing to do for first frame" << std::endl;
    } else {
      for (long unsigned int ipatt=0; ipatt<sizeof(patt_time)/sizeof(patt_time[0]); ipatt++) {
        if (itime == patt_time[ipatt]) info.output_frame->set_adc(info.input_ch, itime, patt_adc[ipatt]);
      }
    }
    uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
    //std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << itime << " \t\tChannel: " << info.input_ch << " \t\tTimeSample: " << itime <<  std::endl;
  }
}

void
pattgen_function_pulse(PattgenInfo& info) 
{
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
        //std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << itime << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
    }
}

void
pattgen_function_square(PattgenInfo& info) 
{
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
      //std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << itime << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
    }
}

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

  void init() {
    tpgpg_map["patt_golden"] = &pattgen_function_golden;
    tpgpg_map["patt_pulse"] = &pattgen_function_pulse;
    tpgpg_map["patt_square"] = &pattgen_function_square;
  } 

  void run_algorithm(std::string& alg_name, PattgenInfo& info) {
    tpgpg_map[alg_name](info);
  }
};

// =================================================================
//                       PATTERN GENERATOR
// =================================================================
void 
execute_tpgpg(PattgenInfo& info, PattgenAlgs& pa) 
{
  // Parse the WIBEth frames
  uint64_t timestamp = info.output_frame->get_timestamp();
  uint64_t expected_timestamp = info.first_timestamp + info.iframe*2048;
  if (info.iframe>0 && timestamp != expected_timestamp) {
    int64_t ts_diff = timestamp - expected_timestamp;
    std::cout << " |______ TIMESTAMP WILL BE OVERWRITTEN! OLD: " << timestamp << ", NEW: " << expected_timestamp <<  std::endl;
    info.output_frame->set_timestamp(expected_timestamp);
  }

  pa.run_algorithm(info.pattern_name, info);
}

// =================================================================
//                       PATTERN GENERATOR VALIDATION
// =================================================================
void hit_finder(std::vector<uint16_t>& adcs, std::vector<std::vector<int>>& out, const int& channel, 
		const int& threshold, const uint64_t timestamp) {

  //int threshold = (int)m_threshold;

      // index of ADCs above threshold
      //int threshold = (int)m_threshold;
      std::cout << "DBG threshold " << threshold << std::endl;
      std::cout << "DBG ADC sample length " << adcs.size() << std::endl;
      for(long unsigned int m = 0; m < adcs.size(); m++){
        std::cout << "DBG ADC sample, value " << m << ", " << adcs[m] << std::endl;
      }

      int m_tov_min = 2;

      bool dbg = false;
      int sz = adcs.size();
      if (sz == 64) {
        for (auto& t : adcs) std::cout << "DBG pattgen 64 adcs " << t << std::endl;
      }

      std::vector<int> igt;
      std::vector<uint16_t>::iterator it_adcs = adcs.begin();
      while ((it_adcs = find_if(it_adcs, adcs.end(), [&threshold](int x){return x > threshold; })) != adcs.end())
      {
        igt.push_back(distance(adcs.begin(), it_adcs));
        // cout << distance(adcs.begin(), it_adcs) << endl;
        it_adcs++;
      }

      if (dbg && sz == 64) std::cout << "DBG pattgen 64 igt " << igt.size() << std::endl;

      // cout << "igt.size = " << igt.size() << endl;
      // cout << (igt.size() < m_tov_min) << endl;

      if(igt.size() < m_tov_min){
        // cout << "m_tov_min condition not fulfilled for whole packet!" << endl;
        return void();
      }

      std::vector<int> igt_diff;
      adjacent_difference (igt.begin(), igt.end(), back_inserter(igt_diff));
      igt_diff.erase(igt_diff.begin());

      // find start and end of hits
      std::vector<int> istart;
      std::vector<int> iend;
      istart.push_back(0);
      std::vector<int>::iterator it_igt = igt_diff.begin();
      while ((it_igt = find_if(it_igt, igt_diff.end(), [ ](int x){return x != 1; })) != igt_diff.end())
      {
        istart.push_back(distance(igt_diff.begin(), it_igt)+1);
        iend.push_back(distance(igt_diff.begin(), it_igt));
        it_igt++;
      }
      iend.push_back(igt.size()-1);

      std::vector<int> start;
      std::vector<int> end;
      std::vector<int> hitcontinue;
      for(long unsigned int i = 0; i < istart.size(); i++){
        start.push_back(igt[istart[i]]);
        end.push_back(igt[iend[i]]);
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

      if (dbg && sz == 64) {
	std::cout << "DBG pattgen 64 hit sizes " << channel << " : " << start.size() << " " << end.size() << " " << peak_times.size() << " " << sums.size() << " " << peak_adcs.size() << " " << hitcontinue.size() << std::endl;
      }

      // check output hits fullfil the m_tov_min condition
      for(long unsigned int k = 0; k < start.size(); k++){
        if (peak_adcs[k] > 16384) {
          if (dbg && sz == 64) std::cout << "DBG pattgen 64 peak adc is too big " << peak_adcs[k] << std::endl;
          continue; // 2**14
        }
        if (end[k]-start[k] >= m_tov_min-1){
          //std::vector<int> aux = {start[k], end[k], peak_times[k], peak_adcs[k], sums[k], hitcontinue[k]};
          //std::vector<int> aux = {k, start[k], end[k], peak_times[k], peak_adcs[k], sums[k], hitcontinue[k]};
          // start, end, peak_time, channel, sum_adc, peak_adc, hitcotninue
	  std::vector<int> aux = {start[k], end[k], peak_times[k], channel, sums[k], peak_adcs[k], hitcontinue[k]};
          out.push_back(aux);

	  std::cout << "# hit " << "[" << adcs.size() << "]"<< k << ", start time " << start[k] << ", end time " << end[k] << ", peak time " << peak_times[k] << ", peak adc " << peak_adcs[k] << ", sum adc " << sums[k] << ", hit continue " << hitcontinue[k] << ", timestamp " << timestamp << std::endl;
          
	  uint64_t ts_tov = (end[k]-start[k])*32;
	  uint64_t ts_start = start[k] * 32 + timestamp;
	  uint64_t ts_peak = peak_times[k] * 32 + timestamp;
	  std::cout << channel << "," << ts_start << "," << ts_tov << "," << ts_peak << "," << sums[k] << "," << peak_adcs[k] << std::endl;
	}
      }
}

	
void
execute_tpgpg_val(PattgenInfo& info)
{
  std::string out_suffix = "";
  std::cout << "Run pedestal subtraction " << std::endl;

  bool do_pedsub = true;  // pedestal subtraction is mandatory but can be switched off only for debugging
  bool do_filter = false; // currently no FIR filter applied 

  // pedsub as stream - and reused headers for writing binary file 
  dunedaq::fddetdataformats::WIBEthFrame* output_frame_pedsub; 
  FrameFile input_file_fake = FrameFile((info.out_prefix+"_wibeth_output.bin").c_str()); 
  std::cout << "Size of the input file wibeth_output.bin: " << input_file_fake.length() << std::endl;
  std::cout << "Number of frames " << info.num_frames << std::endl;

  std::fstream output_file_pedsub;
  if (do_filter) {
    out_suffix = "_fir";
  }
  output_file_pedsub.open(info.out_prefix+"_wibeth_output_pedsub"+out_suffix+".bin", std::ios::app | std::ios::binary);
    // pedestal subtraction
  int16_t median_ssr[64] = { 0 };
  int16_t accum_ssr[64] = { 0 };
  int16_t median = 0;
  int16_t accum = 0;
  // filtering
  /*
  const size_t NTAPS = 8;
  uint16_t absTimeModNTAPS;
  int16_t* prev_samp = new int16_t[64*64];
  for (size_t it = 0; it < 64 * 64; ++it) {
    for (size_t jt = 0; jt < NTAPS; ++jt) {
      prev_samp[it * NTAPS + jt] = 0;
    }
  }
  const uint8_t tap_exponent = 6;
  int16_t multiplier{1 << tap_exponent};
  int16_t adcMax{INT16_MAX / multiplier};
  std::vector<int16_t> tpg_taps = swtpg_wibeth::firwin_int(7, 0.1, multiplier);
  int16_t* taps = new int16_t[tpg_taps.size()];
  if (do_filter) {
    for (size_t it = 0; it < tpg_taps.size(); ++it) {
      std::cout << "DBG fir filter taps " << it << " : " << tpg_taps[it] << std::endl;
      taps[it] = tpg_taps[it];
    }
    std::cout << "DBG fir filter consts: multiplier, adcMax, INT16_MAX: " << multiplier << ", " << adcMax << ", " << INT16_MAX << std::endl;
  }
  */

  for (size_t i=0; i<info.num_frames; i++) {
    std::cout << "========== FRAME_NUM " << i <<  std::endl;
    output_frame_pedsub = input_file_fake.frame(i);
    // pedsub
    for (int ch=0; ch<64; ++ch) {
      if (i==0) {
        median_ssr[ch] = output_frame_pedsub->get_adc(ch, 0);
      }	
      median = median_ssr[ch];
      accum = accum_ssr[ch];
      for (int itime=0; itime<64; ++itime) {
        int16_t sample = output_frame_pedsub->get_adc(ch, itime);
	//std::cout << "DBG pedsub write to file sample 1 " << sample << ", " << median << ", " << accum << std::endl;
	if (ch == info.input_ch) {
	  //std::cout << "DBG frugal before " << median << ", " << accum << ", " << sample << std::endl;
          swtpg_wibeth::frugal_accum_update(median, sample, accum, 10);
	  //std::cout << "DBG frugal after " << median << ", " << accum << ", " << sample << std::endl;
          median_ssr[ch] = median;
          accum_ssr[ch] = accum;
	  sample -= median;
	  //std::cout << "DBG frugal after sample " << sample << std::endl;
	  // filtering
	  /*if (do_filter) {
            int16_t filt = swtpg_wibeth::fir_filter(sample, adcMax, NTAPS, absTimeModNTAPS, taps, prev_samp); 
	    std::cout << "Filtering before / after: " << sample << " / " << filt << " where " << absTimeModNTAPS << std::endl;
	    sample = filt >> tap_exponent;
          }*/
	}
	//std::cout << "DBG pedsub write to file sample 2: " << sample << ", " << median << ", " << accum << ", " << adcMax << std::endl;
	// NB TDAQ ERROR ADC value out of range 
	if (sample > 0 && sample < INT16_MAX) {
	  output_frame_pedsub->set_adc(ch, itime, sample);
	} else {
          output_frame_pedsub->set_adc(ch, itime, 0);
        }
        //std::cout << "Pedsub ADC value: " << sample << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime << " \t\tiCh: " << ch <<  std::endl;
      }
    }
    output_file_pedsub.write(reinterpret_cast<char*>(output_frame_pedsub), sizeof(dunedaq::fddetdataformats::WIBEthFrame) );  
    }
  output_file_pedsub.close();

  // pedestal subtraction - ADCs stored in long vector - works for small number of frames -> text file 
  // todo
  //
  std::cout << "Run hit finding " << std::endl;

  // hit finding - ADCs stored in long vector - works for small number of frames -> text file
  dunedaq::fddetdataformats::WIBEthFrame* input_frame_pedsub; 
  std::string input_file_name = info.out_prefix+"_wibeth_output.bin";

  std::string file_name_hits = info.out_prefix+"_wibeth_output_hits.txt";
  if (do_pedsub) {
    input_file_name = info.out_prefix+"_wibeth_output_pedsub"+out_suffix+".bin";
    file_name_hits = info.out_prefix+"_wibeth_output_pedsub"+out_suffix+"_hits.txt";
  }
  FrameFile input_file_pedsub = FrameFile(input_file_name.c_str()); 
  std::cout << "Size of the input file " << input_file_name << ": " << input_file_pedsub.length() << std::endl;
  std::cout << "Number of frames " << info.num_frames << std::endl;

  std::ofstream output_file_pedsub_hits;
  output_file_pedsub_hits.open(file_name_hits.c_str(), std::ofstream::app);

  for (int ch=0; ch<64; ++ch) {
    if (ch != info.input_ch) continue;
    //std::vector<int16_t> adcs;
    //std::vector<int> adcs;
    std::vector<uint16_t> adcs;
    std::vector<std::vector<int>> tmp_out;
    uint64_t timestamp;
    for (size_t i=0; i<info.num_frames; i++) {
      std::cout << "========== FRAME_NUM " << i << " max " << std::vector<int>().max_size() <<  std::endl;
      //std::cout << "========== FRAME_NUM " << i << " max " << std::vector<std::vector<int> >().max_size() <<  std::endl;
      input_frame_pedsub = input_file_pedsub.frame(i);
      if (i==0) {
	timestamp = input_frame_pedsub->get_timestamp();
      }
      for (int itime=0; itime<64; ++itime) {
	//int16_t adc = input_frame_pedsub->get_adc(ch, itime);
	//int adc = (int)input_frame_pedsub->get_adc(ch, itime);
	uint16_t adc = input_frame_pedsub->get_adc(ch, itime);
	adcs.push_back(adc);
      }
    }
    hit_finder(adcs, tmp_out, ch, info.tpg_threshold, timestamp);
    for (auto& it : tmp_out) {
      // start[k], end[k], peak_times[k], channel, sums[k], peak_adcs[k], hitcontinue[k]
      uint64_t ts_tov = (it[1]-it[0])*32;
      uint64_t ts_start = it[0] * 32 + timestamp;
      uint64_t ts_peak = it[2] * 32 + timestamp;
      output_file_pedsub_hits << it[3] << "," << ts_start << "," << ts_tov << "," << ts_peak << "," << it[4] << "," << it[5] << "\n";
      //++total_hits;
    }
  }
  output_file_pedsub_hits.close();
} // end save_trigprim



#endif // TPGTOOLS_INCLUDE_TPGTOOLS_TPG_PATTERN_GENERATOR_HPP_
