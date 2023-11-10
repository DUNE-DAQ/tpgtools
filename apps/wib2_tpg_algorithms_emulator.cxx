/**
 * @file TestTPGAlgorithms.cxx Main file for testing different tpg algorithms 
 * 
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

// DUNE-DAQ
#include "fddetdataformats/WIB2Frame.hpp"
#include "fddetdataformats/WIBEthFrame.hpp"

#include "iomanager/IOManager.hpp"

#include "readoutlibs/utils/RateLimiter.hpp"
#include "readoutlibs/utils/FileSourceBuffer.hpp"
#include "readoutlibs/utils/BufferedFileWriter.hpp"
#include "fdreadoutlibs/DUNEWIBSuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/wib2/WIB2FrameProcessor.hpp"

#include "fdreadoutlibs/wib2/tpg/ProcessAVX2.hpp"
#include "fdreadoutlibs/wib2/tpg/ProcessRSAVX2.hpp"
#include "fdreadoutlibs/wib2/tpg/ProcessNaive.hpp"
#include "fdreadoutlibs/wib2/tpg/ProcessNaiveRS.hpp"


#include "readoutlibs/models/DefaultRequestHandlerModel.hpp"

#include "triggeralgs/TriggerPrimitive.hpp"


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



// =================================================================
//                       PUBLIC VARIABLES
// =================================================================

struct swtpg_output{
  uint16_t* output_location;
  uint64_t timestamp;
};


unsigned int total_hits = 0;
bool first_hit = true;

dunedaq::fdreadoutlibs::WIB2FrameHandler fh(0);
int duration_test = 120; // default value


std::function<void(swtpg_wib2::ProcessingInfo<swtpg_wib2::NUM_REGISTERS_PER_FRAME>& info, size_t channel_offset)> m_assigned_tpg_algorithm_function;

std::string select_algorithm = "";
std::string select_implementation = "";
bool save_adc_data = false;
bool save_trigprim = false;

// =================================================================
//                       FUNCTIONS and UTILITIES
// =================================================================

// Set CPU affinity of the processing thread
void SetAffinityThread(int executorId) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(executorId, &cpuset);
    int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    if (rc != 0) {
       std::cerr << "Error calling pthread_setaffinity_np Readout: " << rc << "\n";
    }
}

// Function save the TP data to a file 
void save_hit_data( triggeralgs::TriggerPrimitive trigprim, std::string source_name ){
  std::ofstream out_file; 

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);  
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M");
  auto date_time_str = oss.str();

  std::string file_name = "TP_dump_" + source_name + "_" + date_time_str + ".txt";
  out_file.open(file_name.c_str(), std::ofstream::app);

  //offline channel, start time, time over threshold [ns], peak_time, ADC sum, amplitude    
  out_file << trigprim.channel << "," << trigprim.time_start << "," << trigprim.time_over_threshold << "," 
	   << trigprim.time_peak << "," << trigprim.adc_integral << ","  << trigprim.adc_peak << "\n";  

  out_file.close();
}


// Function to save raw ADC data to a file (only for debugging) 
void save_raw_data(swtpg_wib2::MessageRegisters register_array, 
	       uint64_t t0, int channel_number,
           std::string source_name)
{
  std::ofstream out_file;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);  
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M");
  auto date_time_str = oss.str();

  
  std::string file_name;
  if (channel_number == -1) {
    file_name = "all_channels_" + source_name + "_data" + date_time_str + ".txt";
  } else {
    file_name = "Channel_" + std::to_string(channel_number) + "_" + source_name + "_data" + date_time_str + ".txt";
  }
  out_file.open(file_name.c_str(), std::ofstream::app);

  uint64_t t_current= t0 ; 
  
  const uint16_t* input16 = register_array.data();
  for (size_t ichan = 0; ichan < swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER; ++ichan) {
    const int register_index = ichan / swtpg_wib2::SAMPLES_PER_REGISTER;
    if (register_index < 0 || register_index >= swtpg_wib2::NUM_REGISTERS_PER_FRAME)
       continue;

    // Parse only selected channel number. To select all channels choose -1
    if (ichan == channel_number || channel_number == -1) { 
   
      const size_t register_offset = ichan % swtpg_wib2::SAMPLES_PER_REGISTER;
      const size_t register_t0_start = register_index * swtpg_wib2::SAMPLES_PER_REGISTER * swtpg_wib2::FRAMES_PER_MSG;
  
      for (size_t iframe = 0; iframe<swtpg_wib2::FRAMES_PER_MSG; ++iframe) {
    
        const size_t msg_index = iframe / 12;
        const size_t msg_time_offset = iframe % 12;
        // The index in uint16_t of the start of the message we want // NOLINT 
        const size_t msg_start_index = msg_index * swtpg_wib2::ADCS_SIZE / sizeof(uint16_t); // NOLINT
        const size_t offset_within_msg = register_t0_start + swtpg_wib2::SAMPLES_PER_REGISTER * msg_time_offset + register_offset;
        const size_t index = msg_start_index + offset_within_msg;
    
        int16_t adc_value = input16[index];
        //std::cout << adc_value << std::endl;
        out_file << ichan << "," <<  adc_value << "," << t_current << std::endl;
        t_current += 32;
      } 

    }
  }
  out_file.close();


}


// =================================================================
//                       TPG FUNCTIONS
// =================================================================

void extract_hits_naive(uint16_t* output_location, uint64_t timestamp) {

    constexpr int clocksPerTPCTick = 32;
    //uint16_t chan[100], hit_end[100], hit_charge[100], hit_tover[100]; 
    uint16_t chan, hit_end, hit_charge, hit_tover; 
    unsigned int nhits = 0;

    size_t i = 0;
    while (*output_location != swtpg_wib2::MAGIC) {
      chan   = *output_location++;
      hit_end    = *output_location++;
      hit_charge  = *output_location++;
      hit_tover     = *output_location++;


      //if (hit_charge && chan != swtpg_wib2::MAGIC) {
      //  std::cout << "Channel number: " << chan << std::endl;
      //  std::cout << "Hit charge: " << hit_charge << std::endl;
      //}

      
      i += 1;
      uint64_t tp_t_begin =                                                        
        timestamp + clocksPerTPCTick * (int64_t(hit_end ) - hit_tover );       
      uint64_t tp_t_end = timestamp + clocksPerTPCTick * int64_t(hit_end );      

      triggeralgs::TriggerPrimitive trigprim;
      trigprim.time_start = tp_t_begin;
      trigprim.time_peak = (tp_t_begin + tp_t_end) / 2;

      trigprim.time_over_threshold = hit_tover  * clocksPerTPCTick;


      trigprim.channel = chan;
      trigprim.adc_integral = hit_charge ;
      trigprim.adc_peak = hit_charge  / 20;
      trigprim.detid = 666; 
      trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
      trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
      trigprim.version = 1;
    
      if (save_trigprim) {
        save_hit_data(trigprim, "NAIVE");
      }
      ++total_hits;
      
      

    }


}

void extract_hits_avx(uint16_t* output_location, uint64_t timestamp) {

  constexpr int clocksPerTPCTick = 32;
  uint16_t chan[16], hit_end[16], hit_charge[16], hit_tover[16]; 

  while (*output_location != swtpg_wib2::MAGIC) {
    for (int i = 0; i < 16; ++i) {
      chan[i] = *output_location++; 
    }
    for (int i = 0; i < 16; ++i) {
      hit_end[i] = *output_location++; 
    }
    for (int i = 0; i < 16; ++i) {
      hit_charge[i] = *output_location++;
    }
    for (int i = 0; i < 16; ++i) {        
      hit_tover[i] = *output_location++; 
    }  
    
    // Now that we have all the register values in local
    // variables, loop over the register index (ie, channel) and
    // find the channels which actually had a hit, as indicated by
    // nonzero value of hit_charge
    for (int i = 0; i < 16; ++i) {
      if (hit_charge[i] && chan[i] != swtpg_wib2::MAGIC) {

          uint64_t tp_t_begin =                                                        // NOLINT(build/unsigned)
            timestamp + clocksPerTPCTick * (int64_t(hit_end[i]) - hit_tover[i]);       // NOLINT(build/unsigned)
          uint64_t tp_t_end = timestamp + clocksPerTPCTick * int64_t(hit_end[i]);      // NOLINT(build/unsigned)

          // May be needed for TPSet:
          // uint64_t tspan = clocksPerTPCTick * hit_tover[i]; // is/will be this needed?
          //

          // For quick n' dirty debugging: print out time/channel of hits.
          // Can then make a text file suitable for numpy plotting with, eg:
          //
          // sed -n -e 's/.*Hit: \(.*\) \(.*\).*/\1 \2/p' log.txt  > hits.txt
          //
          //TLOG_DEBUG(0) << "Hit: " << tp_t_begin << " " << offline_channel;

          triggeralgs::TriggerPrimitive trigprim;
          trigprim.time_start = tp_t_begin;
          trigprim.time_peak = (tp_t_begin + tp_t_end) / 2;
          trigprim.time_over_threshold = hit_tover[i] * clocksPerTPCTick;
          trigprim.channel = chan[i];
          trigprim.adc_integral = hit_charge[i];
          trigprim.adc_peak = hit_charge[i] / 20;
          trigprim.detid = 666;          
          trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
          trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
          trigprim.version = 1;

          if (save_trigprim){
            save_hit_data(trigprim, "AVX");
          }          



        ++total_hits;
      }
    } // loop over 16 registers   
  } // while not magic    

}



void execute_tpg(const dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter* fp) {

  SetAffinityThread(0);

  // Parse the DUNEWIB frames
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIB2Frame*>((uint8_t*)fp);
  uint64_t timestamp = wfptr->get_timestamp();      
  swtpg_wib2::MessageRegisters registers_array;
  swtpg_wib2::expand_wib2_adcs(fp, &registers_array, 0);

  
  if (first_hit) {                     
    fh.m_tpg_processing_info->setState(registers_array);
    first_hit = false;    

    // Save ADC info
    if (save_adc_data){
      save_raw_data(registers_array, timestamp, -1, select_algorithm + "_" + select_implementation);
    }


  }  
       
  
  fh.m_tpg_processing_info->input = &registers_array;
  uint16_t* destination_ptr = fh.get_hits_dest();
  *destination_ptr = swtpg_wib2::MAGIC;
  fh.m_tpg_processing_info->output = destination_ptr;
  //swtpg_wib2::process_window_avx2(*fh.m_tpg_processing_info, 0);       
  m_assigned_tpg_algorithm_function(*fh.m_tpg_processing_info, 0);
  //m_assigned_tpg_algorithm_function(*fh.m_tpg_processing_info, swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER);
          


  // Insert output of the AVX processing into the swtpg_output 
  swtpg_output swtpg_processing_result = {destination_ptr, timestamp};
    
  if (select_implementation == "AVX") {
    extract_hits_avx(destination_ptr, timestamp);
  } else if(select_implementation == "NAIVE") {  
    extract_hits_naive(destination_ptr, timestamp);
  }
   


}


// =================================================================
//                       MAIN
// =================================================================


int 
main(int argc, char** argv)
{

    CLI::App app{ "Test TPG algorithms" };

    // Set default input frame file
    std::string frame_file_path = "./wib2-frames.bin";
    app.add_option("-f,--frame_file_path", frame_file_path, "Path to the input frame file");

    app.add_option("-a,--algorithm", select_algorithm, "TPG Algorithm (SimpleThreshold / AbsRS)");
  
    app.add_option("-i,--implementation", select_implementation, "TPG implementation (AVX / NAIVE)");

    app.add_option("-d,--duration_test", duration_test, "Duration (in seconds) to run the test");    

    int num_frames = -1;
    app.add_option("-n,--num_frames", num_frames, "Number of frames to read. Default: select all frames.");

    int swtpg_threshold = 100;
    app.add_option("-t,--swtpg_threshold", swtpg_threshold, "Value of the TPG threshold");

    app.add_flag("--save_adc_data", save_adc_data, "Save ADC data");

    app.add_flag("--save_trigprim", save_trigprim, "Save trigger primitive data");


    CLI11_PARSE(app, argc, argv);

    if (select_algorithm == "SimpleThreshold") {
      if (select_implementation == "NAIVE") {
        m_assigned_tpg_algorithm_function = &swtpg_wib2::process_window_naive<swtpg_wib2::NUM_REGISTERS_PER_FRAME>;
        std::cout << "Created an instance of the " << select_algorithm << " algorithm ( " << select_implementation << " )" << std::endl;
      } else if (select_implementation == "AVX") {
        m_assigned_tpg_algorithm_function = &swtpg_wib2::process_window_avx2<swtpg_wib2::NUM_REGISTERS_PER_FRAME>;
        std::cout << "Created an instance of the " << select_algorithm << " algorithm ( " << select_implementation << " )" << std::endl;
      } else {
        std::cout << "Select a valid algorithm implementation. Use --help for further details." << std::endl;
        return 1;
      }
    } else if (select_algorithm == "AbsRS") {
      if (select_implementation == "NAIVE") {        
        //m_assigned_tpg_algorithm_function = &swtpg_wib2::process_window_naive_RS<swtpg_wib2::NUM_REGISTERS_PER_FRAME>;
        std::cout << "Created an instance of the " << select_algorithm << " algorithm ( " << select_implementation << " )" << std::endl;
      } else if (select_implementation == "AVX") {
        m_assigned_tpg_algorithm_function = &swtpg_wib2::process_window_rs_avx2<swtpg_wib2::NUM_REGISTERS_PER_FRAME>;
        std::cout << "Created an instance of the " << select_algorithm << " algorithm ( " << select_implementation << " )" << std::endl;
      } else {
        std::cout << "Select a valid algorithm implementation. Use --help for further details." << std::endl;
        return 1;
      }
    } else {
      std::cout << "Select at least an algorithm. Use --help for further details." << std::endl;
      return 1;
    }


    
    // =================================================================
    //                       READ THE FRAMES.BIN FILE
    // =================================================================
    const std::string input_file = frame_file_path;
    std::unique_ptr<dunedaq::readoutlibs::FileSourceBuffer> m_source_buffer;
    m_source_buffer = std::make_unique<dunedaq::readoutlibs::FileSourceBuffer>(10485100, swtpg_wib2::SUPERCHUNK_FRAME_SIZE);
 
    m_source_buffer->read(input_file);
    auto& source = m_source_buffer->get(); 
    int num_superchunks = m_source_buffer->num_elements(); // file_ size/chunk_size = 180 

    std::cout << "Number of superchunks in the input file: " << num_superchunks << std::endl;

    // Check if the selected number of frames is <= than the ones available in the input file
    if (num_superchunks < num_frames) {
      std::cout << "\n**ERROR**: Select a valid number of frames that is less or equal to the ones available in the input file." << std::endl;
      return 1;
    } else if (num_frames == -1) {
      num_frames = num_superchunks;
    }


    // =================================================================
    //                       Setup the SWTPG
    // =================================================================

    auto limiter = dunedaq::readoutlibs::RateLimiter(166);
    limiter.init();

    fh.initialize(swtpg_threshold);

    // =================================================================
    //                       Process the DUNEWIB superchunks
    // =================================================================
    int superchunk_index = 0; 
    uint64_t frame_repeat = 0;
    auto start_test = std::chrono::high_resolution_clock::now();  

    // Loop over the DUNEWIB superchunks in the file
    while (superchunk_index < num_frames ){      
    //while (true){        
      
      // current superchunk
      auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter*>(source.data() + superchunk_index*swtpg_wib2::SUPERCHUNK_FRAME_SIZE);

      execute_tpg(fp);

      ++superchunk_index;
      //if (superchunk_index % 10 == 0) {
      //  std::cout << "Executing superchunk number " << superchunk_index << " out of " << num_superchunks << std::endl;        
      //}

      // If end of the file is reached, measure the elapsed time
      if (superchunk_index == num_frames) {
        superchunk_index = 0;
	frame_repeat++;
      }
      if (frame_repeat % 500 == 0) {
        // Calculate elapsed time in seconds  
        auto now = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(now - start_test).count();  
        std::cout << "Elapsed time [s]: " << elapsed_seconds << std::endl;
	frame_repeat = 0;

        // stop the testing after a time a condition
        if (elapsed_seconds > duration_test) {
          superchunk_index = num_frames;
        }
	
      }
	

      limiter.limit();

    }

    std::cout << "\n\n===============================" << std::endl;
    std::cout << "Found in total " << total_hits << " hits." << std::endl;
    
    std::cout << "\n\nFinished testing." << std::endl;



}


