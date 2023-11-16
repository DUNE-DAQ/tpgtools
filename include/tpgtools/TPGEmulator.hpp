/**
 * @file TPGEmulator.hpp 
 * Emulator classes for TPG applications
  *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef TPGTOOLS_INCLUDE_TPGTOOLS_TPGEMULATOR_HPP_
#define TPGTOOLS_INCLUDE_TPGTOOLS_TPGEMULATOR_HPP_


#include "fddetdataformats/WIBEthFrame.hpp"
#include "fddetdataformats/WIBEthFrame.hpp"
#include "iomanager/IOManager.hpp"
#include "readoutlibs/utils/RateLimiter.hpp"
#include "readoutlibs/utils/FileSourceBuffer.hpp"
#include "readoutlibs/utils/BufferedFileWriter.hpp"
#include "fdreadoutlibs/DUNEWIBEthTypeAdapter.hpp"
#include "fdreadoutlibs/wibeth/WIBEthFrameProcessor.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessAVX2.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessRSAVX2.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessNaive.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessNaiveRS.hpp"
#include "readoutlibs/models/DefaultRequestHandlerModel.hpp"
#include "detchannelmaps/TPCChannelMap.hpp"
#include "triggeralgs/TriggerPrimitive.hpp"
#include "trgdataformats/TriggerPrimitive.hpp"
#include "hdf5libs/HDF5RawDataFile.hpp"
#include "logging/Logging.hpp"

#include "TPGToolsIssues.hpp"

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;


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
void save_TP_object( triggeralgs::TriggerPrimitive trigprim, std::string algo ){
  std::ofstream out_file; 

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);  
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M");
  auto date_time_str = oss.str();

  std::string file_name = "TP_dump_" + algo + "_" + date_time_str + ".txt";
  out_file.open(file_name.c_str(), std::ofstream::app);

  //offline channel, start time, time over threshold [ns], peak_time, ADC sum, amplitude    
  out_file << trigprim.channel << "," << trigprim.time_start << "," << trigprim.time_over_threshold << "," 
	   << trigprim.time_peak << "," << trigprim.adc_integral << ","  << trigprim.adc_peak <<  ","  << trigprim.type << "\n";  

  out_file.close();
}


// Function to save raw ADC data to a file (only for debugging) 
void save_raw_data(swtpg_wibeth::MessageRegisters register_array, 
	       uint64_t t0, int channel_number,
           std::string algo)
{
  std::ofstream out_file;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);  
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M");
  auto date_time_str = oss.str();

  
  std::string file_name;
  if (channel_number == -1) {
    file_name = "all_channels_" + algo + "_data" + date_time_str + ".txt";
  } else {
    file_name = "Channel_" + std::to_string(channel_number) + "_" + algo + "_data" + date_time_str + ".txt";
  }
  out_file.open(file_name.c_str(), std::ofstream::app);

  uint64_t t_current= t0 ; 
  
  for (size_t ichan = 0; ichan < swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER; ++ichan) {

    const size_t register_index = ichan / swtpg_wibeth::SAMPLES_PER_REGISTER;
    // Parse only selected channel number. To select all channels choose -1
    if (static_cast<int>(ichan) ==channel_number || channel_number == -1) { 
   
      const size_t register_offset = ichan % swtpg_wibeth::SAMPLES_PER_REGISTER;
      const size_t register_t0_start = register_index * swtpg_wibeth::SAMPLES_PER_REGISTER * swtpg_wibeth::FRAMES_PER_MSG;
  
      for (size_t iframe = 0; iframe<swtpg_wibeth::FRAMES_PER_MSG; ++iframe) {
    
        const size_t msg_index = iframe / swtpg_wibeth::FRAMES_PER_MSG; 
        const size_t msg_time_offset = iframe % swtpg_wibeth::FRAMES_PER_MSG;
        // The index in uint16_t of the start of the message we want // NOLINT 
        const size_t msg_start_index = msg_index * (swtpg_wibeth::ADCS_SIZE) / sizeof(uint16_t); // NOLINT
        const size_t offset_within_msg = register_t0_start + swtpg_wibeth::SAMPLES_PER_REGISTER * msg_time_offset + register_offset;
        const size_t index = msg_start_index + offset_within_msg;
    
        int16_t adc_value = register_array.uint16(index);
        out_file << " Time " << iframe << " channel " <<  ichan << " ADC_value " <<  adc_value <<  " timestamp " << t_current << std::endl;

        t_current += 32; 
      } 

    }
  }
  out_file.close();


}

void
save_tp(const dunedaq::trgdataformats::TriggerPrimitive& prim, bool save_trigprim)
{
  std::ofstream out_file; 
  std::ostringstream oss;

  if (save_trigprim) {   
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);  
    
    oss << std::put_time(&tm, "%d-%m-%Y_%H-%M");
    auto date_time_str = oss.str();
  
    std::string file_name = "TriggerPrimitiveRecord_dump_" + date_time_str + ".txt";
    out_file.open(file_name.c_str(), std::ofstream::app);
  
    //offline channel, start time, time over threshold [ns], peak_time, ADC sum, amplitude    
    out_file << prim.channel << "," << prim.time_start << "," << prim.time_over_threshold << "," 
	     << prim.time_peak << "," << prim.adc_integral << ","  << prim.adc_peak << "," << prim.type << "\n";  
  
  }
  out_file.close();
  TLOG_DEBUG(TLVL_BOOKKEEPING) << "Saved TriggerPrimitives output file." ;
  
}


void process_trigger_primitive(std::unique_ptr<dunedaq::daqdataformats::Fragment>&& frag, 
               int TP_index, 
               unsigned int& total_tp_hits, 
               bool save_trigprim)
{
  size_t payload_size = frag->get_size() - sizeof(dunedaq::daqdataformats::FragmentHeader);
  size_t n_tps = payload_size / sizeof(dunedaq::trgdataformats::TriggerPrimitive);
  TLOG_DEBUG(TLVL_BOOKKEEPING) << "Trigger Primitive number " << TP_index << " with SourceID[" << frag->get_element_id() << "] has " << n_tps << " TPs" ;  
  total_tp_hits = total_tp_hits + n_tps ; 
  size_t remainder = payload_size % sizeof(dunedaq::trgdataformats::TriggerPrimitive);
  assert(remainder == 0);
  const dunedaq::trgdataformats::TriggerPrimitive* prim = reinterpret_cast<dunedaq::trgdataformats::TriggerPrimitive*>(frag->get_data());
  
  for (size_t i = 0; i < n_tps; ++i) {
    save_tp(*prim, save_trigprim);
    ++prim;
  }
}



// =================================================================
//                       TPG FUNCTIONS
// =================================================================

void extract_hits_avx(uint16_t* output_location, uint64_t timestamp,
                      std::array<uint, swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER>& register_channels,
                      unsigned int& total_hits,
                      bool save_trigprim) {

  constexpr int clocksPerTPCTick = 32;
  uint16_t chan[16], hit_end[16], hit_charge[16], hit_tover[16]; 

  while (*output_location != swtpg_wibeth::MAGIC) {
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
      if (hit_charge[i] && chan[i] != swtpg_wibeth::MAGIC) {
        //std::cout << "Channel number: " << chan[i] << std::endl;
        //std::cout << "Hit charge: " << hit_charge[i] << std::endl;


        uint64_t tp_t_begin =                                                        // NOLINT(build/unsigned)
          timestamp + clocksPerTPCTick * (int64_t(hit_end[i]) - int64_t(hit_tover[i]));       // NOLINT(build/unsigned)
        uint64_t tp_t_end = timestamp + clocksPerTPCTick * int64_t(hit_end[i]);      // NOLINT(build/unsigned)
        // May be needed for TPSet:
        // uint64_t tspan = clocksPerTPCTick * hit_tover[i]; // is/will be this needed?
        //
        // For quick n' dirty debugging: print out time/channel of hits.
        // Can then make a text file suitable for numpy plotting with, eg:
        //
        //
        //TLOG_DEBUG(0) << "Hit: " << tp_t_begin << " " << offline_channel;
        triggeralgs::TriggerPrimitive trigprim;
        trigprim.time_start = tp_t_begin;
        trigprim.time_peak = (tp_t_begin + tp_t_end) / 2;
        trigprim.time_over_threshold = hit_tover[i] * clocksPerTPCTick;      
        trigprim.channel = register_channels[chan[i]]; //offline channel map
        trigprim.adc_integral = hit_charge[i];
        trigprim.adc_peak = hit_charge[i] / 20;
        trigprim.detid = 666;          
        trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
        trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
        trigprim.version = 1;
        if (save_trigprim){
          save_TP_object(trigprim, "AVX");
        }          

        ++total_hits;
      }
    } // loop over 16 registers 
  } // while not magic   

}


void extract_hits_naive(uint16_t* output_location, uint64_t timestamp,
                      std::array<uint, swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER>& register_channels,
                      unsigned int& total_hits,
                      bool save_trigprim) {

    constexpr int clocksPerTPCTick = 32;
    uint16_t chan, hit_end, hit_charge, hit_tover; 

    size_t i = 0;
    while (*output_location != swtpg_wibeth::MAGIC) {
      chan   = *output_location++;
      hit_end    = *output_location++;
      hit_charge  = *output_location++;
      hit_tover     = *output_location++;
    
      i += 1;
      uint64_t tp_t_begin =                                                        
        timestamp + clocksPerTPCTick * (int64_t(hit_end ) - hit_tover );       
      uint64_t tp_t_end = timestamp + clocksPerTPCTick * int64_t(hit_end );      

      triggeralgs::TriggerPrimitive trigprim;
      trigprim.time_start = tp_t_begin;
      trigprim.time_peak = (tp_t_begin + tp_t_end) / 2;

      trigprim.time_over_threshold = hit_tover  * clocksPerTPCTick;


      trigprim.channel = register_channels[chan]; //offline channel map;
      trigprim.adc_integral = hit_charge ;
      trigprim.adc_peak = hit_charge  / 20;
      trigprim.detid = 666; 
      trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
      trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
      trigprim.version = 1;
    
      if (save_trigprim) {
        save_TP_object(trigprim, "NAIVE");
      }
      ++total_hits;

    }
}




// =================================================================
//                       EMULATOR
// =================================================================


class tpg_emulator {

public:

  tpg_emulator(bool save_adc_data, bool save_trigprim, bool parse_trigger_primitive, std::string select_algorithm, std::string select_channel_map) {

    m_save_adc_data = save_adc_data;
    m_save_trigprim = save_trigprim;
    m_parse_trigger_primitive = parse_trigger_primitive;
    m_select_algorithm = select_algorithm;
    m_select_channel_map = select_channel_map;

  }

  void set_tpg_threshold(int tpg_threshold){
    m_tpg_threshold = tpg_threshold;
  }

  int get_total_hit_number () {
    return m_total_hits;
  } 
  int get_total_hits_trigger_primitive () {
    return m_total_hits_trigger_primitive;
  }   

  void set_CPU_affinity(int core_number) {
    m_CPU_core = core_number;
  }

  void initialize() {

    if (m_select_algorithm == "SimpleThreshold") {
      m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_avx2<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
    } else if (m_select_algorithm == "AbsRS") {
      m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_rs_avx2<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
    } else {
      throw tpgtools::TPGAlgorithmInexistent(ERS_HERE, m_select_algorithm);     
    }

    // Initialize the channel map if a valid name has been selected
    if (m_select_channel_map != "None") {
      TLOG() << "Using channel map: " << m_select_channel_map;
      m_channel_map = dunedaq::detchannelmaps::make_map(m_select_channel_map);
    } else {
      TLOG() << "*** No channel map has been provided. " ;
    }

    // Initialize frame handler
    m_frame_handler.initialize(m_tpg_threshold);

  }


  void execute_tpg(const dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter* fp) {

  // Set CPU affinity of the TPG thread
  SetAffinityThread(m_CPU_core);

  // Parse the WIBEth frames
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>((uint8_t*)fp);
  uint64_t timestamp = wfptr->get_timestamp();     

  // Frame expansion
  swtpg_wibeth::MessageRegisters registers_array;
  swtpg_wibeth::expand_wibeth_adcs(fp, &registers_array);

  
  if (m_first_hit) {   
    m_frame_handler.m_tpg_processing_info->setState(registers_array);
    m_first_hit = false;    
    // Save ADC info
    if (m_save_adc_data){
      save_raw_data(registers_array, timestamp, -1, m_select_algorithm );
    }
    
  }

  m_frame_handler.m_tpg_processing_info->input = &registers_array;
  uint16_t* destination_ptr = m_frame_handler.get_hits_dest();
  *destination_ptr = swtpg_wibeth::MAGIC;
  m_frame_handler.m_tpg_processing_info->output = destination_ptr;
  m_assigned_tpg_algorithm_function(*m_frame_handler.m_tpg_processing_info);
       

  // Parse the output from the TPG    
  extract_hits_avx(m_frame_handler.m_tpg_processing_info->output, timestamp, m_register_channels, m_total_hits, m_save_trigprim);

}


  void process_fragment(std::unique_ptr<dunedaq::daqdataformats::Fragment>&& frag_ptr,
                        int& record_idx_TR, int& record_idx_TP) {

    uint32_t element_id = 0;
    if (frag_ptr->get_fragment_type() == dunedaq::daqdataformats::FragmentType::kWIBEth) {
      
      element_id = frag_ptr->get_element_id().id;
      int num_frames =
        (frag_ptr->get_size() - sizeof(dunedaq::daqdataformats::FragmentHeader)) / sizeof(dunedaq::fddetdataformats::WIBEthFrame);                
      TLOG_DEBUG(TLVL_BOOKKEEPING) << "Trigger Record number " << record_idx_TR << " has "   << num_frames << " frames" ;
      
      for (int i = 0; i < num_frames; ++i) {
        // Read the Trigger Record data as a WIBEth frame
        auto fr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(
          static_cast<char*>(frag_ptr->get_data()) + i * sizeof(dunedaq::fddetdataformats::WIBEthFrame)
        );

        // Execute the TPG algorithm on the WIBEth adapter frames
        auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter*>(fr);


        // Register the offline channel numbers
        // AAA: TODO: find a more elegant way of register the channel map
        if (m_select_channel_map != "None") {
          m_register_channel_map = swtpg_wibeth::get_register_to_offline_channel_map_wibeth(fr, m_channel_map);             
          for (size_t i = 0; i < swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER; ++i) {
            m_register_channels[i] = m_register_channel_map.channel[i];                
          }
        } else {
          // If no channel map is not selected use the values from 0 to 63
          std::iota(m_register_channels.begin(), m_register_channels.end(), 0);  
        }
  
        execute_tpg(fp);


      } // end loop over number of frames      
      // Finished processing all the frames for the given WIBEth fragment
      ++record_idx_TR;
    } // if trigger record is WIBEth type

    else if (frag_ptr->get_fragment_type() == dunedaq::daqdataformats::FragmentType::kTriggerPrimitive) {
      // Parse only the Trigger Primitives with the same ID of the ones with data frames
      // AAA: NOT SURE OF THIS STATEMENT!! (INTRODUCED TO AVOID SAME TPs from multiple trigger id values)
      if (frag_ptr->get_element_id().id == element_id) {
        if (m_parse_trigger_primitive) {
          process_trigger_primitive(std::move(frag_ptr), record_idx_TP, m_total_hits_trigger_primitive, m_save_trigprim);
          ++record_idx_TP;
        }
      }                  
    }
 


  }

private: 
 
  bool m_save_adc_data = false; 
  bool m_save_trigprim = false;  
  bool m_parse_trigger_primitive = false;

  std::string m_select_algorithm = "";
  std::string m_select_channel_map = "None";

  dunedaq::fdreadoutlibs::WIBEthFrameHandler m_frame_handler;


  std::function<void(swtpg_wibeth::ProcessingInfo<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>& info)> m_assigned_tpg_algorithm_function;
  
  // Channel mapping
  std::shared_ptr<dunedaq::detchannelmaps::TPCChannelMap> m_channel_map;  
  // Map from expanded AVX register position to offline channel number
  swtpg_wibeth::RegisterChannelMap m_register_channel_map;   
  // Mapping from expanded AVX register position to offline channel number
  std::array<uint, swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER> m_register_channels = {};


  unsigned int m_total_hits = 0;
  unsigned int m_total_hits_trigger_primitive = 0;
  bool m_first_hit = true;


  int m_tpg_threshold = 500; //default value 
  int m_CPU_core = 0;

  
}; 

#endif // TPGTOOLS_INCLUDE_TPGTOOLS_TPGEMULATOR_HPP_