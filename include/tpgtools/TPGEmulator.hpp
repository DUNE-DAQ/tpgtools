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
#include "TPGUtils.hpp"


using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;





// =================================================================
//                       EMULATOR
// =================================================================





class tpg_emulator_base {

public:

  tpg_emulator_base(bool save_adc_data, bool save_trigprim, bool parse_trigger_primitive, std::string select_algorithm, std::string select_channel_map) {

    m_save_adc_data = save_adc_data;
    m_save_trigprim = save_trigprim;
    m_parse_trigger_primitive = parse_trigger_primitive;
    m_select_algorithm = select_algorithm;
    m_select_channel_map = select_channel_map;

  }

  void register_channel_map(const dunedaq::fddetdataformats::WIBEthFrame* frame) {
    // Register the offline channel numbers
    // AAA: TODO: find a more elegant way of register the channel map
    if (m_select_channel_map != "None") {
      m_register_channel_map = swtpg_wibeth::get_register_to_offline_channel_map_wibeth(frame, m_channel_map);             
      for (size_t i = 0; i < swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER; ++i) {
        m_register_channels[i] = m_register_channel_map.channel[i];                
      }
    } else {
      // If no channel map is not selected use the values from 0 to 63
      std::iota(m_register_channels.begin(), m_register_channels.end(), 0);  
    }
  }

  void set_tpg_threshold(int tpg_threshold){
    m_tpg_threshold = tpg_threshold;
  }

  void set_CPU_affinity(int core_number) {
    m_CPU_core = core_number;
  }  

  unsigned int get_total_hits() {
    return m_total_hits;
  }


  bool m_save_adc_data = false; 
  bool m_save_trigprim = false;  
  bool m_parse_trigger_primitive = false;

  std::string m_select_algorithm = "";
  std::string m_select_channel_map = "None";


  unsigned int m_total_hits = 0;
  unsigned int m_total_hits_trigger_primitive = 0;
  bool m_first_hit = true;


  int m_tpg_threshold = 500; //default value 
  int m_CPU_core = 0;

  // Frame Handler 
  dunedaq::fdreadoutlibs::WIBEthFrameHandler m_frame_handler;

  // TPG algorithm function
  std::function<void(swtpg_wibeth::ProcessingInfo<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>& info)> m_assigned_tpg_algorithm_function;

  // Channel mapping
  std::shared_ptr<dunedaq::detchannelmaps::TPCChannelMap> m_channel_map;  
  // Map from expanded AVX register position to offline channel number
  swtpg_wibeth::RegisterChannelMap m_register_channel_map;   
  // Mapping from expanded AVX register position to offline channel number
  std::array<uint, swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER> m_register_channels = {};



};




class tpg_emulator_avx : public tpg_emulator_base {

public: 

// Inheriting the base class constructor
using tpg_emulator_base::tpg_emulator_base;

void extract_hits(uint16_t* output_location, uint64_t timestamp,
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
        trigprim.channel = m_register_channels[chan[i]]; //offline channel map
        trigprim.adc_integral = hit_charge[i];
        trigprim.adc_peak = hit_charge[i] / 20;
        trigprim.detid = 666;          
        trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
        trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
        trigprim.version = 1;
        if (save_trigprim){
          save_TP_object(trigprim, "AVX");
        }          

        ++m_total_hits;
      }
    } // loop over 16 registers 
  } // while not magic   

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
  extract_hits(m_frame_handler.m_tpg_processing_info->output, timestamp, m_save_trigprim);

 }

 void initialize()  {
    
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

  };





};






class tpg_emulator_naive : public tpg_emulator_base {

public: 

void extract_hits(uint16_t* output_location, uint64_t timestamp,
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


      trigprim.channel = m_register_channels[chan]; //offline channel map;
      trigprim.adc_integral = hit_charge ;
      trigprim.adc_peak = hit_charge  / 20;
      trigprim.detid = 666; 
      trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
      trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
      trigprim.version = 1;
    
      if (save_trigprim) {
        save_TP_object(trigprim, "NAIVE");
      }
      ++m_total_hits;

    }

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
  extract_hits(m_frame_handler.m_tpg_processing_info->output, timestamp, m_save_trigprim);

  }

  void initialize()  {
    
    if (m_select_algorithm == "SimpleThreshold") {
      m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_naive<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
    } else if (m_select_algorithm == "AbsRS") {
      m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_naive_RS<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
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

  };



};


#endif // TPGTOOLS_INCLUDE_TPGTOOLS_TPGEMULATOR_HPP_