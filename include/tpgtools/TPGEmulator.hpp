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


#include "fdreadoutlibs/wibeth/WIBEthFrameProcessor.hpp"

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
  tpg_emulator_base(tpg_emulator_base const&) = delete;
  tpg_emulator_base(tpg_emulator_base&&) = delete;
  tpg_emulator_base& operator=(tpg_emulator_base const&) = delete;
  tpg_emulator_base& operator=(tpg_emulator_base&&) = delete;

  virtual ~tpg_emulator_base() = default;


  virtual void execute_tpg(const dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter* /*fp*/) {};  
  virtual void initialize() {};

  void register_channel_map(const dunedaq::fddetdataformats::WIBEthFrame* frame) {
    // Register the offline channel numbers
    // AAA: TODO: find a more elegant way of register the channel map
    if (!m_select_channel_map.empty()) {
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
  std::string m_select_channel_map = "";


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
                      bool save_trigprim); 

void execute_tpg(const dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter* fp);

void initialize();


};






class tpg_emulator_naive : public tpg_emulator_base {

public: 

// Inheriting the base class constructor
using tpg_emulator_base::tpg_emulator_base;


void extract_hits(uint16_t* output_location, uint64_t timestamp,
                      bool save_trigprim);

void execute_tpg(const dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter* fp);

void initialize();



};


#endif // TPGTOOLS_INCLUDE_TPGTOOLS_TPGEMULATOR_HPP_