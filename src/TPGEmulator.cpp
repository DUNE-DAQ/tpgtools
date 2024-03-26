/**
 * @file TPGEmulator.hpp 
 * Emulator classes for TPG applications
  *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "tpgtools/TPGEmulator.hpp"

#include "fdreadoutlibs/wibeth/tpg/ProcessAVX2.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessAbsRSAVX2.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessStandardRSAVX2.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessNaive.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessNaiveRS.hpp"


// =================================================================
//                       tpg_emulator_base
// =================================================================
// Function to save raw ADC data to a file (only for debugging) 
void tpg_emulator_base::save_raw_data(swtpg_wibeth::MessageRegisters register_array, 
	       uint64_t t0, int channel_number, std::string algo)
{
  std::ofstream out_file;
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);  
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M");
  auto date_time_str = oss.str();

  std::string file_name;
  std::string file_name_TR_info = "";
  if (m_register_TR_record_idx != -1 && m_register_TR_frame_idx != -1) {
    file_name_TR_info = "_tr_" + std::to_string(m_register_TR_record_idx) + "_fr_" + std::to_string(m_register_TR_frame_idx);
  }
  if (channel_number == -1) {
    file_name = "all_channels_" + algo + file_name_TR_info + "_data_" + date_time_str + ".txt";
  } else {
    file_name = "Channel_" + std::to_string(channel_number) + "_" + algo + file_name_TR_info + "_data_" + date_time_str + ".txt";
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





// =================================================================
//                       tpg_emulator_avx
// =================================================================

void tpg_emulator_avx::extract_hits(uint16_t* output_location, uint64_t timestamp,
                      bool save_trigprim, std::string out_suffix="") {

  //constexpr int clocksPerTPCTick = 32;
  //uint16_t chan[16], hit_end[16], hit_charge[16], hit_tover[16];
  constexpr int clocksPerTPCTick = dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter::samples_tick_difference;
  const constexpr std::size_t nreg = swtpg_wibeth::SAMPLES_PER_REGISTER;
  uint16_t chan[nreg], left[nreg], hit_end[nreg], hit_peak_adc[nreg], hit_charge[nreg], hit_tover[nreg], hit_peak_time[nreg];

  while (*output_location != swtpg_wibeth::MAGIC) {
    for (std::size_t i = 0; i < nreg; ++i) {
      chan[i] = *output_location++;
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_end[i] = *output_location++;
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_charge[i] = *output_location++;
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_tover[i] = *output_location++;
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_peak_adc[i] = *output_location++;
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_peak_time[i] = *output_location++;
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      left[i] = *output_location++;
    }    
    
    // Now that we have all the register values in local
    // variables, loop over the register index (ie, channel) and
    // find the channels which actually had a hit, as indicated by
    // nonzero value of hit_charge
    for (std::size_t i = 0; i < nreg; ++i) {
      if (hit_charge[i] && left[i] == std::numeric_limits<std::uint16_t>::max()
          && chan[i] != swtpg_wibeth::MAGIC) {

        uint64_t tp_t_begin = timestamp + clocksPerTPCTick * ((int64_t)hit_end[i] - (int64_t)hit_tover[i]); // NOLINT(build/unsigned)
        uint64_t tp_t_peak  = tp_t_begin + clocksPerTPCTick * hit_peak_time[i]; // NOLINT(build/unsigned)

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
        trigprim.time_peak = tp_t_peak;
        trigprim.time_over_threshold = uint64_t((hit_tover[i] - 1) * clocksPerTPCTick);
        trigprim.channel = m_register_channels[chan[i]];
        trigprim.adc_integral = hit_charge[i];
        trigprim.adc_peak = hit_peak_adc[i];
        trigprim.detid = 666;
        trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
        trigprim.algorithm = m_tp_algo;
        trigprim.version = 1;
        if (save_trigprim){
          save_TP_object(trigprim, "AVX", out_suffix);
        }          

        ++m_total_hits;
      }
    } // loop over 16 registers 
  } // while not magic   

}

void tpg_emulator_avx::execute_tpg(const dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter* fp) {

  // Set CPU affinity of the TPG thread
  SetAffinityThread(m_CPU_core);

  // Parse the WIBEth frames
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>((uint8_t*)fp);
  uint64_t timestamp = wfptr->get_timestamp();     

  // Frame expansion
  swtpg_wibeth::MessageRegisters registers_array;
  swtpg_wibeth::expand_wibeth_adcs(fp, &registers_array);

  
  if (m_first_hit) {   
    m_frame_handler.m_tpg_processing_info->setState(registers_array, m_register_memory_factor);
    m_first_hit = false;    
    // Save ADC info
    if (m_save_adc_data && m_num_frames_to_save == 1) {
      save_raw_data(registers_array, timestamp, -1, m_select_algorithm);
    } 
  }
  if (m_save_adc_data && m_num_frames_to_save == -1) {
    save_raw_data(registers_array, timestamp, -1, m_select_algorithm);
  }

  m_frame_handler.m_tpg_processing_info->input = &registers_array;
  uint16_t* destination_ptr = m_frame_handler.get_hits_dest();
  *destination_ptr = swtpg_wibeth::MAGIC;
  m_frame_handler.m_tpg_processing_info->output = destination_ptr;
  m_assigned_tpg_algorithm_function(*m_frame_handler.m_tpg_processing_info);
       

  // Parse the output from the TPG    
  extract_hits(m_frame_handler.m_tpg_processing_info->output, timestamp, m_save_trigprim, m_out_suffix);

 }

 void tpg_emulator_avx::initialize()  {
    
    if (m_select_algorithm == "SimpleThreshold") {
      m_tp_algo = dunedaq::trgdataformats::TriggerPrimitive::Algorithm::kSimpleThreshold;
      m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_avx2<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;      
    } else if (m_select_algorithm == "AbsRS") {
      m_tp_algo = dunedaq::trgdataformats::TriggerPrimitive::Algorithm::kAbsRunningSum;
      m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_rs_avx2<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
    } else if (m_select_algorithm == "StandardRS") {
      m_tp_algo = dunedaq::trgdataformats::TriggerPrimitive::Algorithm::kRunningSum;
      m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_standard_rs_avx2<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
    } else {
      throw tpgtools::TPGAlgorithmInexistent(ERS_HERE, m_select_algorithm);     
    }

    // Initialize the channel map if a valid name has been selected
    if (!m_select_channel_map.empty()) {
      TLOG() << "Using channel map: " << m_select_channel_map;
      m_channel_map = dunedaq::detchannelmaps::make_map(m_select_channel_map);
    }

    for (size_t i = 0; i < swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER; ++i) {
      m_register_memory_factor[i] = m_tpg_rs_memory_factor;
    }

    // Initialize frame handler
    m_frame_handler.initialize(m_tpg_threshold, m_tpg_rs_memory_factor, m_tpg_rs_scale_factor, m_tpg_frugal_streaming_accumulator_limit);

};


// =================================================================
//                       tpg_emulator_naive
// =================================================================


void tpg_emulator_naive::extract_hits(uint16_t* output_location, uint64_t timestamp,
                      bool save_trigprim, std::string out_suffix="") {


    constexpr int clocksPerTPCTick = dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter::samples_tick_difference;
    const constexpr std::size_t nreg = swtpg_wibeth::SAMPLES_PER_REGISTER;
    uint16_t chan, hit_end, hit_peak_adc, hit_charge, hit_tover, hit_peak_time;

    std::array<int, 16> indices{0, 1, 2, 3, 4, 5, 6, 7, 15, 8, 9, 10, 11, 12, 13, 14};

    size_t i = 0;
    while (*output_location != swtpg_wibeth::MAGIC) {
      chan            = *output_location++;
      hit_end         = *output_location++;
      hit_charge      = *output_location++;
      hit_tover       = *output_location++;
      hit_peak_adc    = *output_location++;
      hit_peak_time   = *output_location++;

      i += 1;
      chan = nreg*(chan/nreg)+indices[chan%nreg];

      uint64_t tp_t_begin = timestamp + clocksPerTPCTick * ((int64_t)hit_end - (int64_t)hit_tover); // NOLINT(build/unsigned)
      uint64_t tp_t_peak  = tp_t_begin + clocksPerTPCTick * hit_peak_time; // NOLINT(build/unsigned)

      triggeralgs::TriggerPrimitive trigprim;
      trigprim.time_start = tp_t_begin;
      trigprim.time_peak = tp_t_peak;

      trigprim.time_over_threshold = uint64_t((hit_tover - 1) * clocksPerTPCTick);


      trigprim.channel = m_register_channels[chan]; //offline channel map;
      trigprim.adc_integral = hit_charge ;
      trigprim.adc_peak = hit_peak_adc;
      trigprim.detid = 666; 
      trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
      trigprim.algorithm = m_tp_algo;
      trigprim.version = 1;
    
      if (save_trigprim) {
        save_TP_object(trigprim, "NAIVE", out_suffix);
      }
      ++m_total_hits;

    }

}

void tpg_emulator_naive::execute_tpg(const dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter* fp) {

  // Set CPU affinity of the TPG thread
  SetAffinityThread(m_CPU_core);

  // Parse the WIBEth frames
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>((uint8_t*)fp);
  uint64_t timestamp = wfptr->get_timestamp();     

  // Frame expansion
  swtpg_wibeth::MessageRegisters registers_array;
  swtpg_wibeth::expand_wibeth_adcs(fp, &registers_array);

  
  if (m_first_hit) {   
    m_frame_handler.m_tpg_processing_info->setState(registers_array, m_register_memory_factor);
    m_first_hit = false;    
    // Save ADC info
    if (m_save_adc_data && m_num_frames_to_save == 1) {
      save_raw_data(registers_array, timestamp, -1, m_select_algorithm);
    }
    
  }
  if (m_save_adc_data && m_num_frames_to_save == -1) {
    save_raw_data(registers_array, timestamp, -1, m_select_algorithm);
  }

  m_frame_handler.m_tpg_processing_info->input = &registers_array;
  uint16_t* destination_ptr = m_frame_handler.get_hits_dest();
  *destination_ptr = swtpg_wibeth::MAGIC;
  m_frame_handler.m_tpg_processing_info->output = destination_ptr;
  m_assigned_tpg_algorithm_function(*m_frame_handler.m_tpg_processing_info);
       

  // Parse the output from the TPG    
  extract_hits(m_frame_handler.m_tpg_processing_info->output, timestamp, m_save_trigprim, m_out_suffix);

}



void tpg_emulator_naive::initialize()  {
  
  if (m_select_algorithm == "SimpleThreshold") {
    m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_naive<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
  } else if (m_select_algorithm == "AbsRS") {
    m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_naive_RS<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
  } else {
    throw tpgtools::TPGAlgorithmInexistent(ERS_HERE, m_select_algorithm);     
  }
  // Initialize the channel map if a valid name has been selected
  if (!m_select_channel_map.empty()) {
    TLOG() << "Using channel map: " << m_select_channel_map;
    m_channel_map = dunedaq::detchannelmaps::make_map(m_select_channel_map);
  } 

  for (size_t i = 0; i < swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER; ++i) {
    m_register_memory_factor[i] = m_tpg_rs_memory_factor;
  }
   
  // Initialize frame handler
  // AAA: fix me with proper TPG configuration parameters
  m_frame_handler.initialize(m_tpg_threshold, m_tpg_rs_memory_factor, m_tpg_rs_scale_factor, m_tpg_frugal_streaming_accumulator_limit);


};
