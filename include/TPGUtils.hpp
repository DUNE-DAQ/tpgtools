/**
 * @file TPGUtils.hpp 
 * Utility functions for the TPG emulator
  *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef TPGTOOLS_INCLUDE_TPGTOOLS_TPGUTILS_HPP_
#define TPGTOOLS_INCLUDE_TPGTOOLS_TPGUTILS_HPP_



#include "fddetdataformats/WIBEthFrame.hpp"


#include "readoutlibs/utils/FileSourceBuffer.hpp"
#include "readoutlibs/utils/BufferedFileWriter.hpp"


#include "fdreadoutlibs/DUNEWIBEthTypeAdapter.hpp"



#include "triggeralgs/TriggerPrimitive.hpp"
#include "trgdataformats/TriggerPrimitive.hpp"
#include "hdf5libs/HDF5RawDataFile.hpp"
#include "logging/Logging.hpp"

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
void save_TP_object( triggeralgs::TriggerPrimitive trigprim, std::string algo, 
		     std::string out_suffix ){
  std::ofstream out_file; 

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);  
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M");
  auto date_time_str = oss.str();

  std::string file_name = "TP_dump_" + algo + "_" + date_time_str + out_suffix + ".txt";
  out_file.open(file_name.c_str(), std::ofstream::app);

  //offline channel, start time, time over threshold [ns], peak_time, ADC sum, amplitude    
  out_file << trigprim.channel << "," << trigprim.time_start << "," << trigprim.time_over_threshold << "," 
	   << trigprim.time_peak << "," << trigprim.adc_integral << ","  << trigprim.adc_peak <<  ","  << trigprim.type << "\n";  

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


#endif // TPGTOOLS_INCLUDE_TPGTOOLS_TPGUTILS_HPP_

