/**
 * @file ReplayTPG.cxx 
 * Replay tpg algorithms on trigger records
 * 
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

// DUNE-DAQ
#include "fddetdataformats/WIBEthFrame.hpp"
#include "fddetdataformats/WIBEthFrame.hpp"
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

#include "tpgtools/TPGEmulator.hpp"

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

using namespace dunedaq::hdf5libs;

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;



// =================================================================
//                       MAIN
// =================================================================
int 
main(int argc, char** argv)
{

    CLI::App app{ "TPG algorithms emulator" };

    // Set default input frame file
    std::string file_path_input = "";
    app.add_option("-f,--file-path-input", file_path_input, "Path to the input file");

    std::string select_algorithm = "";
    app.add_option("-a,--algorithm", select_algorithm, "TPG Algorithm (SimpleThreshold / AbsRS)");
  
    std::string select_channel_map = "None";
    app.add_option("-m,--channel-map", select_channel_map, "Select a valid channel map: None, VDColdboxChannelMap, ProtoDUNESP1ChannelMap, PD2HDChannelMap, HDColdboxChannelMap, FiftyLChannelMap");

    int num_TR_to_read = -1;
    app.add_option("-n,--num-TR-to-read", num_TR_to_read, "Number of Trigger Records to read. Default: select all TRs.");

    int tpg_threshold = 500;
    app.add_option("-t,--tpg-threshold", tpg_threshold, "Value of the TPG threshold. Default value is 500.");

    int core_number = 0;
    app.add_option("-c,--core", core_number, "Set core number of the executing TPG thread. Default value is 0.");

    bool save_adc_data = false;
    app.add_flag("--save-adc-data", save_adc_data, "Save ADC data (first frame only)");

    bool save_trigprim = false;
    app.add_flag("--save-trigprim", save_trigprim, "Save trigger primitive data");

    bool parse_trigger_primitive = false;
    app.add_flag("--parse_trigger_primitive", parse_trigger_primitive, "Parse Trigger Primitive records");


    CLI11_PARSE(app, argc, argv);

    // =================================================================
    //                       Setup the TPG emulator
    // =================================================================

    // AAA: maybe this constructor should accept a config file...
    tpg_emulator emulator(save_adc_data, save_trigprim, parse_trigger_primitive, select_algorithm, select_channel_map) ;
    emulator.set_tpg_threshold(tpg_threshold);
    emulator.set_CPU_affinity(core_number);
    emulator.initialize();

    
    // =================================================================
    //                       READ THE HDF5 FILE
    // =================================================================

    // open our file reading
    const std::string ifile_name = file_path_input;
    if (ifile_name.empty()) {
      throw tpgtools::FileInexistent(ERS_HERE, ifile_name);    
    } 
    HDF5RawDataFile h5_raw_data_file(ifile_name);

    auto recorded_size = h5_raw_data_file.get_attribute<size_t>("recorded_size");
    auto record_type = h5_raw_data_file.get_record_type();

    auto run_number = h5_raw_data_file.get_attribute<unsigned int>("run_number");
    //auto file_index = h5_raw_data_file.get_attribute<unsigned int>("file_index");
  
    auto creation_timestamp = h5_raw_data_file.get_attribute<std::string>("creation_timestamp");
    auto app_name = h5_raw_data_file.get_attribute<std::string>("application_name");


    std::cout << "Run number: " << run_number << std::endl; 
    std::cout << "Recorded size [bytes]: " << recorded_size << std::endl; 
    std::cout << "Recorded type: " << record_type << std::endl; 


    auto records = h5_raw_data_file.get_all_record_ids();

    if (records.empty()) {
      TLOG() << "*** NO TRIGGER RECORDS FOUND" ;    
      return 0;
    }

    int record_idx_TR = 0;
    int record_idx_TP = 0;

    // Measure the time taken to read the whole trigger record file
    auto start_test = std::chrono::high_resolution_clock::now();  

   
    for (auto const& rid : records) {
      for(auto const& frag_dataset : h5_raw_data_file.get_fragment_dataset_paths(rid)) {
        auto frag_ptr = h5_raw_data_file.get_frag_ptr(frag_dataset);

        if (record_idx_TR <= num_TR_to_read || num_TR_to_read == -1 ) {  

          emulator.process_fragment(std::move(frag_ptr), record_idx_TR, record_idx_TP);
        
        } 
        
        
      } // loop over all the fragments in a single trigger record    
    } // loop over all trigger records

    // Calculate elapsed time in seconds 
    // AAA: to be removed if not useufl?
    auto now = std::chrono::high_resolution_clock::now();
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_test).count();  



    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Elapsed time for reading input file [ms]: " << elapsed_milliseconds;
    std::cout << "Found in total " << emulator.get_total_hit_number() << " hits" << std::endl;

    if (parse_trigger_primitive) {
      std::cout << "Found in total  (from Trigger Primitive objects) " << emulator.get_total_hits_trigger_primitive() << " TPs" << std::endl;
    }
    


}

