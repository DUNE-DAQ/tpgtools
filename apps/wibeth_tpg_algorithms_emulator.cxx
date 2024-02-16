/**
 * @file ReplayTPG.cxx 
 * Replay tpg algorithms on trigger records
 * 
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "tpgtools/TPGEmulator.hpp"

// system
#include "CLI/CLI.hpp"
#include <iostream>
#include <fstream>

#include <fmt/core.h>
#include <fmt/format.h>

using namespace dunedaq::hdf5libs;

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;

 
// =================================================================
//                       MAIN
// =================================================================
int 
main(int argc, char** argv)
{

    CLI::App app{ "TPG algorithms emulator using input from Trigger Record files" };

    std::string file_path_input = "";
    app.add_option("-f,--file-path-input", file_path_input, "Path to the input file");

    std::string select_algorithm = "SimpleThreshold";
    app.add_option("-a,--algorithm", select_algorithm, "TPG Algorithm (SimpleThreshold / AbsRS)");

    std::string select_implementation = "AVX";
    app.add_option("-i,--implementation", select_implementation, "TPG implementation (AVX / NAIVE). Default: AVX");

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

    int num_frames_to_save = 1;
    app.add_option("-s,--num-frames-to-save", num_frames_to_save, "Set the number of frames of ADC data from a TR to save: -1 (all) or 1. Default: 1 (first frame only).");

    CLI11_PARSE(app, argc, argv);

    // =================================================================
    //                       Setup the TPG emulator
    // =================================================================


    // Create instance of the TPG emulator implementation    
    std::unique_ptr<tpg_emulator_base> emulator;
    if (select_implementation == "AVX") {
      emulator = std::make_unique<tpg_emulator_avx>(save_adc_data, save_trigprim, parse_trigger_primitive, select_algorithm, select_channel_map);
    } else if (select_implementation == "NAIVE") {
      emulator = std::make_unique<tpg_emulator_naive>(save_adc_data, save_trigprim, parse_trigger_primitive, select_algorithm, select_channel_map);
    } else {
      throw tpgtools::InvalidImplementation(ERS_HERE, select_implementation);  
    }

    emulator->set_tpg_threshold(tpg_threshold);
    emulator->set_CPU_affinity(core_number);
    emulator->set_num_frames_to_save(num_frames_to_save);
    emulator->initialize();

    
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


    fmt::print("Run number {} \n", run_number);
    fmt::print("Recorded size [bytes]: {} \n", recorded_size);
    fmt::print("Recorded type: {} \n", record_type);


    auto records = h5_raw_data_file.get_all_record_ids();

    if (records.empty()) {
      fmt::print("*** NO TRIGGER RECORDS FOUND \n");  
      return 0;
    }

    int record_idx_TR = 0;
    int record_idx_TP = 0;

    unsigned int number_TPs_from_TR = 0;

    // Measure the time taken to read the whole trigger record file
    auto start_test = std::chrono::high_resolution_clock::now();  

   
    for (auto const& rid : records) {
      for(auto const& frag_dataset : h5_raw_data_file.get_fragment_dataset_paths(rid)) {
        auto frag_ptr = h5_raw_data_file.get_frag_ptr(frag_dataset);

        if (record_idx_TR <= num_TR_to_read || num_TR_to_read == -1 ) {  

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
  
              emulator->register_channel_map(fr);
      
              // Execute the TPG algorithm on the WIBEth adapter frames
              auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter*>(fr);
              
              if (num_frames_to_save == -1) { 
	        emulator->register_TR_info(record_idx_TR, i);
              }		
              emulator->execute_tpg(fp);
      
      
            } // end loop over number of frames      
            // Finished processing all the frames for the given WIBEth fragment
            ++record_idx_TR;
          } // if trigger record is WIBEth type
          else if (frag_ptr->get_fragment_type() == dunedaq::daqdataformats::FragmentType::kTriggerPrimitive) {
            // Parse only the Trigger Primitives with the same ID of the ones with data frames
            // AAA: NOT SURE OF THIS STATEMENT!! (INTRODUCED TO AVOID SAME TPs from multiple trigger id values)
            if (frag_ptr->get_element_id().id == element_id) {
              if (parse_trigger_primitive) {
                process_trigger_primitive(std::move(frag_ptr), record_idx_TP, number_TPs_from_TR, save_trigprim);
                ++record_idx_TP;
              }
            }                  
          } // if trigger record is TriggerPrimitive type
        }
        
        
      } // loop over all the fragments in a single trigger record    
    } // loop over all trigger records

    // Calculate elapsed time in seconds 
    // AAA: to be removed if not useufl?
    auto now = std::chrono::high_resolution_clock::now();
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_test).count();  



    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Elapsed time for reading input file [ms]: " << elapsed_milliseconds;
    fmt::print("Found in total {} hits \n", emulator->get_total_hits());


    if (parse_trigger_primitive) {
      fmt::print("Found in total (from Trigger Record) {} TPs \n", number_TPs_from_TR);
    }
    


}

