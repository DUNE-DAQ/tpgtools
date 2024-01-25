/**
 * @file ReplayTPG.cxx 
 * Replay tpg algorithms on trigger records
 * 
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "tpgtools/TPGPatternGenerator.hpp"

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

    CLI::App app{ "TPG pattern generator using as input and storing output as raw WIBEth ADC binary files" };

    std::string file_path_input = "";
    app.add_option("-f,--file-path-input", file_path_input, "Path to the input file");

    int num_frames_to_read = -1;
    app.add_option("-n,--num-frames-to-read", num_frames_to_read, "Number of WIBEth frames to read. Default: select all frames.");

    size_t input_ch = 0;
    app.add_option("-i,--input_channel", input_ch, "Input channel number for adding fake hit. Default: 0");

    int tpg_threshold = 500;
    app.add_option("-t,--tpg-threshold", tpg_threshold, "Value of the TPG threshold. Default value is 500.");

    bool save_adc_data = false;
    app.add_flag("--save-adc-data", save_adc_data, "Save ADC data (first frame only)");

    bool save_trigprim = false;
    app.add_flag("--save-trigprim", save_trigprim, "Save trigger primitive data");

    int time_tick_offset = 1;
    app.add_option("-o,--time-tick-offset", time_tick_offset, "Time tick of pattern start. Default: 1 (max:63)");

    std::string out_suffix = "";
    app.add_option("-s ,--out_suffix", out_suffix, "Append string (suffix) to output hit file name");

    std::string select_pattern = "patt_golden";
    app.add_option("-p,--select-pattern", select_pattern, "Test pattern name (patt_golden, patt_pulse, patt_edge_square, patt_edge_left, patt_edge_right). Default: patt_golden");

    CLI11_PARSE(app, argc, argv);

    // init
    dunedaq::fdreadoutlibs::WIBEthFrameHandler fh;
    PattgenHandler ph;
    PattgenAlgs pa(ph.m_pattgen_info);
    uint64_t first_timestamp = 0;

    // =================================================================
    //                       READ THE FRAMES.BIN FILE
    // =================================================================
    FrameFile input_file = FrameFile(file_path_input.c_str());
    int total_num_frames = input_file.num_frames();

    std::cout << "Input file name: " << file_path_input << std::endl;
    std::cout << "Size of the input file: " << input_file.length() << std::endl;
    std::cout << "Input channel number: " << input_ch << std::endl;
    std::cout << "Number of DUNE WIBEth frames in the input file: " << total_num_frames << std::endl;

    // Check if the selected number of frames is <= than the ones available in the input file
    if (total_num_frames < num_frames_to_read) {
      std::cout << "\n**ERROR**: Select a valid number of frames that is less or equal to the ones available in the input file." << std::endl;
      return 1;
    } else if (num_frames_to_read == -1) {
      num_frames_to_read = total_num_frames;
    }

    //std::unique_ptr<PattgenInfo<> m_pattgen_info;
    //m_pattgen_info = std::make_unique<PattgenInfo<>(500);
    ph.initialize();
    ph.m_pattgen_info->time_tick_offset = time_tick_offset;

    // pattern generation algorithms
    

    // =================================================================
    //                  Process the DUNE WIBEth frames
    // =================================================================
    int wibeth_frame_index = 0;
    first_timestamp = input_file.frame(0)->get_timestamp();
    std::string out_prefix = select_pattern + "_" + std::to_string(time_tick_offset);

    // Loop over the DUNEWIB Ethernet frames in the file
    std::fstream output_file;
    output_file.open(out_prefix+"_wibeth_output.bin", std::ios::app | std::ios::binary);

    ph.m_pattgen_info->num_frames = num_frames_to_read;
    ph.m_pattgen_info->pattern_name = select_pattern;
    ph.m_pattgen_info->first_timestamp = first_timestamp;
    ph.m_pattgen_info->input_ch = input_ch;
    ph.m_pattgen_info->tpg_threshold = tpg_threshold;
    ph.m_pattgen_info->out_prefix = out_prefix;

    while (wibeth_frame_index < num_frames_to_read ){

      // current WIBEth frame
      ph.m_pattgen_info->output_frame = input_file.frame(wibeth_frame_index);
      ph.m_pattgen_info->iframe = wibeth_frame_index;
      execute_tpgpg(*ph.m_pattgen_info, pa);

      ++wibeth_frame_index;

      output_file.write(reinterpret_cast<char*>(ph.m_pattgen_info->output_frame), sizeof(dunedaq::fddetdataformats::WIBEthFrame) );
    }
    output_file.close();

    // =================================================================
    //  Run TPG Hit Finder on pattern-generated DUNE WIBEth frames
    // =================================================================

    if (save_trigprim) {
      execute_tpgpg_val(*ph.m_pattgen_info);
    }



    // =================================================================
    //                       Setup the TPG emulator
    // =================================================================

    /*
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
    emulator->initialize();
    */
    
    // =================================================================
    //                       READ THE HDF5 FILE
    // =================================================================

    /*
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
    */


}

