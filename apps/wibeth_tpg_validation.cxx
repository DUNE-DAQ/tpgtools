/**
 * @file wibeth_tpg_validation.cxx 
 * Application for testing different TPG algorithms (naive/avx) 
 * on known pattern frame files
 * 
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */


#include "tpgtools/TPGEmulator.hpp"

// system
#include "CLI/CLI.hpp"
#include <iostream>

#include <fmt/core.h>
#include <fmt/format.h>


// =================================================================
//                       MAIN
// =================================================================
int 
main(int argc, char** argv)
{

    CLI::App app{ "Test TPG algorithms using binary frame files." };
    // Set default input frame file
    std::string file_path_input = "";
    app.add_option("-f,--file-path-input", file_path_input, "Path to the input file");

    std::string select_algorithm = "SimpleThreshold";
    app.add_option("-a,--algorithm", select_algorithm, "TPG Algorithm (SimpleThreshold / AbsRS). Default: SimpleThreshold");
  
    std::string select_implementation = "AVX";
    app.add_option("-i,--implementation", select_implementation, "TPG implementation (AVX / NAIVE). Default: AVX");

    int tpg_threshold = 500;
    app.add_option("-t,--tpg-threshold", tpg_threshold, "Value of the TPG threshold. Default value is 500.");

    int core_number = 0;
    app.add_option("-c,--core", core_number, "Set core number of the executing TPG thread. Default value is 0.");

    bool save_adc_data = false;
    app.add_flag("--save-adc-data", save_adc_data, "Save ADC data (first frame only)");

    bool save_trigprim = false;
    app.add_flag("--save-trigprim", save_trigprim, "Save trigger primitive data");


    CLI11_PARSE(app, argc, argv);



    
    // =================================================================
    //                       READ THE FRAMES.BIN FILE
    // =================================================================

    const int wibeth_frame_size = dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter::fixed_payload_size;
    std::unique_ptr<dunedaq::readoutlibs::FileSourceBuffer> m_source_buffer;
    // Read only 10 MB of data
    m_source_buffer = std::make_unique<dunedaq::readoutlibs::FileSourceBuffer>(10485100, wibeth_frame_size); 

    m_source_buffer->read(file_path_input);
    auto& source = m_source_buffer->get(); 
    const int total_num_frames = m_source_buffer->num_elements(); // file_ size/chunk_size = 180 

    fmt::print("Number of DUNE WIBEth frames in the input file: {} \n", total_num_frames);

    
    // =================================================================
    //                       Setup the SWTPG
    // =================================================================

    // Create instance of the TPG emulator implementation    
    std::unique_ptr<tpg_emulator_base> emulator;
    if (select_implementation == "AVX") {
      emulator = std::make_unique<tpg_emulator_avx>(save_adc_data, save_trigprim, false, select_algorithm, "");
    } else if (select_implementation == "NAIVE") {
      emulator = std::make_unique<tpg_emulator_naive>(save_adc_data, save_trigprim, false, select_algorithm, "");
    } else {
      throw tpgtools::InvalidImplementation(ERS_HERE, select_implementation);  
    }
    
    emulator->set_tpg_threshold(tpg_threshold);
    emulator->set_CPU_affinity(core_number);
    emulator->initialize();


    // =================================================================
    //                  Process the DUNE WIBEth frames
    // =================================================================
    int wibeth_frame_index = 0; 

    // Loop over the DUNEWIB Ethernet frames in the file
    while (wibeth_frame_index < total_num_frames ){      

      // current WIBEth frame
      auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter*>(source.data() + wibeth_frame_index*wibeth_frame_size);

      emulator->execute_tpg(fp);

      ++wibeth_frame_index;


    }
    fmt::print("\n\n=============================== \n");
    fmt::print("Found in total {} hits. \n", emulator->get_total_hits());
    


}


