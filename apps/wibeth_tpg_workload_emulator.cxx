/**
 * @file wibeth_tpg_workload_emulator.cxx 
 * Emulate workload for TPG algorithms for performance measurements
 * 
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "readoutlibs/utils/RateLimiter.hpp"
#include "readoutlibs/utils/FileSourceBuffer.hpp"
#include "readoutlibs/utils/BufferedFileWriter.hpp"

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

    CLI::App app{ "Generate workload for TPG algorithms using binary frame files." };
    // Set default input frame file
    std::string file_path_input = "";
    app.add_option("-f,--file-path-input", file_path_input, "Path to the input file");

    // Set the output path
    std::string path_output = ".";
    app.add_option("-o,--path-output", path_output, "Path to the output directory. Default: .");     

    std::string select_algorithm = "SimpleThreshold";
    app.add_option("-a,--algorithm", select_algorithm, "TPG Algorithm (SimpleThreshold / AbsRS). Default: SimpleThreshold");
  
    std::string select_implementation = "AVX";
    app.add_option("-i,--implementation", select_implementation, "TPG implementation (AVX / NAIVE). Default: AVX");

    std::string select_channel_map = "None";
    app.add_option("-m,--channel-map", select_channel_map, "Select a valid channel map: None, VDColdboxChannelMap, ProtoDUNESP1ChannelMap, PD2HDChannelMap, HDColdboxChannelMap, FiftyLChannelMap");

    int duration_test = 120; 
    app.add_option("-d,--duration-test", duration_test, "Duration (in seconds) to run the test. Default value is 120.");

    int tpg_threshold = 500;
    app.add_option("-t,--tpg-threshold", tpg_threshold, "Value of the TPG threshold. Default value is 500.");

    int core_number = 0;
    app.add_option("-c,--core", core_number, "Set core number of the executing TPG thread. Default value is 0.");

    bool save_adc_data = false;
    app.add_flag("--save-adc-data", save_adc_data, "Save ADC data (first frame only)");

    bool save_trigprim = false;
    app.add_flag("--save-trigprim", save_trigprim, "Save trigger primitive data");

    // Additional options for validation 
    bool repeat_timer = true;
    app.add_option("-r, --repeat-timer", repeat_timer, "Repeat frame processing until certain time elapsed (true/false). Default: true.");

    std::string out_suffix = "";
    app.add_option("-s ,--out-suffix", out_suffix, "Append string to output hit file name (e.g. __1). Default: empty string).");

    int num_frames_to_read = -1;
    app.add_option("-n,--num-frames-to-read", num_frames_to_read, "Number of frames to read. Default: -1 (select all frames).");


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
    int total_num_frames = m_source_buffer->num_elements(); // file_ size/chunk_size = 180 

    fmt::print("Number of DUNE WIBEth frames in the input file: {} \n", total_num_frames);
    fmt::print("Number of DUNE WIBEth frames to read: {} \n", num_frames_to_read);
    total_num_frames = num_frames_to_read; 

    // =================================================================
    //                       Setup the SWTPG
    // =================================================================

    // Create instance of the TPG emulator implementation    
    std::unique_ptr<tpg_emulator_base> emulator;
    if (select_implementation == "AVX") {
      emulator = std::make_unique<tpg_emulator_avx>(save_adc_data, save_trigprim, false, select_algorithm, select_channel_map);
    } else if (select_implementation == "NAIVE") {
      emulator = std::make_unique<tpg_emulator_naive>(save_adc_data, save_trigprim, false, select_algorithm, select_channel_map);
    } else {
      throw tpgtools::InvalidImplementation(ERS_HERE, select_implementation);  
    }
    
    emulator->set_tpg_threshold(tpg_threshold);
    emulator->set_CPU_affinity(core_number);
    emulator->initialize();
    emulator->set_out_suffix(out_suffix);
    emulator->set_path_output(path_output);


    // Setup the rate limiter for WIBEth frames
    auto limiter = dunedaq::readoutlibs::RateLimiter(31);
    limiter.init();



    // =================================================================
    //                  Process the DUNE WIBEth frames
    // =================================================================
    int wibeth_frame_index = 0; 
    uint64_t frame_repeat_index = 0;
    auto start_test = std::chrono::high_resolution_clock::now();  

    // Loop over the DUNEWIB Ethernet frames in the file
    while (wibeth_frame_index < total_num_frames ){

      // current WIBEth frame
      auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter*>(source.data() + wibeth_frame_index*wibeth_frame_size);
      //auto fr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(
      //          static_cast<char*>((char*)source.data()) +  wibeth_frame_index*sizeof(dunedaq::fddetdataformats::WIBEthFrame));


      //auto fr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(&source.data()[wibeth_frame_index]);
      //auto fr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(fp);
      auto fr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>((uint8_t*)fp);

      emulator->register_channel_map(fr);
      //auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter*>(fr);

      emulator->execute_tpg(fp);

      ++wibeth_frame_index;

      if (!repeat_timer) {
        if (wibeth_frame_index == total_num_frames) {
          continue;
        }
      }

      // If end of the file is reached, restart the index counter
      if (wibeth_frame_index == total_num_frames) {
        wibeth_frame_index = 0;
	      frame_repeat_index++;
      }

      if (frame_repeat_index % 100  == 0) {
        // Calculate elapsed time in seconds  
        auto now = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(now - start_test).count();  
        fmt::print("Elapsed time [s]: {} \n", elapsed_seconds);
	
	      frame_repeat_index = 0;        

        // stop the testing after a time a condition
        if (elapsed_seconds > duration_test) {
          wibeth_frame_index = total_num_frames;
        }
      }

	

      limiter.limit();

    }    
    fmt::print("\n\n=============================== \n");
    fmt::print("Found in total {}  hits. \n", emulator->get_total_hits());

}


