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

    // Set default input frame file
    std::string file_path_input = "";
    app.add_option("-f,--file-path-input", file_path_input, "Path to the input file");

    // Set the output path
    std::string path_output = ".";
    app.add_option("-o,--path-output", path_output, "Path to the output directory. Default: . (i.e. pwd)");     

    int num_frames_to_read = -1;
    app.add_option("-n,--num-frames-to-read", num_frames_to_read, "Number of WIBEth frames to read. Default: select all frames.");

    size_t input_ch = 0;
    app.add_option("-i,--input_channel", input_ch, "Input channel number for adding fake hit. Default: 0. (max:63)");

    int tpg_threshold = 500;
    app.add_option("-t,--tpg-threshold", tpg_threshold, "Value of the TPG threshold. Default value is 500.");

    bool save_adc_data = false;
    app.add_flag("--save-adc-data", save_adc_data, "Save ADC data (first frame only)");

    bool save_trigprim = false;
    app.add_flag("--save-trigprim", save_trigprim, "Save trigger primitive data");

    int clock_tick_offset = 1;
    app.add_option("-c,--clock-tick-offset", clock_tick_offset, "Time tick of pattern start. Default: 1 (max:63).");

    std::string out_suffix = "";
    app.add_option("-s ,--out_suffix", out_suffix, "Append string (suffix) to output hit file name");

    std::string select_pattern = "patt_golden";
    app.add_option("-p,--select-pattern", select_pattern, "Test pattern name (patt_golden, patt_pulse, patt_square, patt_square_left, patt_square_right). Default: patt_golden.");

    bool overwrite_wibeth_header = false;
    app.add_flag("-w,--overwrite-wibeth-header", overwrite_wibeth_header, "Overwrite crate, slot, stream IDs (needed for offline channel map). Default: false.");

    bool verbose = false;
    app.add_flag("-v,--verbose", verbose, "Printout additional information while the application is running. Default: false.");

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

    TLOG() << "Input file name: " << file_path_input;
    TLOG() << "Size of the input file: " << input_file.length();
    TLOG() << "Input channel number: " << input_ch;
    TLOG() << "Number of DUNE WIBEth frames in the input file: " << total_num_frames << std::endl;

    // Check if the selected number of frames is <= than the ones available in the input file
    if (total_num_frames < num_frames_to_read) {
      TLOG() << "\n**ERROR**: Select a valid number of frames that is less or equal to the ones available in the input file.";
      return 1;
    } else if (num_frames_to_read == -1) {
      num_frames_to_read = total_num_frames;
    }

    //std::unique_ptr<PattgenInfo<> m_pattgen_info;
    //m_pattgen_info = std::make_unique<PattgenInfo<>(500);
    ph.initialize();
    ph.m_pattgen_info->clock_tick_offset = clock_tick_offset;

    // pattern generation algorithms
    

    // =================================================================
    //                  Process the DUNE WIBEth frames
    // =================================================================
    int wibeth_frame_index = 0;
    first_timestamp = input_file.frame(0)->get_timestamp();
    std::string out_prefix = select_pattern + "_chan_" + std::to_string(input_ch) + "_tick_" + std::to_string(clock_tick_offset);

    // Loop over the DUNEWIB Ethernet frames in the file
    std::fstream output_file;
    output_file.open(path_output+"/"+out_prefix+"_wibeth_output.bin", std::ios::app | std::ios::binary);

    ph.m_pattgen_info->num_frames = num_frames_to_read;
    ph.m_pattgen_info->pattern_name = select_pattern;
    ph.m_pattgen_info->first_timestamp = first_timestamp;
    ph.m_pattgen_info->input_ch = input_ch;
    ph.m_pattgen_info->tpg_threshold = tpg_threshold;
    ph.m_pattgen_info->out_prefix = out_prefix;
    ph.m_pattgen_info->path_output = path_output;
    ph.m_pattgen_info->overwrite_wibeth_header = overwrite_wibeth_header;
    ph.m_pattgen_info->verbose = verbose;

    while (wibeth_frame_index < num_frames_to_read ){

      // =================================================================
      //  Generate test pattern and store DUNE WIBEth frames
      // =================================================================

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
      execute_tpgpg_validation(*ph.m_pattgen_info);
    }

}

