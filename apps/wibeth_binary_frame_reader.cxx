/**
 * @file wibeth_binary_frame_reader.cxx: reads a WIBEth frame file
 *  and prints the ADC values 
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fddetdataformats/WIBEthFrame.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>


#include <cstring>
#include <immintrin.h>
#include <cstdio> // For printf
#include <array>
#include <chrono>
#include <stdio.h>
#include <stdint.h>


#include "hdf5libs/HDF5RawDataFile.hpp"
#include "hdf5libs/hdf5filelayout/Nljs.hpp"

#include "logging/Logging.hpp"

#include <iostream>
#include <sstream>
#include <string>

using namespace dunedaq::hdf5libs;
using namespace dunedaq::daqdataformats;


void
print_usage()
{
  TLOG() << "Usage: wibeth_binary_frame_reader <input_file_name> ";
}


class FrameFile
{
public:

    FrameFile(const char* filename)
        : m_file(filename, std::ifstream::binary),
          m_buffer(new char[sizeof(dunedaq::fddetdataformats::WIBEthFrame)])
    {
        if(m_file.bad() || m_file.fail() || !m_file.is_open()){
            throw std::runtime_error(std::string("Bad file ")+std::string(filename));
        }
        // Calculate the length of the file
        m_file.seekg(0, m_file.end);
        m_length = m_file.tellg();
        m_file.seekg(0, m_file.beg);
        if(m_length==0){
            throw std::runtime_error("Empty file");
        }
        //if(m_length%sizeof(dunedaq::fddetdataformats::WIBEthFrame)!=0){
        //    throw std::runtime_error("File does not contain an integer number of frames");
        //}
        m_n_frames=m_length/sizeof(dunedaq::fddetdataformats::WIBEthFrame);

	// Reinterpret the frame as WIBEthFrame
	
	m_file.read(m_buffer, m_file.eof());
        m_wibeth_frame = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(m_buffer);
    }

    ~FrameFile()
    {
        m_file.close();
        delete[] m_buffer;
    }

    // Length of the file in bytes
    size_t length() const {return m_length;}
    // Number of frames in the file
    size_t num_frames() const { return m_n_frames; }

    dunedaq::fddetdataformats::WIBEthFrame* get_wibeth_frame() const { return m_wibeth_frame; }

    dunedaq::fddetdataformats::WIBEthFrame* frame(size_t i)
    {
        if(i>=num_frames()) return nullptr;
        // Seek to the right place in the file
        m_file.seekg(i*sizeof(dunedaq::fddetdataformats::WIBEthFrame));
        // Check we didn't go past the end
        if(m_file.bad() || m_file.eof()) return nullptr;
        // Actually read the fragment into the buffer
        m_file.read(m_buffer,sizeof(dunedaq::fddetdataformats::WIBEthFrame));
        if(m_file.bad() || m_file.eof()) return nullptr;
        return reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(m_buffer);
    }
    

    

protected:
    std::ifstream m_file;
    char* m_buffer;
    dunedaq::fddetdataformats::WIBEthFrame* m_wibeth_frame = nullptr;
    size_t m_length;
    size_t m_n_frames;
};




int main(int argc, char** argv)
{
 

  if (argc != 2) {
    print_usage();
    return 1;
  }

  const std::string ifile_name = std::string(argv[1]);

  // Read file
  FrameFile input_file = FrameFile(ifile_name.c_str()); 

  std::cout << "Size of the input file " << input_file.length() << std::endl;
  std::cout << "Number of frames " << input_file.num_frames() << std::endl;
   

  dunedaq::fddetdataformats::WIBEthFrame* output_frame; 
  for (size_t i=0; i<input_file.num_frames(); i++) {
    std::cout << "========== FRAME_NUM " << i <<  std::endl;
    output_frame = input_file.frame(i);
    for (int itime=0; itime<64; ++itime) {
      for (int ch=0; ch<64; ++ch) {
        uint16_t adc_val = output_frame->get_adc(ch, itime);
        std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << ch << " \t\tTimeSample: " << itime <<  std::endl;
      }
   }
  }

  

  
}

