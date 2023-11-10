/**
 * @file WIB2AddFakeHits.cxx: quick and dirty solution to modify an input binary file and add 
 * fake hits at a predefined frame and channel numbers 
 * Usage: ./WIB2AddFakeHits /nfs/home/aabedabu/dunedaq-v3.1.1dunedaq-v3.1.1_pre_source_id/dev/test_files/output_0_3.out 
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fddetdataformats/WIB2Frame.hpp"

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
  TLOG() << "Usage: WIB2AddFakeHits <input_file_name> <input_channel>";
}


class FrameFile
{
public:

    FrameFile(const char* filename)
        : m_file(filename, std::ifstream::binary),
          m_buffer(new char[sizeof(dunedaq::fddetdataformats::WIB2Frame)])
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
        //if(m_length%sizeof(dunedaq::fddetdataformats::WIB2Frame)!=0){
        //    throw std::runtime_error("File does not contain an integer number of frames");
        //}
        m_n_frames=50;//m_length/sizeof(dunedaq::fddetdataformats::WIB2Frame);

	// Reinterpret the frame as WIB2Frame
	
	m_file.read(m_buffer, m_file.eof());
        m_wib2_frame = reinterpret_cast<dunedaq::fddetdataformats::WIB2Frame*>(m_buffer);
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

    dunedaq::fddetdataformats::WIB2Frame* get_wib2_frame() const { return m_wib2_frame; }

    dunedaq::fddetdataformats::WIB2Frame* frame(size_t i)
    {
        if(i>=num_frames()) return nullptr;
        // Seek to the right place in the file
        m_file.seekg(i*sizeof(dunedaq::fddetdataformats::WIB2Frame));
        // Check we didn't go past the end
        if(m_file.bad() || m_file.eof()) return nullptr;
        // Actually read the fragment into the buffer
        m_file.read(m_buffer,sizeof(dunedaq::fddetdataformats::WIB2Frame));
        if(m_file.bad() || m_file.eof()) return nullptr;
        return reinterpret_cast<dunedaq::fddetdataformats::WIB2Frame*>(m_buffer);
    }
    

    

protected:
    std::ifstream m_file;
    char* m_buffer;
    dunedaq::fddetdataformats::WIB2Frame* m_wib2_frame = nullptr;
    size_t m_length;
    size_t m_n_frames;
};




int main(int argc, char** argv)
{
 

  if (argc != 3) {
    print_usage();
    return 1;
  }

  const std::string ifile_name = std::string(argv[1]);

  const size_t input_ch = atoi(argv[2]);

  // Read file
  FrameFile input_file = FrameFile(ifile_name.c_str()); 

  std::cout << "Size of the input file " << input_file.length() << std::endl;
  std::cout << "Number of frames " << input_file.num_frames() << std::endl;
   
    
  
 
  dunedaq::fddetdataformats::WIB2Frame* output_frame; 
  for (size_t i=0; i<input_file.num_frames(); i++) {
     output_frame = input_file.frame(i);
     auto val_to_check = output_frame->get_adc(input_ch);
     std::cout << "Output ADC value: " << val_to_check << " for frame " << i << " and channel number " << input_ch << std::endl;

  } 


}

