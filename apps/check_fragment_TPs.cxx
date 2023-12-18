/**
 * @file check_fragment_TPs.cxx Read TP fragments from file 
 * and check that they have start times within the request window
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "CLI/CLI.hpp"

#include "trgdataformats/TriggerPrimitive.hpp"
#include "hdf5libs/HDF5RawDataFile.hpp"

#include <daqdataformats/Fragment.hpp>
#include <daqdataformats/FragmentHeader.hpp>
#include <daqdataformats/TriggerRecordHeader.hpp>
#include <daqdataformats/Types.hpp>
#include <iostream>
#include <map>
#include <memory>
#include <regex>

// A little struct to hold a TriggerRecordHeader along with the
// corresponding TP fragments
struct MyTriggerRecord
{
  std::unique_ptr<dunedaq::daqdataformats::TriggerRecordHeader> header;
  std::vector<std::unique_ptr<dunedaq::daqdataformats::Fragment>> fragments;
};

int
main(int argc, char** argv)
{
  CLI::App app{ "App description" };

  std::string filename;
  app.add_option("-f,--file", filename, "Input HDF5 file");

  CLI11_PARSE(app, argc, argv);

  // Map from trigger number to struct of header/fragments
  std::map<int, MyTriggerRecord> trigger_records;

  std::regex header_regex("^//TriggerRecord([0-9]+)/TriggerRecordHeader$", std::regex::extended);
  std::regex trigger_fragment_regex("^//TriggerRecord([0-9]+)/Trigger/Element[0-9]+$",
                                    std::regex::extended);

  dunedaq::hdf5libs::HDF5RawDataFile decoder(filename);

  auto trigger_record_ids = decoder.get_all_trigger_record_ids();
  
  // Populate the map with the TRHs and DS fragments
  for (auto trigger_record_id : trigger_record_ids){
    auto trigger_number = trigger_record_id.first;

    trigger_records[trigger_number].header = decoder.get_trh_ptr(trigger_number);
    auto seq_number = trigger_record_id.second;
    
    for(size_t ic=0; ic < trigger_records[trigger_number].header->get_num_requested_components(); ++ic){

      auto const& comp_sourceid = trigger_records[trigger_number].header->at(ic).component;

      if(comp_sourceid.subsystem != dunedaq::daqdataformats::SourceID::Subsystem::kTrigger)
	continue;

      trigger_records[trigger_number].fragments.push_back(decoder.get_frag_ptr(trigger_number,seq_number, comp_sourceid));
    }

  }
  
  int n_failures = 0;

  using dunedaq::trgdataformats::TriggerPrimitive;
  for (auto const& [trigger_number, record] : trigger_records) {
    std::cout << "Trigger number " << trigger_number << " with TRH pointer " << record.header.get() << " and "
              << record.fragments.size() << " fragments" << std::endl;

    // Find the window start and end requested for this trigger
    // record. In principle the request windows for each data
    // selection component could be different, but we'll assume
    // they're the same for now, because matching up the component
    // request to the fragment is too difficult
    size_t n_requests = record.header->get_num_requested_components();
    dunedaq::daqdataformats::timestamp_t window_begin = 0, window_end = 0;
    for (size_t i = 0; i < n_requests; ++i) {
      auto request = record.header->at(i);
      if (request.component.subsystem == dunedaq::daqdataformats::SourceID::Subsystem::kTrigger) {
        window_begin = request.window_begin;
        window_end = request.window_end;
      }
    }
    if (window_begin == 0 || window_end == 0) {
      // We didn't find a component request for a dataselection item, so skip
      continue;
    }
    // Check that each primitive in each fragment falls within the request window
    for (auto const& frag : record.fragments) {
      const size_t n_prim =
        (frag->get_size() - sizeof(dunedaq::daqdataformats::FragmentHeader)) / sizeof(TriggerPrimitive);
      std::cout << "  Fragment has " << n_prim << " primitives" << std::endl;
      const TriggerPrimitive* prim = reinterpret_cast<TriggerPrimitive*>(frag->get_data());
      for (size_t i = 0; i < n_prim; ++i) {
        if (prim->time_start < window_begin || prim->time_start > window_end) {
          std::cout << "Primitive with time_start " << prim->time_start << " is outside request window of ("
                    << window_begin << ", " << window_end << ")" << std::endl;
          ++n_failures;
        }
        ++prim;
      }
    }
  }
  if (n_failures > 0) {
    std::cout << "Found " << n_failures << " TPs outside window in " << trigger_record_ids.size() << " trigger records" << std::endl;
  } else {
    std::cout << "Test passed" << std::endl;
  }
  return (n_failures != 0);
}