/**
 * @file tpg_apps_issues.hpp
 * 
 * ERS issues for TPG applications
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef TPGTOOLS_INCLUDE_TPGTOOLS_TPGISSUES_HPP_
#define TPGTOOLS_INCLUDE_TPGTOOLS_TPGISSUES_HPP_


#include <ers/Issue.hpp>
#include <string>


ERS_DECLARE_ISSUE(tpgtools,
                  TPGAlgorithmInexistent,
                  "The selected algorithm does not exist: " << algorithm_selection << " . Check your command line options and select either SimpleThreshold or AbsRS",
                  ((std::string)algorithm_selection))

ERS_DECLARE_ISSUE(tpgtools,
                  FileInexistent,
                  "The selected input file does not exist. Input file: " << input_file_path << "  Check the path of the input file",
                  ((std::string)input_file_path))

ERS_DECLARE_ISSUE(tpgtools,
                  InvalidImplementation,
                  "The selected TPG algorithm implementation does not exist: " << input_file_path << "  Check your command line options and select either NAIVE or AVX",
                  ((std::string)input_file_path))




#endif // TPGTOOLS_INCLUDE_TPGTOOLS_TPGISSUES_HPP_
