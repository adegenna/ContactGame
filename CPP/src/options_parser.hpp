#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include "options.hpp"

namespace po = boost::program_options;

bool parseOptions(int argc, char *argv[], Options& options){

  std::string input_file;

  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("projDir"          , po::value<std::string>(&options.projDir)   , "Set project directory")
    ("outDir"           , po::value<std::string>(&options.outDir)    , "Set output directory")    
    ("loadDir"          , po::value<std::string>(&options.loadDir)   , "Set load directory")    
    ("inputfile"        , po::value<std::string>(&options.inputfile) , "Set input filename")
    ("dt"               , po::value<double>(&options.dt) , "Set timestep")
    ("integrator"       , po::value<std::string>(&options.integrator) , "Set time integrator type")
    ("tsteps"           , po::value<int>(&options.tsteps)     , "Set number of time steps")
    ("tsave"            , po::value<int>(&options.tsave)      , "Set checkpoint for saving")

    ("config_file"      , po::value<std::string>(&input_file), "Configuration input file.")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
      std::cout << desc << "\n";
      return false;
  }

  std::ifstream ifs(input_file.c_str());
  po::store(po::parse_config_file(ifs, desc), vm);
  po::notify(vm);

  return true;
}

