#pragma once
#include <ostream>
#include <string>

struct Options {
  std::string projDir;
  std::string outDir;
  std::string loadDir;
  std::string inputfile;
  double dt;
  std::string integrator;
  int tsteps;
  int tsave;
  
};

inline std::ostream& operator <<(std::ostream& os, const Options& o){
 os << "projDir    = " << o.projDir  << std::endl
    << "outDir     = " << o.outDir   << std::endl
    << "loadDir    = " << o.loadDir  << std::endl
    << "inputfile  = " << o.inputfile << std::endl
    << "dt         = " << o.dt        << std::endl
    << "integrator = " << o.integrator << std::endl
    << "tsteps     = " << o.tsteps     << std::endl
    << "tsave      = " << o.tsave      << std::endl;

   return os;
}
