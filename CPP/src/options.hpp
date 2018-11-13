#pragma once
#include <ostream>
#include <string>

struct Options {
  std::string inputfile;
  std::string outputfile;
  double dt;
  std::string integrator;
  int tsteps;
  int tsave;
  
};

inline std::ostream& operator <<(std::ostream& os, const Options& o){
 os << "inputfile  = " << o.inputfile << std::endl
    << "outputfile = " << o.outputfile << std::endl
    << "dt         = " << o.dt        << std::endl
    << "integrator = " << o.integrator << std::endl
    << "tsteps     = " << o.tsteps     << std::endl
    << "tsave      = " << o.tsave      << std::endl;

   return os;
}
