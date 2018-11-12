#include "LagrangianSimulation.h"
#include "Inputfile.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <random>
#include <cmath>
#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <omp.h>
#include <sys/time.h>

namespace po = boost::program_options;

// ***************************************************
// CLASS FOR LAGRANGIAN COEFFICIENT INFERENCE
// ***************************************************

LagrangianSimulation::LagrangianSimulation() {

}

LagrangianSimulation::LagrangianSimulation(Options& o) {

  // Parameters
  inputfile_  = o.inputfile;
  projdir_    = o.projDir;
  outdir_     = o.outDir;
  loaddir_    = o.loadDir;
  
}

LagrangianSimulation::~LagrangianSimulation() {

}

void LagrangianSimulation::setupInitialConditions() {
  input_   = load_csv<MatrixXd>(projdir_ + loaddir_ + inputfile_);
  samples_ = input_.rows();
  XY_      = input_.block(0,0,samples_,2);
  UV_      = input_.block(0,2,samples_,2);
  R_       = input_.block(0,4,samples_,1);
  
}

void LagrangianSimulation::updateXY(MatrixXd& DXY) {
  XY_ += DXY;

}

void LagrangianSimulation::writeXY(std::string append) {
  const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
  ofstream xyout(projdir_ + outdir_ + "XY_" + append + ".csv");
  xyout << XY_.format(CSVFormat);
  xyout.close();
  
}