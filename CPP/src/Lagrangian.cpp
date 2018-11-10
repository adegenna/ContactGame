#include "Lagrangian.h"
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

Lagrangian::Lagrangian() {

}

Lagrangian::Lagrangian(Options& o) {

  std::cout << "****************** LAGRANGIAN INITIALIZATION ****************" << endl << endl;
  // Parameters
  inputfile_ = o.inputfile;
  
}

Lagrangian::~Lagrangian() {

}

Lagrangian::setupInitialConditions() {
  input_   = load_csv<MatrixXd>(inputfile_);
  samples_ = input_.rows();
  
}

