#include "LagrangianState.h"
#include "Inputfile.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <random>
#include <cmath>
#include <iostream>
#include <iterator>
#include <omp.h>
#include <sys/time.h>

using namespace std;
using namespace Eigen;

LagrangianState::LagrangianState(const MatrixXd& initialstate) {

  // Parameters
  samples_    = initialstate.rows();
  XY_         = initialstate.block(0,0,samples_,2);
  UV_         = initialstate.block(0,2,samples_,2);
  R_          = initialstate.block(0,4,samples_,1);
}

LagrangianState::~LagrangianState() {

}

void LagrangianState::updateXY(const MatrixXd& DXY) {
  XY_ += DXY;

}

void LagrangianState::writeXY(const std::string& filename) const {
  const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
  ofstream xyout(filename);
  xyout << XY_.format(CSVFormat);
  xyout.close();
  
}
