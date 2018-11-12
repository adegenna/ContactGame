#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <string>
#include <fstream>
#include <istream>
#include <random>
#include <cmath>
#include <Eigen/Dense>
#include <boost/program_options.hpp>
#include "options.hpp"
#include "options_parser.hpp"
#include "LagrangianState.h"
#include "ParticlePhysics.h"
#include "Inputfile.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]) {
  printf("*********** DRIVER PROGRAM FOR LAGRANGIAN PARTICLE SOLVER ***********\n\n");

  // Parse input file options
  Options options;
  if (!parseOptions(argc,argv,options)){
   return 0;
  }
  cout << options << endl;

  // Load state
  MatrixXd input = load_csv<MatrixXd>(options.projDir + options.loadDir + options.inputfile);
  
  // Pass parsed program options to simulation
  LagrangianState simulation(input);
  ParticlePhysics physics(options, simulation);
  
  // Solve
  physics.simulate();

  // Output
  const std::string filename = options.projDir + options.outDir + "XY_final.csv";
  simulation.writeXY(filename);
  
  return 0;
}
