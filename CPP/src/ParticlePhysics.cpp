#include "ParticlePhysics.h"
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
// CLASS FOR PARTICLEPHYSICS
// ***************************************************

ParticlePhysics::ParticlePhysics() {

}

ParticlePhysics::ParticlePhysics(Options& o, LagrangianSimulation& simulation) {
  
  // Parameters
  dt_         = o.dt;
  integrator_ = o.integrator;
  tsteps_     = o.tsteps;
  // Simulation
  simulation_ = simulation;
  
}

ParticlePhysics::~ParticlePhysics() {

}

void ParticlePhysics::eulerDXY(MatrixXd& DXY) {
  int samples = simulation_.getSamples();
  MatrixXd UV = simulation_.getUV();
  DXY.resize(samples,2);
  DXY = UV*dt_;

}

void ParticlePhysics::simulate() {
  MatrixXd DXY;
  for (int i=0; i<tsteps_; i++) {
    eulerDXY(DXY);
    simulation_.updateXY(DXY);
  }

}
