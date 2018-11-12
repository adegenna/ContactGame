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

using namespace Eigen;
namespace po = boost::program_options;

// ***************************************************
// CLASS FOR PARTICLEPHYSICS
// ***************************************************

ParticlePhysics::ParticlePhysics() {

}

ParticlePhysics::ParticlePhysics(Options& o, LagrangianState& simulation) {
  
  // Parameters
  projdir_    = o.projDir;
  outdir_     = o.outDir;
  loaddir_    = o.loadDir;
  dt_         = o.dt;
  integrator_ = o.integrator;
  tsteps_     = o.tsteps;
  tsave_      = o.tsave;
  // Simulation
  simulation_ = &simulation;
  
}

ParticlePhysics::~ParticlePhysics() {
  
}

void ParticlePhysics::eulerDXY(MatrixXd& DXY) {
  MatrixXd UV = simulation_->getUV();
  DXY.resize(samples_,2);
  DXY = UV*dt_;

}

void ParticlePhysics::modelContactForces(int i, int j, VectorXd& dij) {
  // Potential function for collision force calculation
  double distance_ij = dij.norm();
  double eps   = 0.001;
  VectorXd R   = simulation_->getR();
  double delta = std::min( std::abs((R(i)+R(j))-distance_ij) , eps*(R(i)+R(j)) );
  double F     = pow(10.0,5.0)*pow(delta,0.85);
  VectorXd eij = dij/distance_ij;
  // Update forces
  forces_(i,0) -= F*eij(0); forces_(i,1) -= F*eij(1);
  forces_(j,0) += F*eij(0); forces_(j,1) += F*eij(1);

}

void ParticlePhysics::particleContact() {
  // Calculate all particle-particle contacts (N^2 brute force)
  MatrixXd XY = simulation_->getXY();
  MatrixXd UV = simulation_->getUV();
  VectorXd R  = simulation_->getR();
  VectorXd dij(2);
  double distance_ij, ui_tangent, uj_tangent;
  # pragma omp parallel for
  for (int i=0; i<samples_; i++) {
    for (int j=i+1; j<samples_; j++) {
      dij         = XY.row(j)-XY.row(i);
      distance_ij = dij.norm();
      if (distance_ij <= R(i)+R(j)) {
	ui_tangent = UV.row(i).dot(dij)/distance_ij;
	uj_tangent = UV.row(j).dot(dij)/distance_ij;
	if ((ui_tangent - uj_tangent) > 0) {
	  modelContactForces(i,j,dij);
	}
      }
    }
  }

}

void ParticlePhysics::calculateParticleMasses() {
  // Assume all particles have the same density
  double rho = 1.0;
  mass_.resize(samples_);
  VectorXd R = simulation_->getR();
  mass_ = rho*0.5*M_PI*R.array().pow(2);
  
}

void ParticlePhysics::updateParticleVelocities() {  
  double diff   = 0.7;
  MatrixXd DUV(samples_,2);
  DUV.col(0) = diff*dt_*(forces_.col(0).array()/mass_.array());
  DUV.col(1) = diff*dt_*(forces_.col(1).array()/mass_.array());
  simulation_->incrementUV(DUV);
  
}

void ParticlePhysics::simulate() {
  MatrixXd DXY;
  samples_ = simulation_->getSamples();
  calculateParticleMasses();
  for (int i=0; i<tsteps_; i++) {
    forces_ = MatrixXd::Zero(samples_,2);
    // Calculate contacts and forces
    particleContact();
    // Convert forces to acceleration
    updateParticleVelocities();
    // Advance in time
    eulerDXY(DXY);
    // Update state
    simulation_->updateXY(DXY);
    if ( (i+1)%tsave_ == 0) {
      const std::string filename = projdir_ + outdir_ + "XY_" + std::to_string(i+1) + ".csv";
      simulation_->writeXY(filename);
    }
  }

}
