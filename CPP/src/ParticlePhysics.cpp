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

using namespace std;
using namespace Eigen;
namespace po = boost::program_options;

ParticlePhysics::ParticlePhysics(Options& o, LagrangianState& simulation)
  : options_(o), simulation_(&simulation)
{
  samples_ = simulation_->getSamples();
}

ParticlePhysics::~ParticlePhysics() {
  
}

void ParticlePhysics::modelContactForces(int i, int j, VectorXd& dij) {
  // Potential function for collision force calculation
  double distance_ij = dij.norm();
  double eps   = 0.001;
  const VectorXd& R = simulation_->getR();
  double delta = std::min( std::abs((R(i)+R(j))-distance_ij) , eps*(R(i)+R(j)) );
  double F     = pow(10.0,5.0)*pow(delta,0.85);
  VectorXd eij = dij/distance_ij;
  // Update forces
  forces_(i,0) -= F*eij(0); forces_(i,1) -= F*eij(1);
  forces_(j,0) += F*eij(0); forces_(j,1) += F*eij(1);

}

void ParticlePhysics::particleContact() {
  // Calculate all particle-particle contacts (N^2 brute force)
  const MatrixXd& XY = simulation_->getXY();
  const MatrixXd& UV = simulation_->getUV();
  const VectorXd& R  = simulation_->getR();
  VectorXd dij(2);
  double distance_ij, ui_tangent, uj_tangent;
  //# pragma omp parallel for
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
  const VectorXd& R = simulation_->getR();
  mass_ = rho*0.5*M_PI*R.array().pow(2);
  
}

void ParticlePhysics::updateParticleVelocities() {  
  double diff   = 3.77;//0.7;
  MatrixXd DUV(samples_,2);
  DUV.col(0) = diff*options_.dt*(forces_.col(0).array()/mass_.array());
  DUV.col(1) = diff*options_.dt*(forces_.col(1).array()/mass_.array());
  simulation_->incrementUV(DUV);
  
}

double ParticlePhysics::calculateTotalEnergy() {
  double energy = 0.0;
  MatrixXd UV = simulation_->getUV();
  for (int i=0; i<samples_; i++)
    energy += 0.5*mass_(i)*(pow(UV(i,0),2) + pow(UV(i,1),2));
  return energy;
      
}

VectorXd ParticlePhysics::calculateTotalMomentum() {
  VectorXd momentum(2);
  momentum(0) = 0.0; momentum(1) = 0.0;
  MatrixXd UV = simulation_->getUV();
  for (int i=0; i<samples_; i++) {
    momentum(0) += mass_(i)*UV(i,0);
    momentum(1) += mass_(i)*UV(i,1);
  } 
  return momentum;
      
}

void ParticlePhysics::writeEnergyAndMomentum(double energy, VectorXd& momentum, const std::string& filename) const {
  ofstream xyout(filename);
  xyout << energy << ", " << momentum(0) << ", " << momentum(1) << std::endl;
  xyout.close();

}

void ParticlePhysics::zeroForces() {
  forces_  = MatrixXd::Zero(samples_,2);
}
