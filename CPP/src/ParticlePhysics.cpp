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
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <iostream>
#include <iterator>
#include <omp.h>
#include <sys/time.h>

using namespace std;
using namespace Eigen;
namespace po = boost::program_options;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

BOOST_GEOMETRY_REGISTER_POINT_2D(Eigen::Vector2d, double, boost::geometry::cs::cartesian, Eigen::Vector2d::x(), Eigen::Vector2d::y())

ParticlePhysics::ParticlePhysics(Options& o, LagrangianState& simulation)
  : options_(o), simulation_(&simulation)
{
  samples_ = simulation_->getSamples();
  calculateParticleMasses();
}

ParticlePhysics::~ParticlePhysics() {
  
}

void ParticlePhysics::modelContactForces(int i, int j, Vector2d& dij) {
  // Potential function for collision force calculation
  double distance_ij = dij.norm();
  double eps         = 0.001;
  const VectorXd& R  = simulation_->getR();
  double delta = std::min( std::abs((R(i)+R(j))-distance_ij) , eps*(R(i)+R(j)) );
  double F     = pow(10.0,5.0)*pow(delta,0.85);
  VectorXd eij = dij/distance_ij;
  // Update forces
  forces_(i,0) -= F*eij(0); forces_(i,1) -= F*eij(1);
  forces_(j,0) += F*eij(0); forces_(j,1) += F*eij(1);

}

void ParticlePhysics::particleContactRtree() {
  // Calculate all particle-particle contacts (N log N)
  const MatrixXd& XY = simulation_->getXY();
  const MatrixXd& UV = simulation_->getUV();
  const VectorXd& R  = simulation_->getR();
  // Rtree construction
  typedef std::pair<Vector2d, unsigned> value;
  bgi::rtree< value, bgi::quadratic<16> > rtree;
  for (int i=0; i<samples_; i++) {
    Vector2d xy; xy[0] = XY(i,0); xy[1] = XY(i,1);
    rtree.insert(std::make_pair(xy,i));
  }
  int maxN = std::min(samples_,6);
  // knn search
  double epsilon = 1e-8;
  initializeParticleInteractionTracker();
  for (int i=0; i<samples_; i++) {
    std::vector<value> result_n;
    Vector2d xy; xy[0] = XY(i,0); xy[1] = XY(i,1);
    rtree.query(bgi::nearest(xy,maxN), std::back_inserter(result_n));
    for (int j=0; j<maxN; j++) {
      const value& v      = result_n[j];
      int idx             = v.second;
      Vector2d dij; dij[0] = v.first[0]-XY(i,0); dij[1] = v.first[1]-XY(i,1);
      double distance_ij  = dij.norm();
      if ((interactions_(i,idx) == 0) || (interactions_(idx,i) == 0)) {
	double ui_tangent = UV.row(i).dot(dij)/distance_ij;
	double uj_tangent = UV.row(idx).dot(dij)/distance_ij;
	if ((ui_tangent - uj_tangent) > 0) {
	  modelContactForces(i,idx,dij);
	}
	interactions_(i,idx) = 1;
	interactions_(idx,i) = 1;
      }      
    }
  }
}

void ParticlePhysics::initializeParticleInteractionTracker() {
  interactions_ = MatrixXi::Zero(samples_,samples_);
}

void ParticlePhysics::particleContact() {
  // Calculate all particle-particle contacts (N^2 brute force)
  const MatrixXd& XY = simulation_->getXY();
  const MatrixXd& UV = simulation_->getUV();
  const VectorXd& R  = simulation_->getR();
  # pragma omp parallel for
  for (int i=0; i<samples_; i++) {
    for (int j=i+1; j<samples_; j++) {
      Vector2d dij        = XY.row(j)-XY.row(i);
      double distance_ij  = dij.norm();
      if (distance_ij <= R(i)+R(j)) {
	double ui_tangent = UV.row(i).dot(dij)/distance_ij;
	double uj_tangent = UV.row(j).dot(dij)/distance_ij;
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

const MatrixXd& ParticlePhysics::RHS() {
  zeroForces();
  particleContactRtree();
  double diff   = 3.77;//0.7;
  for (int i=0; i<samples_; i++) {
    forces_(i,0) *= diff/mass_(i);
    forces_(i,1) *= diff/mass_(i);
  }
  return forces_;
  
}
