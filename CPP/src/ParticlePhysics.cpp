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

ParticlePhysics::ParticlePhysics(Options& o, LagrangianState& simulation, bool useRtree)
  : options_(o), simulation_(&simulation)
{
  samples_ = simulation_->getSamples();
  if (useRtree){
    contact_model_ = std::make_unique<RTreeContactForceModel>();
  } else {
    contact_model_ = std::make_unique<BruteForceContactForceModel>();
  }
  calculateParticleMasses();
}

ParticlePhysics::~ParticlePhysics() {
  
}


void ParticlePhysics::initializeParticleInteractionTracker() {
  interactions_ = MatrixXi::Zero(samples_,samples_);
}



void
BruteForceContactForceModel::
particleContact(const LagrangianState &simulation,
                Eigen::MatrixXd &forces)
{
  // Calculate all particle-particle contacts (N^2 brute force)
  const MatrixXd& XY = simulation.getXY();
  const MatrixXd& UV = simulation.getUV();
  const VectorXd& R  = simulation.getR();
  const int samples  = simulation.getSamples();

  forces.resizeLike(XY);
  forces.array() = 0;

# pragma omp parallel for
  for (int i=0; i<samples; i++) {
    for (int j=i+1; j<samples; j++) {
      Vector2d dij        = XY.row(j)-XY.row(i);
      double distance_ij  = dij.norm();
      if (distance_ij <= R(i)+R(j)) {
        double ui_tangent = UV.row(i).dot(dij)/distance_ij;
        double uj_tangent = UV.row(j).dot(dij)/distance_ij;
        if ((ui_tangent - uj_tangent) > 0) {
          const auto local_force = modelContactForces(i,j,dij,simulation);

          forces(i,0) -= local_force(0);
          forces(i,1) -= local_force(1);
          forces(j,0) += local_force(0);
          forces(j,1) += local_force(1);
        }
      }
    }
  }
}


void
RTreeContactForceModel::
particleContact(const LagrangianState &simulation,
                Eigen::MatrixXd &forces)
{
  const MatrixXd& XY = simulation.getXY();
  const MatrixXd& UV = simulation.getUV();
  const VectorXd& R  = simulation.getR();
  const int samples  = simulation.getSamples();

  forces.resizeLike(XY);
  forces.array() = 0;

  // Rtree construction
  typedef std::pair<Vector2d, unsigned> value;
  bgi::rtree< value, bgi::quadratic<16> > rtree;
  for (int i=0; i<samples; i++) {
    rtree.insert(std::make_pair(Vector2d(XY(i,0), XY(i,1)), i));
  }


  for (int i=0; i<samples; i++) {

    std::vector<value> contact_points;
    const Vector2d search_point(XY(i,0), XY(i,1));

    // rtree search using a lambda function that returns true for the distance
    // between the points less than the sum of their radii.
    rtree.query(bgi::satisfies([&](value const& v) {
      return bg::distance(v.first, search_point) <= R(i)+R(v.second); }),
      std::back_inserter(contact_points));

    for (const auto& found_point : contact_points){
      const int idx = found_point.second;
      const Vector2d dij( found_point.first[0] - XY(i,0),
                          found_point.first[1] - XY(i,1) );

      const double distance_ij = dij.norm();
      const double ui_tangent  = UV.row(i).dot(dij)/distance_ij;
      const double uj_tangent  = UV.row(idx).dot(dij)/distance_ij;

      if (((ui_tangent - uj_tangent) > 0) && idx > i) {
        const auto local_force = modelContactForces(i,idx,dij,simulation);

        forces(i,0)   -= local_force(0);
        forces(i,1)   -= local_force(1);
        forces(idx,0) += local_force(0);
        forces(idx,1) += local_force(1);
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
//  particleContactRtree();
  contact_model_->particleContact(*simulation_, forces_);
  double diff   = 3.77;//0.7;
  for (int i=0; i<samples_; i++) {
    forces_(i,0) *= diff/mass_(i);
    forces_(i,1) *= diff/mass_(i);
  }
  return forces_;
  
}
