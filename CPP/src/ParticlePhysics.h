#ifndef __PARTICLEPHYSICS_H__
#define __PARTICLEPHYSICS_H__

#include <Eigen/Dense>
#include "options.hpp"
#include "LagrangianState.h"

class ParticlePhysics {

 public:

  ParticlePhysics(Options& options, LagrangianState& simulation);
  ~ParticlePhysics();
  void modelContactForces(int i, int j, Eigen::Vector2d& dij);
  void particleContact();
  void particleContactRtree();
  void calculateParticleMasses();
  double calculateTotalEnergy();
  Eigen::VectorXd calculateTotalMomentum();
  void writeEnergyAndMomentum(double energy, Eigen::VectorXd& momentum, const std::string& filename) const;
  void zeroForces();
  const Eigen::MatrixXd& RHS();
  
 private:
  const Options options_;
  LagrangianState* simulation_;
  int samples_;
  Eigen::MatrixXd forces_;
  Eigen::VectorXd mass_;
  
};


#endif
