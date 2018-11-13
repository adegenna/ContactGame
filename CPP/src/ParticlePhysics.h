#ifndef __PARTICLEPHYSICS_H__
#define __PARTICLEPHYSICS_H__

#include <Eigen/Dense>
#include "options.hpp"
#include "LagrangianState.h"

class ParticlePhysics {

 public:

  ParticlePhysics(Options& options, LagrangianState& simulation);
  ~ParticlePhysics();
  void modelContactForces(int i, int j, Eigen::VectorXd& dij);
  void particleContact();
  void updateParticleVelocities();
  void calculateParticleMasses();
  void eulerDXY(Eigen::MatrixXd& DXY);
  void simulate();
  
 private:
  const Options options_;
  LagrangianState* simulation_;
  int samples_;
  Eigen::MatrixXd forces_;
  Eigen::VectorXd mass_;
  
};


#endif
