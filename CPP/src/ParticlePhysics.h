#ifndef __PARTICLEPHYSICS_H__
#define __PARTICLEPHYSICS_H__

#include <Eigen/Dense>
#include "options.hpp"
#include "LagrangianSimulation.h"


// ***************************************************
// CLASS FOR PHYSICS
// ***************************************************

class ParticlePhysics {

 public:

  ParticlePhysics();
  ParticlePhysics(Options& options, LagrangianSimulation& simulation);
  ~ParticlePhysics();
  void modelContactForces(int i, int j, Eigen::VectorXd& dij);
  void particleContact();
  void updateParticleVelocities();
  void calculateParticleMasses();
  void eulerDXY(Eigen::MatrixXd& DXY);
  void simulate();
  
 private:
  LagrangianSimulation* simulation_;
  double dt_;
  std::string integrator_;
  int tsteps_;
  int samples_;
  int tsave_;
  Eigen::MatrixXd forces_;
  Eigen::VectorXd mass_;
  
  
};


#endif
