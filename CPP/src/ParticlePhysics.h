#ifndef __PARTICLEPHYSICS_H__
#define __PARTICLEPHYSICS_H__

#include <eigen/Eigen/Dense>
#include "options.hpp"
#include "LagrangianSimulation.h"

using namespace std;
using namespace Eigen;

// ***************************************************
// CLASS FOR PHYSICS
// ***************************************************

class ParticlePhysics {

 public:

  ParticlePhysics();
  ParticlePhysics(Options& options, LagrangianSimulation& simulation);
  ~ParticlePhysics();
  void eulerDXY(MatrixXd& DXY);
  void simulate();
  
 private:
  LagrangianSimulation simulation_;
  double dt_;
  std::string integrator_;
  int tsteps_;
  
};


#endif
