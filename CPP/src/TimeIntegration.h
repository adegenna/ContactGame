#ifndef __TIMEINTEGRATION_H__
#define __TIMEINTEGRATION_H__

#include <Eigen/Dense>
#include "options.hpp"
#include "LagrangianState.h"
#include "ParticlePhysics.h"

class TimeIntegration {

 public:
  TimeIntegration(Options& options, ParticlePhysics& physics, LagrangianState& state);
  ~TimeIntegration();
  void euler();
  
 private:
  const Options options_;
  ParticlePhysics* physics_;
  LagrangianState* state_;
  
};


#endif
