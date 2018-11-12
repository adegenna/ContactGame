#ifndef __LAGRANGIANSIMULATION_H__
#define __LAGRANGIANSIMULATION_H__

#include <Eigen/Dense>
#include "options.hpp"


// ***************************************************
// CLASS FOR LAGRANGIAN
// ***************************************************

class LagrangianSimulation {

 public:

  LagrangianSimulation();
  LagrangianSimulation(Options& options);
  ~LagrangianSimulation();
  void setupInitialConditions();
  void updateXY(Eigen::MatrixXd& DXY);
  void writeXY(std::string append);
  void incrementUV(Eigen::MatrixXd& DUV) { UV_ += DUV; }
  int getSamples() { return samples_; }
  Eigen::MatrixXd getXY() { return XY_; }
  Eigen::MatrixXd getUV() { return UV_; }
  Eigen::VectorXd getR()  { return R_;  }
  
 private:

  std::string inputfile_;
  std::string projdir_;
  std::string outdir_;
  std::string loaddir_;
  Eigen::MatrixXd input_;
  Eigen::MatrixXd XY_;
  Eigen::MatrixXd UV_;
  Eigen::VectorXd R_;
  int samples_;
  
};


#endif
