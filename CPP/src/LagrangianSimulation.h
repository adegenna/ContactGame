#ifndef __LAGRANGIANSIMULATION_H__
#define __LAGRANGIANSIMULATION_H__

#include <eigen/Eigen/Dense>
#include "options.hpp"

using namespace std;
using namespace Eigen;

// ***************************************************
// CLASS FOR LAGRANGIAN
// ***************************************************

class LagrangianSimulation {

 public:

  LagrangianSimulation();
  LagrangianSimulation(Options& options);
  ~LagrangianSimulation();
  void setupInitialConditions();
  void updateXY(MatrixXd& DXY);
  void writeXY(std::string append);
  int getSamples() { return samples_; }
  MatrixXd getUV() { return UV_; } 
  
 private:

  std::string inputfile_;
  std::string projdir_;
  std::string outdir_;
  std::string loaddir_;
  MatrixXd input_;
  MatrixXd XY_;
  MatrixXd UV_;
  VectorXd R_;
  int samples_;
  
};


#endif
