#ifndef __LAGRANGIAN_H__
#define __LAGRANGIAN_H__

#include <eigen/Eigen/Dense>
#include "options.hpp"

using namespace std;
using namespace Eigen;

// ***************************************************
// CLASS FOR LAGRANGIAN
// ***************************************************

class Lagrangian {

 public:
  
  Lagrangian();
  Lagrangian(Options& options);
  ~Lagrangian();
  void setupInitialConditions();
  void updateXY(MatrixXd& DXY);
  void writeXY();
  void eulerDXY(MatrixXd& DXY);
  void integrate();
  
 private:

  std::string inputfile_;
  std::string projdir_;
  std::string outdir_;
  std::string loaddir_;
  std::string integrator_;
  MatrixXd input_;
  MatrixXd XY_;
  MatrixXd UV_;
  VectorXd R_;
  int samples_;
  double dt_;
  int tsteps_;
  
};


#endif
