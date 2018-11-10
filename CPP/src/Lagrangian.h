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
  setupInitialConditions();

 private:

  std::string inputfile_;
  MatrixXd input_;
  int samples_;
  
};


#endif
