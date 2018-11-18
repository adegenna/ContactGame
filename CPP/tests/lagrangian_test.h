#ifndef LAGRANGIAN_TEST_H_
#define LAGRANGIAN_TEST_H_

#include "gtest/gtest.h"
#include "math.h"
#include <Eigen/Dense>

class LagrangianTest: public ::testing::Test {
 protected:
  virtual void SetUp() {
    X_.resize(4,2);
    X_ << 10.0, 10.0,
      -10.0, 10.0, 
      -10.0, -10.0,
      10.0, -10.0;
    DXY_.resize(4,2);
    DXY_ << 1.0, 1.0,
      -1.0, 1.0,
      -1.0, -1.0,
      1.0, -1.0;
  }
  
  Eigen::MatrixXd X_;
  Eigen::MatrixXd DXY_;
};

#endif
