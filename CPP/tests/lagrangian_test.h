#ifndef LAGRANGIAN_TEST_H_
#define LAGRANGIAN_TEST_H_

#include "gtest/gtest.h"
#include "math.h"
#include <Eigen/Dense>

class LagrangianTest: public ::testing::Test {
 protected:
  virtual void SetUp() {
    X_.resize(3,2);
    X_ << 1.1, 2.2,
          3.3, 4.4,
          5.5, 6.6;
    DXY_.resize(3,2);
    DXY_ << 11.1, 22.2,
          33.3, 44.4,
          55.5, 66.6;
  }
  
  Eigen::MatrixXd X_;
  Eigen::MatrixXd DXY_;
};

#endif
