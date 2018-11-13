#ifndef PHYSICS_TEST_H_
#define PHYSICS_TEST_H_

#include "gtest/gtest.h"
#include "math.h"
#include <Eigen/Dense>

class PhysicsTest: public ::testing::Test {
 protected:
  virtual void SetUp() {
    Xfinal_.resize(4,2);
    Xfinal_ << 11.0, 11.0,
      -11.0, 11.0,
      -11.0, -11.0,
      11.0, -11.0;
  }

  Eigen::MatrixXd Xfinal_;
  
};

#endif
