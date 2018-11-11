#ifndef PHYSICS_TEST_H_
#define PHYSICS_TEST_H_

#include "gtest/gtest.h"
#include "math.h"
#include <eigen/Eigen/Dense>

class PhysicsTest: public ::testing::Test {
 protected:
  virtual void SetUp() {
    Xfinal_.resize(3,2);
    Xfinal_ << 1.1, 2.2,
      3.3, 4.4,
      5.5, 6.6;
  }

  Eigen::MatrixXd Xfinal_;
  
};

#endif
