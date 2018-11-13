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
    XYfinalBilliards_.resize(11,2);
    XYfinalBilliards_ << 0, 1.88628,
      -3.20498, 3.07456,
      3.20498, 3.07456,
      -4.07322, 5.20091,
      0, 4.44799,
      4.07322, 5.20091,
      -4.61692, 7.9969,
      -1.00566, 6.86574,
      1.00566, 6.86574,
      4.61692, 7.9969,
      0, -1.96947;
  }

  Eigen::MatrixXd Xfinal_;
  Eigen::MatrixXd XYfinalBilliards_;
  
};

#endif
