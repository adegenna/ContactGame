#ifndef INPUTS_TEST_H_
#define INPUTS_TEST_H_

#include "gtest/gtest.h"
#include "math.h"
#include <Eigen/Dense>

class InputsTest: public ::testing::Test {
 protected:
  virtual void SetUp() {
    X_.resize(4,5);
    X_ << 10.0, 10.0, 1.0, 1.0, 1.0,
      -10.0, 10.0, -1.0, 1.0, 1.0,
      -10.0, -10.0, -1.0, -1.0, 1.0,
      10.0, -10.0, 1.0, -1.0, 1.0;
  }
  
  Eigen::MatrixXd X_;
};

#endif
