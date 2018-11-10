#ifndef INPUTS_TEST_H_
#define INPUTS_TEST_H_

#include "gtest/gtest.h"
#include "math.h"

class InputsTest: public ::testing::Test {
 protected:
  virtual void SetUp() {
    X_.resize(3,2);
    X_ << 1.1, 2.2,
          3.3, 4.4,
          5.5, 6.6;
  }
  
  Eigen::MatrixXd X_;
};

#endif
