#include "inputs_test.h"
#include "../src/Inputfile.hpp"

using namespace Eigen;

TEST_F(InputsTest, InputReading) {
  // Test reading from input file
  MatrixXd A = load_csv<MatrixXd>(std::string(SRCDIR)+"tests/testinput.csv");
  ASSERT_TRUE(A.block(0,0,3,2).isApprox(X_));
  
}
