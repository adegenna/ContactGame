#include <boost/program_options.hpp>
#include "lagrangian_test.h"
#include "../src/Inputfile.hpp"
#include "../src/LagrangianState.h"

using namespace Eigen;

TEST_F(LagrangianTest, testUpdateXY) {
  Options options;
  options.inputfile  = std::string(SRCDIR)+"tests/testinput.csv";
  options.outputfile = "final.csv";

  LagrangianState solver(load_csv<MatrixXd>(options.inputfile));

  solver.updateXY(DXY_);
  solver.writeXY(options.outputfile);

  // Read the output
  MatrixXd out;
  out   = load_csv<MatrixXd>(options.outputfile);
  
  ASSERT_TRUE(out.isApprox(X_+DXY_));
  
}
