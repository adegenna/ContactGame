#include <boost/program_options.hpp>
#include "lagrangian_test.h"
#include "../src/Inputfile.hpp"
#include "../src/LagrangianSimulation.h"

using namespace Eigen;

TEST_F(LagrangianTest, testUpdateXY) {
  Options options;
  options.projDir   = "/home/adegennaro/ContactGame/CPP/";
  options.outDir    = "tests/";
  options.loadDir   = "tests/";
  options.inputfile = "testinput.csv";

  LagrangianSimulation solver(options);

  solver.setupInitialConditions();
  solver.updateXY(DXY_);
  solver.writeXY();

  // Read the output
  MatrixXd out;
  out   = load_csv<MatrixXd>(options.projDir + options.outDir + "XY.csv");
  
  ASSERT_TRUE(out.isApprox(X_+DXY_));
  
}
