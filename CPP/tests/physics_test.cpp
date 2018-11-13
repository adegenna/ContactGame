#include <boost/program_options.hpp>
#include "physics_test.h"
#include "../src/Inputfile.hpp"
#include "../src/LagrangianState.h"
#include "../src/ParticlePhysics.h"

using namespace Eigen;

TEST_F(PhysicsTest, testEulerIntegrator) {

  Options options;
  options.projDir   = "/home/adegennaro/ContactGame/CPP/";
  options.outDir    = "tests/";
  options.loadDir   = "tests/";
  options.inputfile = "testinput.csv";
  options.dt        = 0.1;
  options.tsteps    = 10;
  
  // Setup
  LagrangianState simulation(load_csv<MatrixXd>(options.projDir + options.loadDir + options.inputfile));
  ParticlePhysics physics(options, simulation);
  
  // Solve
  physics.simulate();

  // Output
  simulation.writeXY("final");
  
  // Read the output
  MatrixXd out;
  out   = load_csv<MatrixXd>(options.projDir + options.outDir + "XY_final.csv");
  
  ASSERT_TRUE(out.isApprox(Xfinal_));
  
}
