#include <boost/program_options.hpp>
#include "physics_test.h"
#include "../src/Inputfile.hpp"
#include "../src/LagrangianState.h"
#include "../src/ParticlePhysics.h"

using namespace Eigen;

TEST_F(PhysicsTest, testEulerIntegrator) {
  
  Options options;
  options.inputfile  = std::string(SRCDIR)+"tests/testinput.csv";
  options.outputfile = "final.csv";
  options.dt         = 0.1;
  options.tsteps     = 10;
  options.tsave      = 20;
  
  // Setup
  LagrangianState simulation(load_csv<MatrixXd>(options.inputfile));
  ParticlePhysics physics(options, simulation);
  
  // Solve
  physics.simulate();

  // Output
  simulation.writeXY(options.outputfile);
  
  // Read the output
  MatrixXd out;
  out   = load_csv<MatrixXd>(options.outputfile);
  
  ASSERT_TRUE(out.isApprox(Xfinal_));
  
}

TEST_F(PhysicsTest, testBilliards) {

  Options options;
  options.inputfile  = std::string(SRCDIR)+"tests/billiards.csv";
  options.outputfile = "billiardsfinal";
  options.dt         = 0.001;
  options.tsteps     = 20000;
  options.tsave      = 500;
  
  // Setup
  LagrangianState simulation(load_csv<MatrixXd>(options.inputfile));
  ParticlePhysics physics(options, simulation);
  
  // Solve
  physics.simulate();

  // Output
  simulation.writeXY(options.outputfile);
  
  // Read the output
  MatrixXd out;
  out   = load_csv<MatrixXd>(options.outputfile);
  
  ASSERT_TRUE(out.isApprox(XYfinalBilliards_));
  
}
