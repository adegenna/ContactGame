#include <boost/program_options.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>
#include <random>
#include "quadtree_test.h"
#include "../src/Inputfile.hpp"
#include "../src/LagrangianState.h"
#include "../src/ParticlePhysics.h"
#include "../src/TimeIntegration.h"
#include "../src/Quadtree/Bucket.h"

using namespace Eigen;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

BOOST_GEOMETRY_REGISTER_POINT_2D(Eigen::Vector2d, double, boost::geometry::cs::cartesian, Eigen::Vector2d::x(), Eigen::Vector2d::y())

TEST_F(QuadtreeTest, testQuadtreeEigen) {
  
  Options options;
  options.inputfile  = std::string(SRCDIR)+"tests/billiards.csv";
  options.outputfile = "final.csv";
  options.dt         = 0.1;
  options.tsteps     = 10;
  options.tsave      = 20;
  
  // Setup
  //typedef bg::model::box<Vector2d> box;
  typedef std::pair<Vector2d, unsigned> value;
  MatrixXd XYUVR = load_csv<MatrixXd>(options.inputfile);
  
  // create the rtree using default constructor
  bgi::rtree< value, bgi::quadratic<16> > rtree;
  // create some values
  for (int i=0; i<XYUVR.rows(); i++) {
    // create a box
    Vector2d a; a[0] = XYUVR(i,0); a[1] = XYUVR(i,1);
    // insert new value
    rtree.insert(std::make_pair(a, i));
  }
  // find nearest values to a point
  std::vector<value> result_n;
  Vector2d q; q << 0. , 0.;  
  rtree.query(bgi::nearest(q, 1), std::back_inserter(result_n));

  // display results
  // std::cout << "knn query point:" << std::endl;
  // Vector2d a; a << 0. , 0.;
  // std::cout << bg::wkt<Vector2d>(a) << std::endl;
  // std::cout << "knn query result:" << std::endl;
  // BOOST_FOREACH(value const& v, result_n)
  //   std::cout << bg::wkt<Vector2d>(v.first) << " - " << v.second << std::endl;

  // Actual answer
  Vector2d xy; xy << XYUVR(0,0) , XYUVR(0,1);
  value act = std::make_pair(xy,0);
  
  ASSERT_TRUE(result_n[0].first.isApprox(xy));
  
}

TEST_F(QuadtreeTest, testQuadtreeCustom) {
  
  Options options;
  options.inputfile  = std::string(SRCDIR)+"tests/billiards.csv";
  options.outputfile = "final.csv";
  options.dt         = 0.1;
  options.tsteps     = 10;
  options.tsave      = 20;
  const std::string outdir = std::string(SRCDIR)+"build";
  
  // Initialize first bucket
  double SW[2] = {0.,0.};
  double SE[2] = {1.,0.};
  double NW[2] = {0.,1.};
  double NE[2] = {1.,1.};
  Bucket* QT = new Bucket(&SW[0],&SE[0],&NW[0],&NE[0]);
  QT->setOutDir(outdir);
  
  // Test set/get points
  const int nrolls=1000;  // number of experiments

  std::default_random_engine generator;
  std::normal_distribution<double> distX(0.3,.1);
  std::normal_distribution<double> distY(0.4,0.2);
  FILE* fout = fopen("../build/DataSet.dat","w");
  
  double sampsX[nrolls];
  double sampsY[nrolls];
  for (int i=0; i<nrolls; i++) {
    double nx = distX(generator);
    double ny = distY(generator);
    sampsX[i] = nx;
    sampsY[i] = ny;   
    fprintf(fout,"%f\t%f\n",nx,ny);
  }

  // Divide buckets
  QT->calcQuadTree(&sampsX[0],&sampsY[0],nrolls);

  // Search for a query point
  double Xq = 0.41;
  double Yq = 0.5;
  double Xnn, Ynn;
  int indnn;
  Eigen::MatrixXd xyi;
  QT->nnSearch(&Xq,&Yq,&Xnn,&Ynn,&indnn);
  QT->knnSearch(&Xq,&Yq,xyi);
  
  printf("Xq = %f, Yq = %f\nXnn = %f, Ynn = %f\n",Xq,Yq,Xnn,Ynn);
  std::cout << "K-nearest neighbors: " << std::endl;
  std::cout << xyi;
  
  // Delete allocated memory
  delete QT;
  fclose(fout);

  EXPECT_TRUE(std::abs(Xq-Xnn) < 0.1);
  EXPECT_TRUE(std::abs(Yq-Ynn) < 0.1);
  
}


