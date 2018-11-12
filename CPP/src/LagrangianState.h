#ifndef __LAGRANGIANSTATE_H__
#define __LAGRANGIANSTATE_H__

#include <Eigen/Dense>
#include "options.hpp"

class LagrangianState {

 public:

  LagrangianState(const Eigen::MatrixXd& initialstate);
  ~LagrangianState();
  void updateXY(const Eigen::MatrixXd& DXY);
  void writeXY(const std::string& filename) const;
  void incrementUV(const Eigen::MatrixXd& DUV) { UV_ += DUV; }
  int getSamples() { return samples_; }
  const Eigen::MatrixXd& getXY() const { return XY_; }
  const Eigen::MatrixXd& getUV() const { return UV_; }
  const Eigen::VectorXd& getR()  const { return R_;  }
  
 private:

  std::string inputfile_;
  std::string projdir_;
  std::string outdir_;
  std::string loaddir_;
  Eigen::MatrixXd input_;
  Eigen::MatrixXd XY_;
  Eigen::MatrixXd UV_;
  Eigen::VectorXd R_;
  int samples_;
  
};


#endif
