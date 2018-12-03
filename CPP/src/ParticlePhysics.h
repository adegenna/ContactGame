#ifndef __PARTICLEPHYSICS_H__
#define __PARTICLEPHYSICS_H__

#include <Eigen/Dense>
#include "options.hpp"
#include "LagrangianState.h"
#include <iostream>
#include <memory>


class ContactForceModel{
public:
    virtual void particleContact(const LagrangianState& simulation,
                                 Eigen::MatrixXd& forces) = 0;
    virtual void particleWallContact(const LagrangianState& simulation,
				     const LagrangianState& boundary,
				     Eigen::MatrixXd& forces) = 0;

protected:
    Eigen::Vector2d modelContactForces(double Ri,
                                       double Rj,
                                       const Eigen::Vector2d& dij,
                                       const LagrangianState& simulation) const
    {
        double distance_ij = dij.norm();
        double eps         = 0.001;
        double delta = std::min( std::abs((Ri+Rj)-distance_ij) , eps*(Ri+Rj) );
        double F     = pow(10.0,5.0)*pow(delta,0.85);
        Eigen::VectorXd eij = dij/distance_ij;
	double diff  = 3.77; //0.7;
        return Eigen::Vector2d(diff*F*eij(0), diff*F*eij(1));
    }
};

class BruteForceContactForceModel : public ContactForceModel {
public:
    void particleContact(const LagrangianState &simulation,
                         Eigen::MatrixXd &forces) override;
    void particleWallContact(const LagrangianState& simulation,
			     const LagrangianState& boundary,
			     Eigen::MatrixXd& forces) override;
};

class RTreeContactForceModel : public ContactForceModel {
public:
    void particleContact(const LagrangianState &simulation,
                         Eigen::MatrixXd &forces) override;
    void particleWallContact(const LagrangianState& simulation,
			     const LagrangianState& boundary,
			     Eigen::MatrixXd& forces) override;
};

class ParticlePhysics {

 public:

  ParticlePhysics(Options& options, LagrangianState& simulation, bool useRtree = false);
  ParticlePhysics(Options& options, LagrangianState& simulation, LagrangianState& boundaryParticles, bool useRtree = false);
  ~ParticlePhysics();
//  void particleContactRtree();
  void calculateParticleMasses();
  double calculateTotalEnergy();
  Eigen::VectorXd calculateTotalMomentum();
  void writeEnergyAndMomentum(double energy, Eigen::VectorXd& momentum, const std::string& filename) const;
  void zeroForces();
  void initializeParticleInteractionTracker();
  const Eigen::MatrixXd& RHS();
  const Eigen::MatrixXd& getForces() const { return forces_; }
  
 private:
  const Options options_;
  LagrangianState* simulation_;
  LagrangianState* boundaryParticles_;
  int samples_;
  std::unique_ptr<ContactForceModel> contact_model_;
  Eigen::MatrixXd forces_;
  Eigen::VectorXd mass_;
  Eigen::MatrixXi interactions_;
  
};


#endif
