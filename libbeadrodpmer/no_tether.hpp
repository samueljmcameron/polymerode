
#ifndef BEADRODPMER_NO_TETHER_HPP
#define BEADRODPMER_NO_TETHER_HPP

#include "polymer.hpp"

namespace BeadRodPmer {
class NoTether : public Polymer {
public:
  // constructor
  NoTether(const std::vector<std::string> &);

  virtual int single_step(double, double,
			  const std::vector<Eigen::Vector3d> &,
			  int itermax = 20, int numtries = 5,
			  bool throw_exception=true) override ;


  virtual void setup() override;

  
  void test_jacob(int,double) ;


  int correct_tension(double,int itermax = 20,double tol = 1e-14) ;





  virtual void compute_noise() override;
  void compute_tension() ;

  virtual void compute_effective_kappa() override;

private:
  
  virtual void set_G();
  virtual void set_Hhat();
  virtual void set_dCdlambda();


  
};
};
#endif
