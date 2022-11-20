
#ifndef BEADRODPMER_NO_TETHER_HPP
#define BEADRODPMER_NO_TETHER_HPP

#include "polymer.hpp"

namespace BeadRodPmer {
class NoTether : public Polymer {
public:
  // constructor
  NoTether(const std::vector<std::string> &);

  int single_step(double, double,
		  const std::vector<std::vector<double>> &,
		  int itermax = 20, int numtries = 5,
		  bool throw_exception=true) override ;

    
  void set_G() override ;
  
  void set_Hhat() override;

  void set_dCdlambda() override;
  
  void test_jacob(int,double) ;


  int correct_tension(double,int itermax = 20,double tol = 1e-14) ;





  void compute_noise() override;
  void compute_tension() ;

  void compute_effective_kappa() override;

private:
  

  void init_atoms_rand();
  void init_atoms_line();
  void init_atoms_caret();


  
};
};
#endif
