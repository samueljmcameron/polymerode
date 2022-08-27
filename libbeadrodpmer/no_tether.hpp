
#ifndef BEADRODPMER_NO_TETHER_HPP
#define BEADRODPMER_NO_TETHER_HPP

#include "polymer.hpp"

namespace BeadRodPmer {
class NoTether : public Polymer {
public:
  // constructor
  NoTether(const std::vector<std::string> &);
  ~NoTether();  



  void single_step(double, double,
		   const std::vector<std::vector<double>> &,
		   int itermax = 20, int numtries = 5);

  
  void init_atoms();
  
  void set_G();
  void update_G();


  
  void set_Hhat();
  void update_Hhat();




  void set_dCdlambda();
  void update_dCdlambda(double);
  void test_jacob(int,double);  


  int correct_tension(double,int itermax = 20,double tol = 1e-14);





  void compute_noise();
  void compute_tension();

  
  void initial_integrate(double);
  void final_integrate(double);

  void calculate_constraint_errors();  

  void compute_effective_kappa();

private:
  

  std::vector<T> init_G_coeffsmatrix();

  

  std::vector<T> init_Hhat_coeffsmatrix();  

  std::vector<T> init_dCdlambda_coeffsmatrix();

  void set_rhs_of_G();
  void set_rhs_of_Hhat();
  
};
};
#endif
