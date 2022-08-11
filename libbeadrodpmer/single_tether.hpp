
#ifndef BEADRODPMER_SINGLE_TETHER_HPP
#define BEADRODPMER_SINGLE_TETHER_HPP

#include "polymer.hpp"

namespace BeadRodPmer {
class SingleTether : public Polymer {
public:
  // constructor
  SingleTether(const std::vector<std::string> &);
  ~SingleTether();  
  
  


  void set_G();
  void update_G();


  
  void set_Hhat();
  void update_Hhat();




  void set_dCdlambda();
  void update_dCdlambda(double);
  void test_jacob(int,double,const Eigen::Vector3d &);  


  int correct_tension(double,const Eigen::Vector3d&,
		      int itermax = 20, double tol = 1e-14);





  void compute_noise();
  void compute_tension(const Eigen::Vector3d &);

  
  void initial_integrate(double);
  void final_integrate(double);

  void calculate_constraint_errors(const Eigen::Vector3d &);  

  void compute_effective_kappa();

private:
  

  std::vector<T> init_G_coeffsmatrix();

  

  std::vector<T> init_Hhat_coeffsmatrix();  

  std::vector<T> init_dCdlambda_coeffsmatrix();

  void set_rhs_of_G();
  void set_rhs_of_Hhat(const Eigen::Vector3d &);
  
};
};
#endif
