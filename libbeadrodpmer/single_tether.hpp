
#ifndef BEADRODPMER_SINGLE_TETHER_HPP
#define BEADRODPMER_SINGLE_TETHER_HPP

#include "polymer.hpp"
#include <vector>
#include <functional>

namespace BeadRodPmer {
class SingleTether : public Polymer {
public:
  // constructor
  SingleTether(const std::vector<std::string> &);
  

  int single_step(double,double,
		  const std::vector<std::vector<double>> &,
		  int itermax = 20, int numtries = 5,
		  bool throw_exception=true) override;

  
  int single_step(double,double,
		  const std::vector<std::vector<double>> &,
		  std::function<Eigen::Vector3d (double)>,
		  std::function<Eigen::Vector3d (double)>,
		  int itermax = 20, int numtries = 5,
		  bool throw_exception=true) ;
  void set_G() override;



  
  void set_Hhat() override;



  void set_dCdlambda() override;


  void test_jacob(int,double,const Eigen::Vector3d &);  

  void compute_noise() override;


  void compute_tension(const Eigen::Vector3d &);
  int correct_tension(double,const Eigen::Vector3d&,
		      int itermax = 20, double tol = 1e-14) ;

  void calculate_constraint_errors(const Eigen::Vector3d &) ;  

  void compute_effective_kappa() override;

private:
  

  std::vector<T> init_G_coeffsmatrix();

  

  std::vector<T> init_Hhat_coeffsmatrix();  

  std::vector<T> init_dCdlambda_coeffsmatrix();

  void set_rhs_of_G();
  void set_rhs_of_Hhat(const Eigen::Vector3d &) ;

  void update_G();
  void update_Hhat();  
  void update_dCdlambda(double);
  
};
};
#endif
