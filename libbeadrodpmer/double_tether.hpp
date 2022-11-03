
#ifndef BEADRODPMER_DOUBLE_TETHER_HPP
#define BEADRODPMER_DOUBLE_TETHER_HPP

#include "polymer.hpp"

namespace BeadRodPmer {
class DoubleTether : public Polymer {
public:
  // constructor
  DoubleTether(const std::vector<std::string> &);
  
  int single_step(double , double ,
		  const std::vector<std::vector<double>> & ,
		  int itermax = 20,int numtries = 5,
		  bool throw_exception = true) override;
 

  int single_step(double , double ,
		  const std::vector<std::vector<double>> & ,
		  std::function<Eigen::Vector3d (double)>,
		  std::function<Eigen::Vector3d (double)>,
		  std::function<Eigen::Vector3d (double)>,
		  std::function<Eigen::Vector3d (double)>,
		  int itermax = 20, int numtries = 5,
		  bool throw_exception = true);
  
  void set_G() override;
  void update_G() override;

  void set_Hhat() override;
  void update_Hhat() override;




  void set_dCdlambda() override;
  void update_dCdlambda(double) override;
  void test_jacob(int,double,const Eigen::Vector3d &, const Eigen::Vector3d &) ;  


  int correct_tension(double,const Eigen::Vector3d&,const Eigen::Vector3d&,
		      int itermax = 20, double tol = 1e-14) ;





  void compute_noise() override;
  void compute_tension(const Eigen::Vector3d &,
		       const Eigen::Vector3d &) ;

  
  void initial_integrate(double) override;
  void final_integrate(double) override;

  void calculate_constraint_errors(const Eigen::Vector3d &,
				   const Eigen::Vector3d &) ;  

  void compute_effective_kappa() override;

private:
  

  std::vector<T> init_G_coeffsmatrix() override;

  

  std::vector<T> init_Hhat_coeffsmatrix() override;  

  std::vector<T> init_dCdlambda_coeffsmatrix() override;

  void set_rhs_of_G() override;
  void set_rhs_of_Hhat(const Eigen::Vector3d &,
		       const Eigen::Vector3d &) ;
  
};
};
#endif
