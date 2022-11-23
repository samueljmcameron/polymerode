
#ifndef BEADRODPMER_DOUBLE_TETHER_HPP
#define BEADRODPMER_DOUBLE_TETHER_HPP

#include "single_tether.hpp"

namespace BeadRodPmer {
class DoubleTether : public SingleTether {
public:
  // constructor
  DoubleTether(const std::vector<std::string> &);
  
  virtual int single_step(double , double ,
			  const std::vector<Eigen::Vector3d> & ,
			  int itermax = 20,int numtries = 5,
			  bool throw_exception = true) final;
 

  int single_step(double , double ,
		  const std::vector<Eigen::Vector3d> & ,
		  std::function<Eigen::Vector3d (double)>,
		  std::function<Eigen::Vector3d (double)>,
		  std::function<Eigen::Vector3d (double)>,
		  std::function<Eigen::Vector3d (double)>,
		  int itermax = 20, int numtries = 5,
		  bool throw_exception = true);

  virtual void compute_noise() final;

  
  void test_jacob(int,double,const Eigen::Vector3d &, const Eigen::Vector3d &) ;  


  int correct_tension(double,const Eigen::Vector3d&,const Eigen::Vector3d&,
		      int itermax = 20, double tol = 1e-14) ;



  void compute_tension(const Eigen::Vector3d &,
		       const Eigen::Vector3d &) ;

  
  void calculate_constraint_errors(const Eigen::Vector3d &,
				   const Eigen::Vector3d &) ;  

  virtual void compute_effective_kappa() final;

private:
  

  virtual std::vector<T> init_G_coeffsmatrix() final;
  virtual std::vector<T> init_Hhat_coeffsmatrix() final;  
  virtual std::vector<T> init_dCdlambda_coeffsmatrix();

  virtual void update_G() final;  
  virtual void update_Hhat() final;
  virtual void update_dCdlambda(double) final;

  

  virtual void set_rhs_of_G() final;
  void set_rhs_of_Hhat(const Eigen::Vector3d &,
		       const Eigen::Vector3d &) ;
  
};
};
#endif
