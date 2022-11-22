
#ifndef BEADRODPMER_SINGLE_TETHER_HPP
#define BEADRODPMER_SINGLE_TETHER_HPP

#include "no_tether.hpp"
#include <vector>
#include <functional>

namespace BeadRodPmer {
class SingleTether : public NoTether {
public:
  // constructor
  SingleTether(const std::vector<std::string> &);
  

  virtual int single_step(double,double,
			  const std::vector<Eigen::Vector3d> &,
			  int itermax = 20, int numtries = 5,
			  bool throw_exception=true) override;

  
  int single_step(double,double,
		  const std::vector<Eigen::Vector3d> &,
		  std::function<Eigen::Vector3d (double)>,
		  std::function<Eigen::Vector3d (double)>,
		  int itermax = 20, int numtries = 5,
		  bool throw_exception=true) ;


  virtual void compute_noise() override;
  
  void test_jacob(int,double,const Eigen::Vector3d &);  




  void compute_tension(const Eigen::Vector3d &);
  int correct_tension(double,const Eigen::Vector3d&,
		      int itermax = 20, double tol = 1e-14) ;

  void calculate_constraint_errors(const Eigen::Vector3d &) ;  


  void set_bdets_and_tdets();
  virtual void compute_effective_kappa() override;
  
  virtual std::vector<T> init_G_coeffsmatrix();
  virtual std::vector<T> init_Hhat_coeffsmatrix();  
  virtual std::vector<T> init_dCdlambda_coeffsmatrix();

  virtual void update_G();
  virtual void update_Hhat();  
  virtual void update_dCdlambda(double);

  
  virtual void set_rhs_of_G();
  
private:
  

  void set_rhs_of_Hhat(const Eigen::Vector3d &) ;

  virtual void set_G() override;
  virtual void set_Hhat() override;
  virtual void set_dCdlambda() override;
  
};
};
#endif
