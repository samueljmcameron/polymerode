
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
  


  virtual void first_step(Eigen::Ref<Eigen::Matrix3Xd> ,
			  Eigen::Ref<Eigen::Matrix3Xd> ,double) override;

  virtual int second_step(Eigen::Ref<Eigen::Matrix3Xd> ,
			  Eigen::Ref<Eigen::Matrix3Xd> ,double,int) override;



  void first_step(Eigen::Ref<Eigen::Matrix3Xd> ,
			  Eigen::Ref<Eigen::Matrix3Xd> ,double,
			  const Eigen::Ref<const Eigen::Vector3d> &);

  int second_step(Eigen::Ref<Eigen::Matrix3Xd> ,
			  Eigen::Ref<Eigen::Matrix3Xd> ,double,int,
			  const Eigen::Ref<const Eigen::Vector3d> &,
			  const Eigen::Ref<const Eigen::Vector3d> &);

  virtual void compute_noise() override;
  
  void test_jacob(int,double,const Eigen::Vector3d &);  




  void compute_tension(const Eigen::Vector3d &,
		       const Eigen::Ref<const Eigen::Matrix3Xd>&);
  int correct_tension(Eigen::Ref<Eigen::Matrix3Xd>,
		      const Eigen::Ref<const Eigen::Matrix3Xd>&,
		      double,const Eigen::Vector3d&,
		      int itermax = 20, double tol = 1e-14) ;

  void calculate_constraint_errors(const Eigen::Vector3d &,
				   const Eigen::Ref<const Eigen::Matrix3Xd>&) ;  


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
  
  void set_rhs_of_Hhat(const Eigen::Vector3d &,
		       const Eigen::Ref<const Eigen::Matrix3Xd>&) ;

  virtual void set_G() override;
  virtual void set_Hhat() override;
  virtual void set_dCdlambda() override;
  
};
};
#endif
