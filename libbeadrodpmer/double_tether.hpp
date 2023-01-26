
#ifndef BEADRODPMER_DOUBLE_TETHER_HPP
#define BEADRODPMER_DOUBLE_TETHER_HPP

#include "single_tether.hpp"

namespace BeadRodPmer {
class DoubleTether : public SingleTether {
public:
  // constructor
  DoubleTether(const std::vector<std::string> &);
  
  virtual void first_step(Eigen::Ref<Eigen::Matrix3Xd> ,
			  Eigen::Ref<Eigen::Matrix3Xd> ,double) override;

  virtual int second_step(Eigen::Ref<Eigen::Matrix3Xd> ,
			  Eigen::Ref<Eigen::Matrix3Xd> ,double,int) override;

  
  void first_step(Eigen::Ref<Eigen::Matrix3Xd> ,
		  Eigen::Ref<Eigen::Matrix3Xd> ,double,
		  const Eigen::Ref<const Eigen::Vector3d> &,
		  const Eigen::Ref<const Eigen::Vector3d> &);

  int second_step(Eigen::Ref<Eigen::Matrix3Xd> ,
		  Eigen::Ref<Eigen::Matrix3Xd> ,double,int,
		  const Eigen::Ref<const Eigen::Vector3d> &,
		  const Eigen::Ref<const Eigen::Vector3d> &,
		  const Eigen::Ref<const Eigen::Vector3d> &,
		  const Eigen::Ref<const Eigen::Vector3d> &);


  virtual void compute_noise() final;

  
  void test_jacob(int,double,const Eigen::Vector3d &, const Eigen::Vector3d &) ;  


  int correct_tension(Eigen::Ref<Eigen::Matrix3Xd>,
		      const Eigen::Ref<const Eigen::Matrix3Xd>&,
		      double,const Eigen::Vector3d&,const Eigen::Vector3d&,
		      int itermax = 20, double tol = 1e-14) ;



  void compute_tension(const Eigen::Vector3d &,
		       const Eigen::Vector3d &,
		       const Eigen::Ref<const Eigen::Matrix3Xd>&) ;

  
  void calculate_constraint_errors(const Eigen::Vector3d &,
				   const Eigen::Vector3d &,
				   const Eigen::Ref<const Eigen::Matrix3Xd>&) ;  

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
		       const Eigen::Vector3d &,
		       const Eigen::Ref<const Eigen::Matrix3Xd>&) ;
  
};
};
#endif
