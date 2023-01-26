
#ifndef BEADRODPMER_NO_TETHER_HPP
#define BEADRODPMER_NO_TETHER_HPP

#include "polymer.hpp"

#include <memory>

namespace BeadRodPmer {
class NoTether : public Polymer {
public:
  // constructor
  NoTether(const std::vector<std::string> &);

  virtual void setup(const Eigen::Ref<const Eigen::Matrix3Xd> &) override;

  virtual void first_step(Eigen::Ref<Eigen::Matrix3Xd> ,
			  Eigen::Ref<Eigen::Matrix3Xd> ,double) override;

  virtual int second_step(Eigen::Ref<Eigen::Matrix3Xd> ,
			  Eigen::Ref<Eigen::Matrix3Xd> ,double,int) override;

  
  void test_jacob(int,double) ;

  
protected:

  SpMat Gmunu;      // geometric tensor
  Eigen::VectorXd rhs_of_G;  // rhs_of_G vector
  Eigen::VectorXd dummy_for_noise;


  SpMat Hhat;      // dynamical (hydrodynamic) tensor  
  Eigen::VectorXd rhs_of_Hhat;  // rhs_of_H vector AND dC vector
  Eigen::VectorXd tension;


  
  SpMat dCdlambda; // matrix for storing Jacobian of dC=0 Newton solver  
  Eigen::VectorXd negative_tension_change;


  std::shared_ptr<Eigen::SimplicialLDLT< SpMat, Eigen::Lower >> Gmunu_solver;
  std::shared_ptr<Eigen::SimplicialLDLT< SpMat, Eigen::Lower >> Hhat_solver ;

  std::shared_ptr<Eigen::SparseLU< SpMat >> jacob_solver;

  Eigen::VectorXd tDets;
  Eigen::VectorXd bDets;


  Eigen::VectorXd k_effs;
  Eigen::VectorXd end_inverses;

  
  int correct_tension(Eigen::Ref<Eigen::Matrix3Xd>,const Eigen::Ref<const Eigen::Matrix3Xd> &,
		      double,int itermax = 20,double tol = 1e-14) ;

  void compute_uc_forces(Eigen::Ref<Eigen::Matrix3Xd>);
  double Hhat_endblocks(int,int,int);
  double Hhat_leftside(int);
  double dCdlambda_leftside(int);
  double Hhat_bottomside(int);
  double dCdlambda_bottomside(int);



  void set_rhs_of_G(int);

  void set_rhs_of_Hhat(int,const Eigen::Ref<const Eigen::Matrix3Xd> &);
  void init_G_coeffsmatrix(int ,std::vector<T> &);
  void update_G(int);
  void update_Hhat(int);
  void update_dCdlambda(double,int) ;
  void init_Hhat_coeffsmatrix(int, std::vector<T>  &);

  
  void set_bdets_and_tdets(int);

  
  void init_dCdlambda_coeffsmatrix(int,std::vector<T> &);

  void update_noise(int);

  
  
  void initial_integrate(Eigen::Ref<Eigen::Matrix3Xd>,
			 const Eigen::Ref<const Eigen::Matrix3Xd>&,
			 double,int,PTYPE);
  void final_integrate(Eigen::Ref<Eigen::Matrix3Xd>,
		       const Eigen::Ref<const Eigen::Matrix3Xd>&,
		       double,int,PTYPE);

  void calculate_constraint_errors(int offset,const Eigen::Ref<const Eigen::Matrix3Xd>&);

  virtual void compute_noise();
  void compute_tension(const Eigen::Ref<const Eigen::Matrix3Xd>&) ;

  virtual void compute_effective_kappa();

  Eigen::Matrix3Xd tmp_xs;

private:
  
  virtual void set_G();
  virtual void set_Hhat();
  virtual void set_dCdlambda();


  double Hhat_diag_val(int);
  double Hhat_loweroff_val(int);
  
  double dCdlambda_diag_val(int);
  double dCdlambda_loweroff_val(int);
  double dCdlambda_upperoff_val(int);

  Eigen::Matrix3Xd tmp_bonds;
 
  
};
};
#endif
