
#ifndef BEADRODPMER_NO_TETHER_HPP
#define BEADRODPMER_NO_TETHER_HPP

#include "polymer.hpp"

namespace BeadRodPmer {
class NoTether : public Polymer {
public:
  // constructor
  NoTether(const std::vector<std::string> &);

  virtual int single_step(double, double,
			  const std::vector<Eigen::Vector3d> &,
			  int itermax = 20, int numtries = 5,
			  bool throw_exception=true) override ;


  virtual void setup() override;

  
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


  Eigen::SimplicialLDLT< SpMat, Eigen::Lower > Gmunu_solver;
  Eigen::SimplicialLDLT< SpMat, Eigen::Lower > Hhat_solver;

  Eigen::SparseLU< SpMat > jacob_solver;

  Eigen::VectorXd tDets;
  Eigen::VectorXd bDets;
  
  int correct_tension(double,int itermax = 20,double tol = 1e-14) ;

  void compute_uc_forces();
  double Hhat_endblocks(int,int,int);
  double Hhat_leftside(int);
  double dCdlambda_leftside(int);
  double Hhat_bottomside(int);
  double dCdlambda_bottomside(int);



  void set_rhs_of_G(int);

  void set_rhs_of_Hhat(int);
  void init_G_coeffsmatrix(int ,std::vector<T> &);
  void update_G(int);
  void update_Hhat(int);
  void update_dCdlambda(double,int) ;
  void init_Hhat_coeffsmatrix(int, std::vector<T>  &);

  
  void set_bdets_and_tdets(int);

  
  void init_dCdlambda_coeffsmatrix(int,std::vector<T> &);

  void update_noise(int);

  
  
  void initial_integrate(double,int,PTYPE);
  void final_integrate(double,int,PTYPE);

  void calculate_constraint_errors(int offset);




  virtual void compute_noise();
  void compute_tension() ;

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
