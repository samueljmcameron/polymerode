
#ifndef POLYMER_HPP
#define POLYMER_HPP

#include "atom.hpp"
#include "bond.hpp"


#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <vector>
#include <string>
#include <random>


typedef Eigen::SparseMatrix<double> SpMat; // declares column-major
typedef Eigen::Triplet<double> T;

class Polymer {
public:
  // constructor
  Polymer(std::vector<std::string>);
  ~Polymer();  
  
  std::vector<Eigen::Vector3d> Rtmp;    // temporary storate for midstep

  
  Eigen::VectorXd dummy_for_noise, tension,negative_tension_change;
  Eigen::VectorXd constraint_errors;
  double omega, rad;

  void rescale_positions(bool);

  virtual void init_atoms();
  void compute_tangents_and_friction();

  void set_unprojected_noise(double);

  
  void test_jacob(int,double,double);  

  Eigen::VectorXd rhs_of_G;  // rhs_of_G vector
  void set_rhs_of_G();

  SpMat Gmunu;      // geometric tensor
  void set_G();
  void update_G();


  void compute_uc_forces();

  
  Eigen::VectorXd rhs_of_Hhat;  // rhs_of_H vector AND dC vector
  void set_rhs_of_Hhat(double );
  
  SpMat Hhat;      // dynamical (hydrodynamic) tensor
  void set_Hhat();
  void update_Hhat();




  
  Eigen::VectorXd rhs_tmp;  // temp vector for updating dC


  SpMat dCdlambda; // matrix for storing Jacobian of dC=0 Newton solver
  SpMat mat_tmp;  // temporary matrix for storing parts of dCdlambda which change upon iteration
  void set_dCdlambda();
  void update_dCdlambda(double);



  void correct_tension(double,double,  double tol = 1e-14);


  std::vector<Atom> atoms;
  std::vector<Bond> bonds,tmpbonds;



  void compute_noise();
  void compute_tension(double);

  
  void single_inv_friction(int);

  void initial_integrate(double);
  void final_integrate(double);

  void calculate_constraint_errors(double);  


  
  int get_Nbeads() const;
  double get_temp() const;
  double get_zpara() const;
  double get_zperp() const;
  double get_bondlength() const;
  double get_timescale(double) const;
  
  Eigen::SimplicialLDLT< SpMat, Eigen::Lower > Gmunu_solver;
  Eigen::SimplicialLDLT< SpMat, Eigen::Lower > Hhat_solver;

  Eigen::SparseLU< SpMat > jacob_solver;
  Eigen::VectorXd costhetas; // costhetas[i] = u[i+2].u[i+1]

  std::vector<Eigen::VectorXd> tension_corrections;
  
protected:
  int Nbeads;           // number of polymer beads
  double bondlength;    // length of rods connecting beads
  double zpara;         // parallel coefficient of friction
  double zperp;         // perpendicular coefficient of friction
  double temp;
  double kappa;         // bare bending energy
  
  Eigen::Vector3d x0,xN;         // location of tethered bead



private:
  
  int seed;
  std::mt19937 gen;

  std::uniform_real_distribution<double> dist;

  std::vector<T> init_G_coeffsmatrix();

  
  void compute_effective_kappa();

  std::vector<T> init_Hhat_coeffsmatrix();  

  std::vector<T> init_dCdlambda_coeffsmatrix();

  
  double Hhat_diag_val(int);
  double dCdlambda_diag_val(int);
  
  double Hhat_loweroff_val(int);
  double dCdlambda_loweroff_val(int);
  double dCdlambda_upperoff_val(int);

  double Hhat_endblocks(int,int,int);
  double Hhat_leftside(int);
  double dCdlambda_leftside(int);
  double Hhat_bottomside(int);
  double dCdlambda_bottomside(int);



  Eigen::VectorXd tDets;
  Eigen::VectorXd bDets;


  Eigen::VectorXd k_effs;
  Eigen::VectorXd end_inverses;




};

#endif
