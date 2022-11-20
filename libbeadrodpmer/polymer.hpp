
#ifndef BEADRODPMER_POLYMER_HPP
#define BEADRODPMER_POLYMER_HPP

#include "atom.hpp"
#include "bond.hpp"


#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <vector>
#include <string>
#include <random>
#include <functional>

namespace BeadRodPmer {
typedef Eigen::SparseMatrix<double> SpMat; // declares column-major
typedef Eigen::Triplet<double> T;


  

class Polymer {
public:
  // constructor
  Polymer(const std::vector<std::string> &);
  enum PTYPE {NONE,SINGLE,DOUBLE};

  void compute_tangents_and_friction();
  void set_unprojected_noise(double);
  void compute_uc_forces();

  void single_inv_friction(int);

  void add_external_force(const std::vector<double> &,int);

  int get_Nbeads() const;
  double get_temp() const;
  double get_zpara() const;
  double get_zperp() const;
  double get_bondlength() const;
  double get_timescale(double) const;

  Eigen::Vector3d get_x0() const;
  Eigen::Vector3d get_xN() const;

  std::vector<Atom> atoms;
  std::vector<Bond> bonds;
  
  std::vector<Eigen::Vector3d> Rtmp;    // temporary storate for midstep

  
  std::vector<int> nuc_beads;
  std::vector<double> nuc_strengths;
  std::vector<double> nuc_maxs;
  std::vector<double> nuc_widths;

  SpMat Gmunu;      // geometric tensor
  Eigen::VectorXd rhs_of_G;  // rhs_of_G vector
  Eigen::VectorXd dummy_for_noise;


  SpMat Hhat;      // dynamical (hydrodynamic) tensor  
  Eigen::VectorXd rhs_of_Hhat;  // rhs_of_H vector AND dC vector
  Eigen::VectorXd tension;


  
  SpMat dCdlambda; // matrix for storing Jacobian of dC=0 Newton solver  
  Eigen::VectorXd constraint_errors;
  Eigen::VectorXd negative_tension_change;


  Eigen::SimplicialLDLT< SpMat, Eigen::Lower > Gmunu_solver;
  Eigen::SimplicialLDLT< SpMat, Eigen::Lower > Hhat_solver;

  Eigen::SparseLU< SpMat > jacob_solver;
  Eigen::VectorXd costhetas; // costhetas[i] = u[i+2].u[i+1]


  int equilibration_steps;
  double initspringK;   // spring constant for initializing double tether
  double initdt;        // time step for initializing double tether
  double inittolerance; // how far from specified tether is acceptable when




  virtual int single_step(double, double,
			  const std::vector<std::vector<double>> &,
			  int, int , bool ) = 0;


    
  virtual void set_G() = 0;

  virtual void set_Hhat() = 0;

  virtual void set_dCdlambda() = 0;


protected:
  int Nbeads;           // number of polymer beads
  double bondlength;    // length of rods connecting beads
  double zpara;         // parallel coefficient of friction
  double zperp;         // perpendicular coefficient of friction
  double temp;
  double kappa;         // bare bending energy

                                  //  initializing double tether
  
  Eigen::Vector3d x0,xN;         // desired location of tethered beads

  bool flag_x0, flag_xN, flag_initdoubleteth;

  std::mt19937 gen;

  std::uniform_real_distribution<double> dist;
  
  std::vector<Bond> tmpbonds;




  virtual void compute_noise() = 0;
  


  virtual void compute_effective_kappa() = 0;
  
  double Hhat_endblocks(int,int,int);
  double Hhat_leftside(int);
  double dCdlambda_leftside(int);
  double Hhat_bottomside(int);
  double dCdlambda_bottomside(int);

  void rescale_positions(bool);

  Eigen::VectorXd tDets;
  Eigen::VectorXd bDets;


  Eigen::VectorXd k_effs;
  Eigen::VectorXd end_inverses;

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

private:

  double Hhat_diag_val(int);
  double Hhat_loweroff_val(int);
  
  double dCdlambda_diag_val(int);
  double dCdlambda_loweroff_val(int);
  double dCdlambda_upperoff_val(int);

};
};
#endif
