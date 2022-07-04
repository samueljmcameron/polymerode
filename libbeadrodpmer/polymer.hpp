
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

namespace BeadRodPmer {
typedef Eigen::SparseMatrix<double> SpMat; // declares column-major
typedef Eigen::Triplet<double> T;


  

class Polymer {
public:
  // constructor
  Polymer(const std::vector<std::string> &);
  ~Polymer();  


  virtual void init_atoms();
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

  std::vector<Atom> atoms;
  std::vector<Bond> bonds;
  
  std::vector<Eigen::Vector3d> Rtmp;    // temporary storate for midstep

  
  std::vector<int> nuc_beads;


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

  
protected:
  int Nbeads;           // number of polymer beads
  double bondlength;    // length of rods connecting beads
  double zpara;         // parallel coefficient of friction
  double zperp;         // perpendicular coefficient of friction
  double temp;
  double kappa;         // bare bending energy
  
  Eigen::Vector3d x0,xN;         // location of tethered bead
  

  std::mt19937 gen;

  std::uniform_real_distribution<double> dist;
  
  std::vector<Bond> tmpbonds;
  
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

  void rescale_positions(bool);

  Eigen::VectorXd tDets;
  Eigen::VectorXd bDets;


  Eigen::VectorXd k_effs;
  Eigen::VectorXd end_inverses;


};
};
#endif
