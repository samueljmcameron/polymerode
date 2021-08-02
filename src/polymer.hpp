
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

  Eigen::VectorXd rhs_of_G;  // rhs_of_G vector
  Eigen::VectorXd rhs_of_Hhat;  // rhs_of_H vector
  Eigen::VectorXd dummy_for_noise, tension;



  void rescale_positions(bool);

  virtual void init_atoms();
  void compute_tangents_and_rods_and_friction();

  
  SpMat Gmunu;      // geometric tensor
  void set_G();
  void update_G();
  void set_rhs_of_G();
  
  SpMat Hhat;      // dynamical tensor
  void set_Hhat();
  void update_Hhat();
  void set_rhs_of_Hhat();
  

  std::vector<Atom> atoms;
  std::vector<Bond> bonds;

  void set_unprojected_noise(double);

  void compute_noise();
  void compute_tension();
  void compute_uc_forces();
  
  void single_inv_friction(int);

  void initial_integrate(double);
  void final_integrate(double);

  
  int get_Nbeads();
  double get_temp();
  double get_zpara();
  double get_zperp();
  double get_bondlength();
  double get_timescale(double);
  
  Eigen::SimplicialLDLT< SpMat, Eigen::Lower > solver;
  
protected:
  int Nbeads;           // number of polymer beads
  double bondlength;    // length of rods connecting beads
  double zpara;         // parallel coefficient of friction
  double zperp;         // perpendicular coefficient of friction
  double temp;
  double kappa;         // bare bending energy
  
  Eigen::Vector3d x0;         // location of tethered bead



private:
  
  int seed;
  std::mt19937 gen;

  std::uniform_real_distribution<double> dist;



  
  std::vector<T> init_G_coeffsmatrix();

  double Hhat_diag_val(int);
  double Hhat_loweroff_val(int);
  std::vector<T> init_Hhat_coeffsmatrix();

  Eigen::VectorXd tDets;
  Eigen::VectorXd bDets;

  Eigen::VectorXd costhetas; // costhetas[i] = u[i+1].u[i]
  Eigen::VectorXd k_effs;

  void compute_effective_kappa();


};

#endif
