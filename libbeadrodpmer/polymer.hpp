
#ifndef BEADRODPMER_POLYMER_HPP
#define BEADRODPMER_POLYMER_HPP

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

  void single_inv_friction(int);

  void add_external_force(const Eigen::Vector3d &,int);

  int get_Nbeads() const;
  double get_temp() const;
  double get_zpara() const;
  double get_zperp() const;
  double get_bondlength() const;
  double get_timescale(double) const;

  Eigen::Vector3d get_x0() const;
  Eigen::Vector3d get_xN() const;

  std::vector<int> nuc_beads;
  std::vector<double> nuc_strengths;
  std::vector<double> nuc_maxs;
  std::vector<double> nuc_widths;



  int equilibration_steps;
  double initspringK;   // spring constant for initializing double tether
  double initdt;        // time step for initializing double tether
  double inittolerance; // how far from specified tether is acceptable when


  Eigen::Matrix3Xd xs;

  virtual int single_step(double, double,
			  const std::vector<Eigen::Vector3d> &,
			  int, int , bool ) = 0;



  virtual void setup() = 0;


  Eigen::VectorXd constraint_errors;
  
  void init_atoms_rand();
  void init_atoms_line();
  void init_atoms_caret();
  
protected:
  int Nbeads;           // number of polymer beads
  double bondlength;    // length of rods connecting beads
  double zpara;         // parallel coefficient of friction
  double zperp;         // perpendicular coefficient of friction
  double temp;
  double kappa;         // bare bending energy


  Eigen::VectorXd costhetas; // costhetas[i] = u[i+2].u[i+1]

  

                                  //  initializing double tether
  
  Eigen::Vector3d x0,xN;         // desired location of tethered beads

  bool flag_x0, flag_xN, flag_initdoubleteth;

  std::mt19937 gen;

  std::uniform_real_distribution<double> dist;





  Eigen::Matrix3Xd Fpots,t_forces;
  Eigen::Matrix3Xd noises,unprojected_noises;
  Eigen::Matrix3Xd tangents;

  Eigen::Matrix3Xd bonds;
  
  std::vector<Eigen::Matrix3d> frictions;


  

private:

  void resize(int Nbeads) {
    xs.resize(Eigen::NoChange,Nbeads);
    Fpots.resize(Eigen::NoChange,Nbeads);
    t_forces.resize(Eigen::NoChange,Nbeads);
    noises.resize(Eigen::NoChange,Nbeads);
    unprojected_noises.resize(Eigen::NoChange,Nbeads);
    tangents.resize(Eigen::NoChange,Nbeads);
    bonds.resize(Eigen::NoChange,Nbeads-1);
    frictions.resize(Nbeads);
  }

  


  
};
};
#endif
