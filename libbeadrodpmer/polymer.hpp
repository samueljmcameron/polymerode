
#ifndef BEADRODPMER_POLYMER_HPP
#define BEADRODPMER_POLYMER_HPP

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
  enum PTYPE {NONE,SINGLE,DOUBLE};

  void compute_tangents_and_friction(const Eigen::Ref<const Eigen::Matrix3Xd> &xs);
  void set_unprojected_noise(double);

  void single_inv_friction(int);

  int get_Nbeads() const;
  double get_temp() const;
  double get_zpara() const;
  double get_zperp() const;
  double get_bondlength() const;
  double get_bending() const;
  int get_seed() const;
  double get_timescale(double) const;


  virtual void setup(const Eigen::Ref<const Eigen::Matrix3Xd> &) = 0;

  virtual void first_step(Eigen::Ref<Eigen::Matrix3Xd> ,
			  Eigen::Ref<Eigen::Matrix3Xd> ,double) = 0;

  virtual int second_step(Eigen::Ref<Eigen::Matrix3Xd> ,
			  Eigen::Ref<Eigen::Matrix3Xd> ,double,int) = 0;

  
  Eigen::VectorXd constraint_errors;
  Eigen::Matrix3Xd bonds;
  
protected:
  int Nbeads;           // number of polymer beads
  double bondlength;    // length of rods connecting beads
  double zpara;         // parallel coefficient of friction
  double zperp;         // perpendicular coefficient of friction
  double temp;
  double kappa;         // bare bending energy


  Eigen::VectorXd costhetas; // costhetas[i] = u[i+2].u[i+1]

  std::mt19937 gen;

  std::uniform_real_distribution<double> dist;





  Eigen::Matrix3Xd t_forces;
  Eigen::Matrix3Xd noises,unprojected_noises;
  Eigen::Matrix3Xd tangents;

  
  std::vector<Eigen::Matrix3d> frictions;


  

private:

  int seed;
  void resize(int Nbeads) {
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
