
#ifndef BASE_HPP
#define BASE_HPP

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <vector>
#include <random>

typedef Eigen::SparseMatrix<double> SpMat; // declares column-major
typedef Eigen::Triplet<double> T;

class Base {
public:
  // constructor
  Base(int,double,double,double,double,double,double,
       Eigen::Vector3d,int);
  ~Base();  
  
  Eigen::VectorXd R_x,R_y,R_z;    // components of polymer beads
  Eigen::VectorXd Rtmp_x,Rtmp_y,Rtmp_z;    // components of polymer beads


  
  Eigen::VectorXd u_x,u_y,u_z;    // components of bond tangents
  Eigen::VectorXd tu_x,tu_y,tu_z; // components of bead tangents
  Eigen::VectorXd F_x,F_y,F_z;    // components of polymer bead forces

  
  Eigen::VectorXd eta_x,eta_y,eta_z; // components of noise
  Eigen::VectorXd unproj_eta_x,unproj_eta_y,unproj_eta_z; // unprojected noise.
  Eigen::VectorXd rhs_of_G;  // rhs_of_G vector
  Eigen::VectorXd rhs_of_Hhat;  // rhs_of_H vector
  

  Eigen::VectorXd dummy_for_noise, tension;



  void rescale_positions(bool);

  virtual void init_R();
  void compute_tangents();

  
  SpMat Gmunu;      // geometric tensor
  void set_G();
  void update_G();
  void set_rhs_of_G();
  
  SpMat Hhat;      // dynamical tensor
  void set_Hhat();
  void update_Hhat();
  void set_rhs_of_Hhat();
  


  void set_unprojected_noise();

  void compute_noise();
  void compute_tension();
  void compute_Rdot(double);
  void compute_uc_forces();
  
  std::vector<Eigen::Matrix3d> zetainvs;  
  void set_zetainv();

  void initial_integrate();
  void final_integrate();
  
  Eigen::SimplicialLDLT< SpMat, Eigen::Lower > solver;
  
protected:
  const int Nbeads;           // number of polymer beads
  const double rodlength;     // length of rods connecting beads
  const double Delta_t;       // time discretisation
  const double zpara;         // parallel coefficient of friction
  const double zperp;         // perpendicular coefficient of friction
  const double kappa;         // bare bending energy
  
  Eigen::Vector3d x0;         // location of tethered bead

  double temp;
  int seed;

private:

  std::mt19937 gen;
  std::uniform_real_distribution<double> dist;

  
  std::vector<T> init_G_coeffsmatrix();

  
  Eigen::Matrix3d zetainv_i;
  void set_zetainv_i(int i);


  double Force_mid(Eigen::VectorXd &,int );
  double Force_0(Eigen::VectorXd &);
  double Force_1(Eigen::VectorXd &);
  double Force_Nm2(Eigen::VectorXd&);
  double Force_Nm1(Eigen::VectorXd&);

  double Hhat_diag_val(int);
  double Hhat_loweroff_val(int);
  std::vector<T> init_Hhat_coeffsmatrix();

};

#endif
