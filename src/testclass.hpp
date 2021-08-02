
#ifndef TESTCLASS_HPP
#define TESTCLASS_HPP

#include "matrix.hpp"
#include "atom.hpp"
#include "bond.hpp"

#include <vector>



class TestClass {
public:

  typedef Eigen::SparseMatrix<double> SpMat; // declares column-major  
  // constructor
  TestClass(std::vector<Atom>,std::vector<Bond>);
  ~TestClass();  
  SpMat compute_Ghat();
  SpMat compute_Hhat();
  Eigen::MatrixXd invert_matrix(SpMat&);
  Eigen::Matrix3d projection(int, int);
  Eigen::VectorXd compute_boldn(int,int);

  Matrix <Eigen::Matrix3d> geo_projection(Eigen::MatrixXd&,
					  Matrix<Eigen::Vector3d>&);

  Matrix <Eigen::Matrix3d> dyn_projection(Eigen::MatrixXd&,
					  Matrix<Eigen::Vector3d>&);


  
  Matrix<Eigen::Vector3d> set_boldn_matrix();
  std::vector<Atom> atoms;
  std::vector<Bond> bonds;

private:
    Eigen::SimplicialLDLT< SpMat, Eigen::Lower > solver;

  int Nbeads;

  template<int>
  SpMat compute_matHat();

  template<int>
  Matrix <Eigen::Matrix3d> projection(Eigen::MatrixXd&,
				      Matrix<Eigen::Vector3d>&);



  
};
#endif
