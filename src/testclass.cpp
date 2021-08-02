#include "testclass.hpp"

#include <iostream>
#include <cmath>
#include <stdexcept>


typedef Eigen::Triplet<double> Triplet;

/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
TestClass::TestClass(std::vector<Atom> as, std::vector<Bond> bs)
{
  atoms = as;
  bonds = bs;
  Nbeads = atoms.size();
  
}

/* -------------------------------------------------------------------------- */
/* Destructor */
/* -------------------------------------------------------------------------- */
TestClass::~TestClass()
{
}


template<int HHAT>
TestClass::SpMat TestClass::compute_matHat()
{
  std::vector<Triplet> coefficients;

  Matrix<Eigen::Vector3d> boldn = set_boldn_matrix();


  for (int mu = 0; mu < Nbeads-1; mu++ ){
    for (int rho = 0; rho < Nbeads-1; rho ++ ) {
      double sum = 0;
      for (int i = 0; i < Nbeads; i++) {
	if (HHAT) {
	  sum += boldn(i,mu).dot(atoms[i].friction*boldn(i,rho));
	} else {
	  sum += boldn(i,mu).dot(boldn(i,rho));
	}
	
      }

      coefficients.push_back(Triplet(mu,rho,sum));
    }
  }
  
  SpMat mat(Nbeads-1,Nbeads-1);
  
  mat.setFromTriplets(coefficients.begin(),coefficients.end());
  
  return mat;
}

TestClass::SpMat TestClass::compute_Ghat()
{
  return compute_matHat<0>();
}

TestClass::SpMat TestClass::compute_Hhat()
{
  return compute_matHat<1>();
}

Eigen::MatrixXd TestClass::invert_matrix(TestClass::SpMat& mat)
{
  solver.analyzePattern(mat);
  solver.factorize(mat);

  SpMat Identity(Nbeads-1,Nbeads-1);
  Identity.setIdentity();

  return solver.solve(Identity);

}


Matrix <Eigen::Matrix3d> TestClass::geo_projection(Eigen::MatrixXd& inverse,
						   Matrix<Eigen::Vector3d>& boldn)
{
  return projection<0>(inverse,boldn);
}

Matrix <Eigen::Matrix3d> TestClass::dyn_projection(Eigen::MatrixXd& inverse,
						   Matrix<Eigen::Vector3d>& boldn)
{
  return projection<1>(inverse,boldn);
}

template<int DYNAMIC>
Matrix<Eigen::Matrix3d> TestClass::projection(Eigen::MatrixXd& inverse,
					      Matrix<Eigen::Vector3d>& boldn)

{

  Matrix<Eigen::Matrix3d> mat(Nbeads,Nbeads);
  for (int i = 0; i < Nbeads; i++ ) {
    for (int j = 0; j < Nbeads; j++) {
      

      Eigen::Matrix3d ptemp;

      ptemp.setZero();
  
      if (i == j) {
	ptemp.setIdentity();
      }


      
      for (int mu = 0; mu < Nbeads -1 ; mu++) {
	for (int rho = 0; rho < Nbeads -1 ; rho ++) {

	  if (DYNAMIC) {
	    
	    ptemp += -inverse(mu,rho)*atoms[i].friction*boldn(i,mu)*boldn(j,rho).transpose();
	    
	  } else {
	    ptemp += -inverse(mu,rho)*boldn(i,mu)*boldn(j,rho).transpose();
	  }
      
	  
	}
      }

      mat(i,j) = ptemp;
      
    }
  }
  
  return mat;

}


Matrix<Eigen::Vector3d> TestClass::set_boldn_matrix()
{
  Matrix<Eigen::Vector3d> boldn(Nbeads, Nbeads-1);
  
  for (int i = 0; i < Nbeads; i++) {
    for (int mu = 0; mu < Nbeads-1; mu++) {
      if (i == mu) {
	boldn(i,mu) = -1*bonds[mu].rod;
      } else if (i == mu+1) {
	boldn(i,mu) = bonds[mu].rod;
      } else {
	boldn(i,mu).setZero();
      }
    }
  }
  return boldn;
      
}
