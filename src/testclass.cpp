#include "testclass.hpp"

#include <iostream>
#include <cmath>
#include <stdexcept>


/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
TestClass::TestClass(int Nin, double rodin, double Delta_tin, double zparain,
		     double zperpin, double tempin,double kappain,
		     Eigen::Vector3d x0in,  int seedin)
  : Base(Nin,rodin, Delta_tin,zparain,zperpin, tempin,kappain, x0in ,seedin)

{

  
}

/* -------------------------------------------------------------------------- */
/* Destructor */
/* -------------------------------------------------------------------------- */
TestClass::~TestClass()
{
}


Eigen::Matrix3d TestClass::projection(int i, int j)
{

  
  solver.factorize(Gmunu);

  SpMat Identity(Nbeads-1,Nbeads-1);
  Identity.setIdentity();

  Eigen::MatrixXd Ginv = solver.solve(Identity);


  Eigen::VectorXd nx_i = compute_boldn_x(i);
  Eigen::VectorXd nx_j = compute_boldn_x(j);

  Eigen::VectorXd ny_i = compute_boldn_y(i);
  Eigen::VectorXd ny_j = compute_boldn_y(j);

  Eigen::VectorXd nz_i = compute_boldn_z(i);
  Eigen::VectorXd nz_j = compute_boldn_z(j);
  

  Eigen::Matrix3d ptemp;

  ptemp.setZero();

  if (i == j) {
    ptemp(0,0) = ptemp(1,1) = ptemp(2,2) = 1.0;
  }
  
  ptemp(0,0) -= nx_i.dot(Ginv*nx_j);
  ptemp(0,1) -= nx_i.dot(Ginv*ny_j);
  ptemp(0,2) -= nx_i.dot(Ginv*nz_j);

  ptemp(1,0) -= ny_i.dot(Ginv*nx_j);
  ptemp(1,1) -= ny_i.dot(Ginv*ny_j);
  ptemp(1,2) -= ny_i.dot(Ginv*nz_j);

  ptemp(2,0) -= nz_i.dot(Ginv*nx_j);
  ptemp(2,1) -= nz_i.dot(Ginv*ny_j);
  ptemp(2,2) -= nz_i.dot(Ginv*nz_j);
  
  
  return ptemp;
  

}

Eigen::VectorXd TestClass::compute_boldn_x(int i)
{


  Eigen::VectorXd bmn(Nbeads-1);

  bmn.setZero();


  if (i > 0) {
    bmn(i-1) = u_x(i-1);
  }
  if (i < Nbeads-1) {
    bmn(i) = -u_x(i);
  }

  return bmn;
  

}



Eigen::VectorXd TestClass::compute_boldn_y(int i)
{


  Eigen::VectorXd bmn(Nbeads-1);

  bmn.setZero();


  if (i > 0) {
    bmn(i-1) = u_y(i-1);
  }
  if (i < Nbeads-1) {
    bmn(i) = -u_y(i);
  }

  return bmn;
  

}

Eigen::VectorXd TestClass::compute_boldn_z(int i)
{


  Eigen::VectorXd bmn(Nbeads-1);

  bmn.setZero();


  if (i > 0) {
    bmn(i-1) = u_z(i-1);
  }
  if (i < Nbeads-1) {
    bmn(i) = -u_z(i);
  }

  return bmn;
  

}

