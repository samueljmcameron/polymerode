
#include "testclass.hpp"
#include <Eigen/SparseLU>

#include <iostream>
#include <fstream>
#include <string>

#include <chrono>

int main(int argc, char* argv[])
{

  int N = 9;
  double rodin = 0.1;
  double dt = 0.00001;
  double zpara = 2.0;
  double zperp = 1.0;
  double kbT = 1.0;
  double kappa = 0.1;

  int seed = 124890;
  
  Eigen::Vector3d x0(0.0,0.0,0.0);


  TestClass pmer(N,rodin,dt,zpara,zperp,kbT,kappa,x0,seed);

  double norm;
  for (int mu = 0; mu < N-1; mu++) {
    norm = sqrt((pmer.R_x(mu+1)-pmer.R_x(mu))*(pmer.R_x(mu+1)-pmer.R_x(mu))
		+(pmer.R_y(mu+1)-pmer.R_y(mu))*(pmer.R_y(mu+1)-pmer.R_y(mu))
		+(pmer.R_z(mu+1)-pmer.R_z(mu))*(pmer.R_z(mu+1)-pmer.R_z(mu)));

    std::cout << "mu = " << mu << ", |u_mu|/a = " << norm/rodin << std::endl;
  }


  Eigen::VectorXd R_x = pmer.R_x;
  Eigen::VectorXd R_y = pmer.R_y;
  Eigen::VectorXd R_z = pmer.R_z;
  
  std::cout << dt*kbT/(rodin*zperp) << std::endl;
  
  pmer.compute_tangents();
  
  pmer.set_zetainv();
  pmer.set_unprojected_noise();
  pmer.set_G();
  pmer.set_Hhat();
  pmer.compute_noise();
  pmer.compute_uc_forces();
  pmer.compute_tension();
  pmer.initial_integrate();


  pmer.compute_tangents();


  for (int mu = 0; mu < N-1; mu++) {
    norm = sqrt((pmer.R_x(mu+1)-pmer.R_x(mu))*(pmer.R_x(mu+1)-pmer.R_x(mu))
		+(pmer.R_y(mu+1)-pmer.R_y(mu))*(pmer.R_y(mu+1)-pmer.R_y(mu))
		+(pmer.R_z(mu+1)-pmer.R_z(mu))*(pmer.R_z(mu+1)-pmer.R_z(mu)));

    std::cout << "mu = " << mu << ", |u_mu|/a = " << norm/rodin << std::endl;
  }

  
  pmer.set_zetainv();

  pmer.update_Hhat();
  pmer.compute_uc_forces();
  pmer.compute_tension();
  pmer.final_integrate();
  


  Eigen::VectorXd Rdot_x = pmer.R_x-R_x;
  Eigen::VectorXd Rdot_y = pmer.R_y-R_y;
  Eigen::VectorXd Rdot_z = pmer.R_z-R_z;

  double tmp;
  for (int mu = 0; mu < N-1; mu++) {

    tmp = (pmer.u_x(mu)*(Rdot_x(mu+1)-Rdot_x(mu))
	   +pmer.u_y(mu)*(Rdot_y(mu+1)-Rdot_y(mu))
	   +pmer.u_z(mu)*(Rdot_z(mu+1)-Rdot_z(mu)));

    std::cout << "mu = " << mu << ", dC_mu/dt = " << tmp/dt << std::endl;
  }
  pmer.compute_tangents();

  for (int mu = 0; mu < N-1; mu++) {
    norm = sqrt((pmer.R_x(mu+1)-pmer.R_x(mu))*(pmer.R_x(mu+1)-pmer.R_x(mu))
		+(pmer.R_y(mu+1)-pmer.R_y(mu))*(pmer.R_y(mu+1)-pmer.R_y(mu))
		+(pmer.R_z(mu+1)-pmer.R_z(mu))*(pmer.R_z(mu+1)-pmer.R_z(mu)));

    std::cout << "mu = " << mu << ", |u_mu|/a = " << norm/rodin << std::endl;
  }
  return 0;
  
}

