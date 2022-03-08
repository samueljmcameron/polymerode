#include "testclass.hpp"

#include <iostream>
#include <fstream>
#include <string>

#include <chrono>

int main(int argc, char* argv[])
{

  int N = 9;
  double rodin = 0.1;
  double dt = 0.1;
  double zpara = 2.0;
  double zperp = 1.0;
  double kbT = 0.0;
  double kappa = 1.0;

  int seed = 124890;
  
  Eigen::Vector3d x0(0.0,0.0,0.0);


  TestClass pmer(N,rodin,dt,zpara,zperp,kbT,kappa,x0,seed);
  pmer.compute_tangents_and_rods_and_friction();
  double norm;
  for (int mu = 0; mu < N-1; mu++) {
    norm = sqrt((pmer.atoms[mu+1].R-pmer.atoms[mu].R).dot(pmer.atoms[mu+1].R-pmer.atoms[mu].R));

    std::cout << "mu = " << mu << ", |u_mu|/a = " << norm/rodin << std::endl;
  }

  std::vector<Eigen::Vector3d> R_previous_time(N);
  std::vector<Eigen::Vector3d> bond_previous_time(N-1);

  for (int i = 0; i < N-1; i++) {
    R_previous_time[i] = pmer.atoms[i].R;
    bond_previous_time[i] = pmer.bonds[i].rod;
  }
  R_previous_time[N-1] = pmer.atoms[N-1].R;
  std::cout << dt*kbT/(rodin*zperp) << std::endl;
  

  
  pmer.set_unprojected_noise();
  pmer.set_G();
  pmer.set_Hhat();
  pmer.compute_noise();
  pmer.compute_uc_forces();


  
  pmer.compute_tension();
  pmer.initial_integrate();


  pmer.compute_tangents_and_rods_and_friction();
  
  for (int mu = 0; mu < N-1; mu++) {
    norm = sqrt((pmer.atoms[mu+1].R-pmer.atoms[mu].R).dot(pmer.atoms[mu+1].R-pmer.atoms[mu].R));
    
    std::cout << "mu = " << mu << ", |u_mu|/a = " << norm/rodin << std::endl;
  }
  
  for (int i = 0; i < N; i++) {
    std::cout << "i = " << i << ", F[i] = " << pmer.atoms[i].F << std::endl;
  }
  
  pmer.update_Hhat();
  pmer.compute_uc_forces();
  pmer.compute_tension();
  pmer.final_integrate();
  

  double tmp;

  for (int mu = 0; mu < N-1; mu++) {

    tmp = pmer.bonds[mu].rod.dot(pmer.atoms[mu+1].R-R_previous_time[mu+1]
				 -(pmer.atoms[mu].R-R_previous_time[mu]));

    std::cout << "mu = " << mu << ", dC_mu/dt = " << tmp/dt << std::endl;

  }

  pmer.compute_tangents_and_rods_and_friction();    
  for (int mu = 0; mu < N-1; mu++) {
    norm = sqrt((pmer.atoms[mu+1].R-pmer.atoms[mu].R).dot(pmer.atoms[mu+1].R-pmer.atoms[mu].R));

    std::cout << "mu = " << mu << ", |u_mu|/a = " << norm/rodin << std::endl;
  }

  return 0;
  
}

