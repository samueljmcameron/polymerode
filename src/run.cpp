#include "run.hpp"

#include <iostream>
#include <fstream>

void run(GlobalParams& gp, Polymer& pmer)
{

  int numsteps = gp.steps;

  double dt = gp.timestep;

  int dump_every = gp.dump_every;


  std::ofstream myfile;
  myfile.open(gp.dump_file);


  pmer.compute_tangents_and_rods_and_friction();

  Eigen::Vector3d startpoint(pmer.atoms[0].R);
  Eigen::Vector3d endpoint(pmer.atoms[pmer.get_Nbeads()-1].R);
  
  // save at time t = 0 always.
  
  if (myfile.is_open()) {
    
    myfile << "TIMESTEP: " << 0 << std::endl;
    int N = pmer.get_Nbeads();
    Eigen::MatrixXd msave(N,6);
    
    Eigen::Vector3d dummy(0,0,0);
    
    for (int i = 0 ; i < N-1; i++) {
      
      msave.row(i) << pmer.atoms[i].R.transpose(), pmer.bonds[i].rod.transpose();
    }
    msave.row(N-1) << pmer.atoms[N-1].R.transpose(), dummy.transpose();
    
    myfile << msave << std::endl;
  }
  int max_iters = 5;
  
  
    
  pmer.set_unprojected_noise(dt);
  pmer.set_G();

  
  pmer.set_Hhat();

  
  pmer.compute_noise();
  pmer.compute_uc_forces();
  pmer.compute_tension();
  pmer.initial_integrate(dt);

  for (int iters = 0; iters < max_iters; iters ++ ) {
  
    pmer.compute_tangents_and_rods_and_friction();

    pmer.update_Hhat();
    pmer.compute_uc_forces();
    pmer.compute_tension();
    pmer.final_integrate(dt);
    
    
  }
  
  pmer.compute_tangents_and_rods_and_friction();

  // save first step only if saving on every step.
  if (myfile.is_open() && dump_every == 1) {
    
    myfile << "TIMESTEP: " << 1 << std::endl;
    int N = pmer.get_Nbeads();
    Eigen::MatrixXd msave(N,6);
    
    Eigen::Vector3d dummy(0,0,0);
    
    for (int i = 0 ; i < N-1; i++) {
      
      msave.row(i) << pmer.atoms[i].R.transpose(), pmer.bonds[i].rod.transpose();
    }
    msave.row(N-1) << pmer.atoms[N-1].R.transpose(), dummy.transpose();
    
    myfile << msave << std::endl;
  }
  
  // continue the remaining numstep - 2 steps
  for (int i = 2; i <= numsteps; i++) {
    

    pmer.set_unprojected_noise(dt);
    pmer.update_G();
    pmer.update_Hhat();
    pmer.compute_noise();
    pmer.compute_uc_forces();
    pmer.compute_tension();
    pmer.initial_integrate(dt);
    



    for (int iters = 0; iters < max_iters; iters ++ ) {
      
      pmer.compute_tangents_and_rods_and_friction();
      
      pmer.update_Hhat();
      pmer.compute_uc_forces();
      pmer.compute_tension();
      pmer.final_integrate(dt);
      
      
    }
    pmer.compute_tangents_and_rods_and_friction();  
    
    if (myfile.is_open() && i % dump_every == 0) {

      myfile << "TIMESTEP: " << i << std::endl;
      int N = pmer.get_Nbeads();
      Eigen::MatrixXd msave(N,6);
      
      Eigen::Vector3d dummy(0,0,0);

      for (int i = 0 ; i < N-1; i++) {

	msave.row(i) << pmer.atoms[i].R.transpose(), pmer.bonds[i].rod.transpose();
      }
      msave.row(N-1) << pmer.atoms[N-1].R.transpose(), dummy.transpose();
      
      myfile << msave << std::endl;
    }

      
  }

  
  double norm;
  
  std::cout << "difference in start position = " << pmer.atoms[0].R-startpoint << std::endl;
  for (int mu = 0; mu < pmer.get_Nbeads()-1; mu++) {
    norm = sqrt((pmer.atoms[mu+1].R-pmer.atoms[mu].R).dot(pmer.atoms[mu+1].R-pmer.atoms[mu].R));
    
    std::cout << "mu = " << mu << ", |u_mu|/a = " << norm/pmer.get_bondlength() << std::endl;
  }
  
  std::cout << "difference in end position = "
	    << pmer.atoms[pmer.get_Nbeads()-1].R-endpoint << std::endl;
  
  myfile.close();
  return;
}
