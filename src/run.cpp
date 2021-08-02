#include "run.hpp"

#include <fstream>

void run(GlobalParams& gp, Polymer& pmer)
{

  int numsteps = gp.steps;

  double dt = gp.timestep;

  int dump_every = gp.dump_every;


  std::ofstream myfile;
  myfile.open(gp.dump_file);


  pmer.compute_tangents_and_rods_and_friction();

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
  
  
  
    
  pmer.set_unprojected_noise(dt);
  pmer.set_G();
  pmer.set_Hhat();
  pmer.compute_noise();
  pmer.compute_uc_forces();
  pmer.compute_tension();
  pmer.initial_integrate(dt);
    
    
  pmer.compute_tangents_and_rods_and_friction();
    
  pmer.update_Hhat();
  pmer.compute_uc_forces();
  pmer.compute_tension();
  pmer.final_integrate(dt);
    
    
  pmer.compute_tangents_and_rods_and_friction();
  pmer.rescale_positions(true);


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
    
    
    pmer.compute_tangents_and_rods_and_friction();    
      
    pmer.update_Hhat();
    pmer.compute_uc_forces();
    pmer.compute_tension();
    
    pmer.final_integrate(dt);
    pmer.compute_tangents_and_rods_and_friction();


    pmer.rescale_positions(true);
    
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
    

  myfile.close();
  return;
}
