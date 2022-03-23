#include "run.hpp"
#include "iovtk.hpp"

#include <iostream>


void run(GlobalParams& gp, Polymer& pmer)
{

  int numsteps = gp.steps;

  double dt = gp.timestep;

  int dump_every = gp.dump_every;

  double t = 0;

  pmer.compute_tangents_and_rods_and_friction();

  Eigen::Vector3d startpoint(pmer.atoms[0].R);
  Eigen::Vector3d endpoint(pmer.atoms[pmer.get_Nbeads()-1].R);
  
  // save at time t = 0 always.

  std::string collection_name = gp.dump_file + std::string(".pvd");
  ioVTK::writeVTKcollectionHeader(collection_name);

  
  std::string fname = gp.dump_file + std::string("_") + std::to_string(0) + std::string(".vtp");
  ioVTK::writeVTKPolyData(fname,pmer);
  ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  
  int max_iters = 0;


  
  pmer.set_unprojected_noise(dt);
  pmer.set_G();

  
  pmer.set_Hhat();

  
  pmer.compute_noise();
  pmer.compute_uc_forces();
  pmer.compute_tension(t+dt/2);
  pmer.initial_integrate(dt);
  
  pmer.compute_tangents_and_rods_and_friction();
  
  pmer.update_Hhat();
  pmer.compute_uc_forces();
  pmer.compute_tension(t+dt);
  pmer.final_integrate(dt);
  
  for (int iters = 0; iters < max_iters; iters ++ ) {
    
    pmer.compute_tangents_and_rods_and_friction();
    pmer.update_Hhat();
    pmer.compute_uc_forces();
    pmer.correct_tension();
    pmer.final_integrate(dt);    
    
  }
  
  pmer.compute_tangents_and_rods_and_friction();
  t += dt;

  if (dump_every == 1) {


    fname = gp.dump_file + std::string("_") + std::to_string(1) + std::string(".vtp");
  
    ioVTK::writeVTKPolyData(fname,pmer);

    ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  }
  // continue the remaining numstep - 2 steps
  for (int i = 2; i <= numsteps; i++) {
    

    pmer.set_unprojected_noise(dt);
    pmer.update_G();
    pmer.update_Hhat();
    pmer.compute_noise();
    pmer.compute_uc_forces();
    pmer.compute_tension(t+dt/2);
    pmer.initial_integrate(dt);
    
    pmer.compute_tangents_and_rods_and_friction();
    
    pmer.update_Hhat();
    pmer.compute_uc_forces();
    pmer.compute_tension(t+dt);
    pmer.final_integrate(dt);

    for (int iters = 0; iters < max_iters; iters ++ ) {

      pmer.compute_tangents_and_rods_and_friction();
      pmer.update_Hhat();
      pmer.compute_uc_forces();
      pmer.correct_tension();
      pmer.final_integrate(dt);    
      
    }
    pmer.compute_tangents_and_rods_and_friction();  
    t += dt;

    if (i % dump_every == 0) {
      fname = gp.dump_file + std::string("_") + std::to_string(i) + std::string(".vtp");
  
      ioVTK::writeVTKPolyData(fname,pmer);

      ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);
    }
      
  }


  ioVTK::writeVTKcollectionFooter(collection_name);
  
  double norm;

  double omega = 100.0;
  double rad = 0.1;

  startpoint(0) += rad*(cos(omega*t)-1);
  startpoint(1) += rad*sin(omega*t);

  
  std::cout << "difference in start position = " << pmer.atoms[0].R-startpoint << std::endl;
  for (int mu = 0; mu < pmer.get_Nbeads()-1; mu++) {
    norm = (pmer.atoms[mu+1].R-pmer.atoms[mu].R).norm();
    
    std::cout << "mu = " << mu << ", |u_mu|/a = " << norm/pmer.get_bondlength() << std::endl;
  }
  
  std::cout << "difference in end position = "
	    << pmer.atoms[pmer.get_Nbeads()-1].R-endpoint << std::endl;
  

  return;
}
