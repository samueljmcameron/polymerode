#include "run.hpp"
#include "iovtk.hpp"

#include <iostream>


struct atEnd
{
  Eigen::Vector3d X_of_t;
  Eigen::Vector3d dXdt;

};

void update_bead1(atEnd & bead1,const Eigen::Vector3d &x0,double t)
{

  double omega = 100;
  double rad = 0.1;
  bead1.X_of_t(0) = x0(0) + rad*(cos(omega*t)-1);
  bead1.X_of_t(1) = x0(1) + rad*sin(omega*t);
  bead1.X_of_t(2) = x0(2);

  bead1.dXdt(0) = -rad*omega*sin(omega*t);
  bead1.dXdt(1) = rad*omega*cos(omega*t);
  bead1.dXdt(2) = 0;

  return;
  
}

void update_beadN(atEnd & beadN,const Eigen::Vector3d &xN,double t)
{
  beadN.X_of_t = xN;
  beadN.dXdt.setZero();
}



void single_step(double & t,BeadRodPmer::DoubleTether & pmer,
		 atEnd & bead1,atEnd & beadN,const Eigen::Vector3d &x0,
		 const Eigen::Vector3d &xN,double dt)
{
  pmer.set_unprojected_noise(dt);
  pmer.update_G();
  pmer.update_Hhat();
  pmer.compute_noise();
  pmer.compute_effective_kappa();
  pmer.compute_uc_forces();
  
  update_bead1(bead1,x0,t+dt/2);
  update_beadN(beadN,xN,t+dt/2);

  
  pmer.compute_tension(bead1.dXdt,beadN.dXdt);
  pmer.initial_integrate(dt);
  

  
  pmer.update_Hhat();
  pmer.compute_effective_kappa();
  pmer.compute_uc_forces();
  
  update_bead1(bead1,x0,t+dt);
  update_beadN(beadN,xN,t+dt);
  
  pmer.compute_tension(bead1.dXdt,beadN.dXdt);
  pmer.correct_tension(dt,bead1.X_of_t,beadN.X_of_t,1e-8);
  pmer.final_integrate(dt);
  
  pmer.compute_tangents_and_friction();
  t += dt;
  return;
}




void run(BeadRodPmer::GlobalParams& gp, BeadRodPmer::DoubleTether& pmer)
{

  int numsteps = gp.steps;

  double dt = gp.timestep;

  int dump_every = gp.dump_every;

  double t = 0;

  Eigen::Vector3d startpoint(pmer.atoms[0].R);
  Eigen::Vector3d endpoint(pmer.atoms[pmer.get_Nbeads()-1].R);


  atEnd bead1,beadN;
  
  // save at time t = 0 always.

  std::string collection_name = gp.dump_file + std::string(".pvd");
  BeadRodPmer::ioVTK::writeVTKcollectionHeader(collection_name);

  
  std::string fname = gp.dump_file + std::string("_") + std::to_string(0) + std::string(".vtp");
  BeadRodPmer::ioVTK::writeVTKPolyData(fname,pmer);
  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  pmer.compute_tangents_and_friction();
  
  pmer.set_Hhat();
  pmer.set_dCdlambda();
  pmer.set_G();
  
  single_step(t,pmer,bead1,beadN,startpoint,endpoint,dt);

  if (dump_every == 1) {

    fname = gp.dump_file + std::string("_") + std::to_string(1) + std::string(".vtp");
  
    BeadRodPmer::ioVTK::writeVTKPolyData(fname,pmer);

    BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  }

  // continue the remaining numstep - 2 steps
  for (int i = 2; i <= numsteps; i++) {
    std::cout << "t = " << t << std::endl;

    single_step(t,pmer,bead1,beadN,startpoint,endpoint,dt);
    
    if (i % dump_every == 0) {
      fname = gp.dump_file + std::string("_") + std::to_string(i) + std::string(".vtp");
  
      BeadRodPmer::ioVTK::writeVTKPolyData(fname,pmer);

      BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);
      
    }
      
  }


  BeadRodPmer::ioVTK::writeVTKcollectionFooter(collection_name);

  std::cout << "Constraints are different than zero: " << std::endl;
  std::cout << pmer.constraint_errors << std::endl;
  

  return;
}


void single_step(double & t,BeadRodPmer::SingleTether & pmer,
		 atEnd & bead1, const Eigen::Vector3d &x0,double dt)
{
  pmer.set_unprojected_noise(dt);
  pmer.update_G();
  pmer.update_Hhat();
  pmer.compute_noise();
  pmer.compute_effective_kappa();
  pmer.compute_uc_forces();
  
  update_bead1(bead1,x0,t+dt/2);
  
  pmer.compute_tension(bead1.dXdt);
  pmer.initial_integrate(dt);
  

  
  pmer.update_Hhat();
  pmer.compute_effective_kappa();
  pmer.compute_uc_forces();
  
  update_bead1(bead1,x0,t+dt);

  
  pmer.compute_tension(bead1.dXdt);
  pmer.correct_tension(dt,bead1.X_of_t,1e-8);
  pmer.final_integrate(dt);
  
  pmer.compute_tangents_and_friction();
  t += dt;
  return;
}


void run(BeadRodPmer::GlobalParams& gp, BeadRodPmer::SingleTether& pmer)
{

  int numsteps = gp.steps;

  double dt = gp.timestep;

  int dump_every = gp.dump_every;

  double t = 0;

  Eigen::Vector3d startpoint(pmer.atoms[0].R);


  atEnd bead1;
  
  // save at time t = 0 always.

  std::string collection_name = gp.dump_file + std::string(".pvd");
  BeadRodPmer::ioVTK::writeVTKcollectionHeader(collection_name);

  
  std::string fname = gp.dump_file + std::string("_") + std::to_string(0) + std::string(".vtp");
  BeadRodPmer::ioVTK::writeVTKPolyData(fname,pmer);
  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  pmer.compute_tangents_and_friction();
  
  pmer.set_Hhat();
  pmer.set_dCdlambda();
  pmer.set_G();
  single_step(t,pmer,bead1,startpoint,dt);  

  if (dump_every == 1) {

    fname = gp.dump_file + std::string("_") + std::to_string(1) + std::string(".vtp");
  
    BeadRodPmer::ioVTK::writeVTKPolyData(fname,pmer);

    BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  }

  // continue the remaining numstep - 2 steps
  for (int i = 2; i <= numsteps; i++) {
    std::cout << "t = " << t << std::endl;
    
    single_step(t,pmer,bead1,startpoint,dt);  

    if (i % dump_every == 0) {
      fname = gp.dump_file + std::string("_") + std::to_string(i) + std::string(".vtp");
  
      BeadRodPmer::ioVTK::writeVTKPolyData(fname,pmer);

      BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);
      
    }
      
  }


  BeadRodPmer::ioVTK::writeVTKcollectionFooter(collection_name);

  std::cout << "Constraints are different than zero: " << std::endl;
  std::cout << pmer.constraint_errors << std::endl;
  

  return;
}


void single_step(double & t, BeadRodPmer::NoTether & pmer,double dt)
{
  pmer.set_unprojected_noise(dt);
  pmer.update_G();
  pmer.update_Hhat();
  pmer.compute_noise();
  pmer.compute_effective_kappa();
  pmer.compute_uc_forces();
  

  
  pmer.compute_tension();
  pmer.initial_integrate(dt);
  

  
  pmer.update_Hhat();
  pmer.compute_effective_kappa();
  pmer.compute_uc_forces();
 

  
  pmer.compute_tension();
  pmer.correct_tension(dt,1e-8);
  pmer.final_integrate(dt);
  
  pmer.compute_tangents_and_friction();
  t += dt;
  return;
}


void run(BeadRodPmer::GlobalParams& gp, BeadRodPmer::NoTether& pmer)
{

  int numsteps = gp.steps;

  double dt = gp.timestep;

  int dump_every = gp.dump_every;

  double t = 0;


  atEnd bead1;
  
  // save at time t = 0 always.

  std::string collection_name = gp.dump_file + std::string(".pvd");
  BeadRodPmer::ioVTK::writeVTKcollectionHeader(collection_name);

  
  std::string fname = gp.dump_file + std::string("_") + std::to_string(0) + std::string(".vtp");
  BeadRodPmer::ioVTK::writeVTKPolyData(fname,pmer);
  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  pmer.compute_tangents_and_friction();
  
  pmer.set_Hhat();
  pmer.set_dCdlambda();
  pmer.set_G();
  single_step(t,pmer,dt);  

  if (dump_every == 1) {

    fname = gp.dump_file + std::string("_") + std::to_string(1) + std::string(".vtp");
  
    BeadRodPmer::ioVTK::writeVTKPolyData(fname,pmer);

    BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  }

  // continue the remaining numstep - 2 steps
  for (int i = 2; i <= numsteps; i++) {
    std::cout << "t = " << t << std::endl;
    
    single_step(t,pmer,dt);  

    if (i % dump_every == 0) {
      fname = gp.dump_file + std::string("_") + std::to_string(i) + std::string(".vtp");
  
      BeadRodPmer::ioVTK::writeVTKPolyData(fname,pmer);

      BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);
      
    }
      
  }


  BeadRodPmer::ioVTK::writeVTKcollectionFooter(collection_name);

  std::cout << "Constraints are different than zero: " << std::endl;
  std::cout << pmer.constraint_errors << std::endl;
  

  return;
}
