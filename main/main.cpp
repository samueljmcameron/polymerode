
#include "PolyConfig.h"
#include "testclass.hpp"
#include <Eigen/SparseLU>

#include <iostream>
#include <fstream>
#include <string>

#include <chrono>

int main(int argc, char* argv[])
{
  std::cout << argv[0] << " version " << poly_VERSION_MAJOR << "."
	    << poly_VERSION_MINOR << std::endl;

  int nummols = 40;
  int numsteps = 1000000;

  
  int N = 9;
  double rodin = 0.1;
  double dt = 0.0001;
  double zpara = 1.0;
  double zperp = 2.0;
  double kbT = 1.0;
  double kappa = 0.0;
  int seed;
  Eigen::Vector3d x0(0.0,0.0,0.0);


  std::ofstream myfiles[nummols];

  
  for (int molnum = 0; molnum < nummols; molnum ++ ) {

    std::ostringstream filename;

    filename << "mol" << molnum << ".txt";

    seed = std::rand();
  



    TestClass pmer(N,rodin,dt,zpara,zperp,kbT,kappa,x0,seed);


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
    
    pmer.set_zetainv();
    
    pmer.update_Hhat();
    pmer.compute_uc_forces();
    pmer.compute_tension();
    pmer.final_integrate();
    
    
    pmer.compute_tangents();
    pmer.rescale_positions(true);
    
    for (int i = 0; i < numsteps; i++) {
  
      
      pmer.set_zetainv();
      pmer.set_unprojected_noise();
      pmer.update_G();
      pmer.update_Hhat();
      pmer.compute_noise();
      pmer.compute_uc_forces();
      pmer.compute_tension();
      
      pmer.initial_integrate();
      pmer.compute_tangents();    
      
      
      
      pmer.set_zetainv();
      pmer.update_Hhat();
      pmer.compute_uc_forces();
      pmer.compute_tension();

      pmer.final_integrate();
      pmer.compute_tangents();


      pmer.rescale_positions(true);

      
    }
    
    myfiles[molnum].open(filename.str());
    if (myfiles[molnum].is_open()) {
      Eigen::MatrixXd msave(pmer.R_x.size(),6);
      
      Eigen::VectorXd dummy_x(pmer.u_x.size()+1);
      Eigen::VectorXd dummy_y(pmer.u_y.size()+1);
      Eigen::VectorXd dummy_z(pmer.u_z.size()+1);
      
      Eigen::VectorXd dumb(1);
      dumb(0) = 0.0;
      
      dummy_x << pmer.u_x, dumb;
      dummy_y << pmer.u_y, dumb;
      dummy_z << pmer.u_z, dumb;
      
      msave << pmer.R_x, pmer.R_y, pmer.R_z, dummy_x, dummy_y, dummy_z;
      
      myfiles[molnum] << msave;
    }
    myfiles[molnum].close();
  }
  return 0;
  
}

