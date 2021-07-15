
#include "testclass.hpp"
#include <Eigen/SparseLU>

#include <iostream>
#include <fstream>
#include <string>

#include <chrono>

int main(int argc, char* argv[])
{
  int N = 10;
  double rodin = 0.1;
  double dt = 0.01;
  double zpara = 1.0;
  double zperp = 0.5;
  double kbT = 1.0;
  double kappa = 1.0;

  int seed = 124890;
  
  Eigen::Vector3d x0(0.0,0.0,0.0);


  TestClass pmer(N,rodin,dt,zpara,zperp,kbT,kappa,x0,seed);

  pmer.compute_tangents();
  pmer.set_G();
  pmer.set_unprojected_noise();
  std::cout << "Unprojected noise: " << pmer.eta_x << std::endl;

  Eigen::VectorXd noise_x = pmer.eta_x;
  Eigen::VectorXd noise_y = pmer.eta_y;
  Eigen::VectorXd noise_z = pmer.eta_z;
  
  pmer.compute_noise();

  std::cout << "Dummy vector: " << pmer.dummy_for_noise << std::endl;

  std::cout << "Projected noise: " << pmer.eta_x << std::endl;


  Eigen::VectorXd projnoise_x = noise_x;
  Eigen::VectorXd projnoise_y = noise_y;
  Eigen::VectorXd projnoise_z = noise_z;
  Eigen::Matrix3d projection;
  for (int i = 0; i < N; i++) {
    projnoise_x(i) = 0.0;
    projnoise_y(i) = 0.0;
    projnoise_z(i) = 0.0;
    for (int j = 0; j < N; j++) {

      projection = pmer.projection(i,j);
      projnoise_x(i) += projection(0,0)*noise_x(j)+projection(0,1)*noise_y(j)+projection(0,2)*noise_z(j);
      projnoise_y(i) += projection(1,0)*noise_x(j)+projection(1,1)*noise_y(j)+projection(1,2)*noise_z(j);
      projnoise_z(i) += projection(2,0)*noise_x(j)+projection(2,1)*noise_y(j)+projection(2,2)*noise_z(j);

      
    }
  }
  std::cout << "Unprojected noise in x: " << noise_x << std::endl;
  std::cout << "Projected noise in x: " << pmer.eta_x << std::endl;
  std::cout << "Matrix calculated projection noise in x: " << projnoise_x << std::endl;

  std::cout << "Unprojected noise in y: " << noise_y << std::endl;
  std::cout << "Projected noise in y: " << pmer.eta_y << std::endl;
  std::cout << "Matrix calculated projection noise in y: " << projnoise_y << std::endl;
  
  std::cout << "Unprojected noise in z: " << noise_z << std::endl;
  std::cout << "Projected noise in z: " << pmer.eta_z << std::endl;
  std::cout << "Matrix calculated projection noise in z: " << projnoise_z << std::endl;
  
  
  

  std::cout << "Dummy vector differences: " << pmer.dummy_for_noise(2)-pmer.dummy_for_noise(1)
	    << std::endl << pmer.dummy_for_noise(3)-pmer.dummy_for_noise(2)
    	    << std::endl << pmer.dummy_for_noise(4)-pmer.dummy_for_noise(3)
	    << std::endl << pmer.dummy_for_noise(5)-pmer.dummy_for_noise(4)
	    << std::endl;

  double tmp;

  std::cout << "Test explicitly if noise is projected correctly: " << std::endl;
  for (int mu = 0; mu < N-1; mu++) {
   tmp = ((pmer.eta_x(mu+1)-pmer.eta_x(mu))*pmer.u_x(mu)
	   +(pmer.eta_y(mu+1)-pmer.eta_y(mu))*pmer.u_y(mu)
	   +(pmer.eta_z(mu+1)-pmer.eta_z(mu))*pmer.u_z(mu));


    std::cout << "mu = " << mu << ", eqn gives " << tmp << std::endl;
  }

  /*
  std::ofstream file ("test.txt");
  if (file.is_open()) {
    Eigen::MatrixXd msave(pmer.eta_x.size(),3);

    msave << pmer.eta_x, pmer.eta_y, pmer.eta_z;

    file << msave;
  }
  */
  
  return 0;
  
}

