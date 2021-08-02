
#include "polymer.hpp"
#include "testclass.hpp"
#include "matrix.hpp"
#include "input.hpp"
#include <Eigen/SparseLU>

#include <iostream>
#include <fstream>
#include <string>

#include <chrono>

int main(int argc, char* argv[])
{

  std::ifstream infile;

  double dt = 0.0;
  

  int iarg = 1;  
  while(iarg < argc) {
    if (strcmp(argv[iarg],"-in") == 0) {
      if (iarg+1 == argc) {
	std::cerr << "Error: input flag '-in' specified, but no file given."
		  << std::endl;
	return EXIT_FAILURE;
      }
      infile.open(argv[iarg+1]);
      iarg += 2;
      
    } else if (strcmp(argv[iarg],"-var") == 0) {
      
      if (iarg + 2 >= argc) {
	std::cerr << "Error: invalid command line variable specification."
		  << std::endl;
	return EXIT_FAILURE;
      }
      variables[argv[iarg+1]] = argv[iarg+2];
      iarg += 3;
    } else {
      std::cerr << "Error: invalid command line variable specification."
		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  if (not infile.is_open()) {
    std::cerr < "Error: need to specify input file." << std::endl;
    return EXIT_FAILURE;
  }


  std::string line;  

  std::vector<std::string> splitvec;
  
  line = line.substr(0,line.find_first_of("#"));
  splitvec = input::split_line(line);
  

  while ( std::getline(infile,line) &&
	  splitvec[0] != "build_polymer")  {

    std::getline(infile,line);
    line = line.substr(0,line.find_first_of("#"));
    splitvec = input::split_line(line);
    
  }

  
  splitvec.erase(splitvec.begin());

  Polymer pmer(splitvec);

  pmer.compute_tangents_and_rods_and_friction();

  TestClass test(pmer.atoms,pmer.bonds);


  /* Compare the matrix G in the two cases. */
  
  std::cout << "Comparison of G :" << std::endl;
  
  pmer.set_G();
  TestClass::SpMat Ghat = test.compute_Ghat();
  std::cout << pmer.Gmunu - Ghat << std::endl;


  /* Compare the geometrically projected noise in the two cases. */
  
  std::cout << "Geometrically projected noises comparison: " << std::endl;
  
  pmer.set_unprojected_noise();
  pmer.compute_noise();

  Eigen::MatrixXd Ginv = test.invert_matrix(Ghat);  
  Matrix<Eigen::Vector3d> boldn = test.set_boldn_matrix();
  Matrix<Eigen::Matrix3d> projection = test.geo_projection(Ginv,boldn);
  
  for (int i = 0; i < N; i++) {
    Eigen::Vector3d noise(0,0,0);
    for (int j = 0; j < N; j++) {
      noise += projection(i,j)*pmer.atoms[j].unprojected_noise;
    }
    std::cout << "Atom " << i << ": "
	      << noise  - pmer.atoms[i].noise
	      << std::endl;
  }
  
  std::cout << "Explicitly check noise projection: " << std::endl;


  for (int mu = 0; mu < N-1; mu++) {
    double sum = 0;
    for (int i = 0; i < N; i++) {
      sum += pmer.atoms[i].noise.dot(boldn(i,mu));
    }
    std::cout << "Bond " << mu << ": "
	      << sum << std::endl;
  }

  /* Now check RHS of g eqn. */

  
  pmer.set_rhs_of_G();
  
  std::cout << "RHS of G eqn: " << std::endl;

  for (int mu = 0; mu < N -1; mu++ ) {
    double sum = 0;
    for (int j = 0; j < N; j++) {
      sum += boldn(j,mu).dot(pmer.atoms[j].unprojected_noise);
    }
    std::cout << "Bond " << mu << ": " << pmer.rhs_of_G(mu) - sum
	      << std::endl;
  }
  
  

  /* Compare the H matrix in the two cases. */
  
  std::cout << "Comparison of H :" << std::endl;

  pmer.set_Hhat();

  TestClass::SpMat Hhat = test.compute_Hhat();

  std::cout << pmer.Hhat - Hhat << std::endl;

  
  /* Then check the RHS of the H eqn . */

  pmer.compute_uc_forces();
  pmer.set_rhs_of_Hhat();

  std::cout << "RHS of H eqn: " << std::endl;

  for (int mu = 0; mu < N -1; mu++ ) {
    double sum = 0;
    for (int j = 0; j < N; j++) {
      sum += boldn(j,mu).dot(pmer.atoms[j].friction*pmer.atoms[j].F);
    }
    std::cout << "Bond " << mu << ": " << pmer.rhs_of_Hhat(mu) - sum
	      << std::endl;
  }
  
  
  /* Compare the dynamic projection in the two cases. */
  pmer.compute_tension();
  pmer.initial_integrate();

  Eigen::MatrixXd Hinv = test.invert_matrix(Hhat);  
  projection = test.dyn_projection(Ginv,boldn);
  std::vector<Eigen::Vector3d> Rdots(N);
  
  for (int i = 0; i < N; i++) {
    Rdots[i].setZero();
    for (int j = 0; j < N; j++) {
      Rdots[i] += projection(i,j)*pmer.atoms[j].friction*pmer.atoms[j].F;
    }
  }
  
  std::cout << "Constraint checks: " << std::endl;

  
  for (int mu = 0; mu < N-1; mu++) {
    double sum = 0;
    for (int i = 0; i < N; i++) {
      sum += boldn(i,mu).dot(Rdots[i]);
    }
    std::cout << "Bond " << mu << "constraint: "
	      << sum << std::endl;
  }

  std::cout << "R comparison in the two cases: " << std::endl;

  std::vector<Eigen::Vector3d> Rtmps(N);
  for (int i = 0; i < N; i++) {
    Rtmps[i] = pmer.Rtmp[i]+dt/2.0*Rdots[i];
    std::cout //<< pmer.atoms[i].R << " " << Rtmps[i] << " "
	      << pmer.atoms[i].R -Rtmps[i] << std::endl;
  }

  std::cout << dt*kbT/(rodin*zperp) << std::endl;


  pmer.compute_tangents_and_rods_and_friction();

  test.atoms = pmer.atoms;
  test.bonds = pmer.bonds;

  pmer.update_Hhat();

  //  std::cout << pmer.Hhat << std::endl;
  //  std::cout << test.compute_Hhat() << std::endl;
  pmer.compute_uc_forces();
  pmer.compute_tension();
  pmer.final_integrate();


  std::cout << "new bond lengths final value: " << std::endl;


  for (int mu = 0; mu < N-1; mu++) {
    std::cout << sqrt((pmer.atoms[mu+1].R-pmer.atoms[mu].R).dot(pmer.atoms[mu+1].R-pmer.atoms[mu].R))/rodin
	      << std::endl;
  }
    
  return 0;
  
}

