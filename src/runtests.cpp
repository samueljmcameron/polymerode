#include "runtests.hpp"

#include "testclass.hpp"


void test_noise(GlobalParams &gp, Polymer &pmer)
{


  int N = pmer.get_Nbeads();
  double dt = gp.timestep;

  pmer.compute_tangents_and_rods_and_friction();

  TestClass test(pmer.atoms,pmer.bonds);


  /* Compare the matrix G in the two cases. */
  
  std::cout << "Comparison of G :" << std::endl;
  
  pmer.set_G();
  TestClass::SpMat Ghat = test.compute_Ghat();
  std::cout << pmer.Gmunu - Ghat << std::endl;


  /* Compare the geometrically projected noise in the two cases. */
  
  std::cout << "Geometrically projected noises comparison: " << std::endl;
  
  pmer.set_unprojected_noise(dt);
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
  pmer.set_rhs_of_Hhat(dt/2);

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
  pmer.compute_tension(dt/2);
  pmer.initial_integrate(dt);

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

  std::cout << pmer.get_timescale(dt) << std::endl;


  pmer.compute_tangents_and_rods_and_friction();

  test.atoms = pmer.atoms;
  test.bonds = pmer.bonds;

  pmer.update_Hhat();

  //  std::cout << pmer.Hhat << std::endl;
  //  std::cout << test.compute_Hhat() << std::endl;
  pmer.compute_uc_forces();
  pmer.compute_tension(dt);
  pmer.final_integrate(dt);


  std::cout << "new bond lengths final value: " << std::endl;


  for (int mu = 0; mu < N-1; mu++) {
    std::cout << sqrt((pmer.atoms[mu+1].R-pmer.atoms[mu].R).dot(pmer.atoms[mu+1].R-pmer.atoms[mu].R))/pmer.get_bondlength()
	      << std::endl;
  }
    
  return;
  
}

