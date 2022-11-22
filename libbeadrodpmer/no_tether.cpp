#include "no_tether.hpp"
#include "input.hpp"

#include <iostream>
#include <cmath>
#include <stdexcept>


namespace BeadRodPmer {

/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
NoTether::NoTether(const std::vector<std::string> & splitvec)
  : Polymer(splitvec)
{

	      
  rhs_of_G.resize(Nbeads-1);
  dummy_for_noise.resize(Nbeads-1);
  Gmunu.resize(Nbeads-1,Nbeads-1);



  rhs_of_Hhat.resize(Nbeads-1);
  tension.resize(Nbeads-1);
  Hhat.resize(Nbeads-1,Nbeads-1);


  
  dCdlambda.resize(Nbeads-1,Nbeads-1);

  


  constraint_errors.resize(Nbeads-1);

  negative_tension_change.resize(Nbeads-1);
  
  tDets.resize(Nbeads);
  bDets.resize(Nbeads);

  

  // set ^ configuration, with the bottoms of the ^ being x0 and xN
  init_atoms_caret();

}



void NoTether::setup() {
  compute_tangents_and_friction();
  set_G();
  set_Hhat();
  set_dCdlambda();
}

int NoTether::single_step(double t, double dt,
			  const std::vector<Eigen::Vector3d> & dFdX_i,
			  int itermax, int numtries,bool throw_exception)
{

  // get here if this step has been restarted numtry times
  if (numtries == 0) {
    if (throw_exception)
      throw std::runtime_error("could not solve constraint equation for no tether polymer.");
    else
      return -1;
  }

  
  set_unprojected_noise(dt);
  update_G(0);
  update_Hhat(0);
  compute_noise();
  compute_effective_kappa();
  compute_uc_forces();

  for (int index = 0; index < nuc_beads.size(); index ++ )
    add_external_force(dFdX_i[index],nuc_beads[index]);
  
  compute_tension();
  initial_integrate(dt,0,NONE);
  

  
  update_Hhat(0);
  compute_effective_kappa();
  compute_uc_forces();
  
  for (int index = 0; index < nuc_beads.size(); index ++ )
    add_external_force(dFdX_i[index],nuc_beads[index]);
  
  
  compute_tension();

  
  int iterations = correct_tension(dt,itermax,1e-8);

  if (iterations > itermax) {
    numtries -= 1;
    std::cout << "too many iterations when correcting tension at time " << t
	      <<  ", retrying the step with new noise ( " << numtries 
	      << " attempts left). " << std::endl;
    for (int i = 0; i < get_Nbeads(); i++) 
      xs.col(i) = tmp_xs.col(i);

    compute_tangents_and_friction();
    return single_step(t,dt,dFdX_i,itermax,numtries,throw_exception);
  } else {
    final_integrate(dt,0,NONE);
  
    compute_tangents_and_friction();
  }

  return 0;
}







/* -------------------------------------------------------------------------- */
/* Initialse G in G*eta = P (only call once). */
/* -------------------------------------------------------------------------- */
void NoTether::set_G()
{

  std::vector<T> coefficients;
  init_G_coeffsmatrix(0,coefficients);

  Gmunu.setFromTriplets(coefficients.begin(),coefficients.end());

  Gmunu_solver.analyzePattern(Gmunu);

  return;

}




void NoTether::compute_effective_kappa()
{


  bDets(Nbeads-1) = 1.0;
  tDets(0) = 1.0;

  bDets(Nbeads-2) = 2;
  tDets(1) = 2;

  set_bdets_and_tdets(0);

  double gDet = tDets(Nbeads-1);
  
  for (int i = 0; i < Nbeads-2; i++) {

    k_effs(i) = (kappa - temp*bondlength*costhetas(i)*tDets(i)*bDets(i+2)/gDet
		 )/(bondlength*bondlength);
  }


  
  
  return;

}





/* -------------------------------------------------------------------------- */
/* Initialise matrix H hat. */
/* -------------------------------------------------------------------------- */
void NoTether::set_Hhat()
{

  std::vector<T> coefficients;
  init_Hhat_coeffsmatrix(0,coefficients);

  Hhat.setFromTriplets(coefficients.begin(),coefficients.end());


  Hhat_solver.analyzePattern(Hhat);

  return;

}





/* -------------------------------------------------------------------------- */
/* Initialise both the dCdlambda and mat_tmp matrices. */
/* -------------------------------------------------------------------------- */
void NoTether::set_dCdlambda()
{

  std::vector<T> coefficients;
  
  init_dCdlambda_coeffsmatrix(0,coefficients);

  dCdlambda.setFromTriplets(coefficients.begin(),coefficients.end());

  jacob_solver.analyzePattern(dCdlambda);

  return;

}


void NoTether::compute_noise()
{


  set_rhs_of_G(0);
  Gmunu_solver.factorize(Gmunu);

  dummy_for_noise =  Gmunu_solver.solve(rhs_of_G);

  update_noise(0);
  
}


void NoTether::compute_tension()
{



  set_rhs_of_Hhat(0);

  Hhat_solver.factorize(Hhat);
  tension =  Hhat_solver.solve(rhs_of_Hhat);
  
}

void NoTether::test_jacob(int mu,double Delta_t)
{


  // start by initializing the vector dC at whatever tensions
  
  
  calculate_constraint_errors(0);


  update_dCdlambda(Delta_t,0);

  
  //  std::cout << "exact calc: " << std::endl;
  //  std::cout << dCdlambda.col(mu) << std::endl;

  Eigen::VectorXd zeroval_change = constraint_errors;

 
  double tens_change_old = 0;
  for (double tens_change = 0.00001; tens_change > 0.000009; tens_change /= 2) {

    tension(mu) -= tens_change_old;
    tension(mu) += tens_change;
    
    final_integrate(Delta_t,0,NONE);

    calculate_constraint_errors(0);
    
    std::cout << "numerical estimate with dlamda =  " << tens_change
	      << "and lambda = " << tension(mu) << " is = " << std::endl;
    std::cout << (constraint_errors-zeroval_change)/tens_change-dCdlambda.col(mu)
	      << std::endl;

    tens_change_old = tens_change;

  }
  
  return;
}  

int NoTether::correct_tension(double Delta_t,int itermax,double tolerance)
{

  // set C_mu and dC_mu/dlambda_nu
  final_integrate(Delta_t,0,NONE);
  calculate_constraint_errors(0);
  update_dCdlambda(Delta_t,0);

  
  //and then solve

  jacob_solver.factorize(dCdlambda);


  if (jacob_solver.info() != Eigen::Success)  {
    std::cout << "Matrix factorization failed in tension correction on first step."
	      << std::endl;
    return itermax + 1;
  }

  negative_tension_change = jacob_solver.solve(constraint_errors);
  
  tension = tension - negative_tension_change;

  int count = 0;

  while (negative_tension_change.norm() > tolerance && count <= itermax) {



    final_integrate(Delta_t,0,NONE);
    calculate_constraint_errors(0);
    update_dCdlambda(Delta_t,0);

  
    jacob_solver.factorize(dCdlambda);


    if (jacob_solver.info() != Eigen::Success)  {
      std::cout << "Matrix factorization failed in tension correction."
		<< std::endl;
      return itermax + 1;
    }

    
    negative_tension_change = jacob_solver.solve(constraint_errors);

    
    tension = tension - negative_tension_change;
    count += 1;
  }
  return count;
  
}
 
}
