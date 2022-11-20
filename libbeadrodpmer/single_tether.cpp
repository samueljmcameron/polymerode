#include "single_tether.hpp"
#include "input.hpp"
#include "initialise.hpp"

#include <iostream>
#include <cmath>
#include <stdexcept>


using namespace BeadRodPmer;
/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
SingleTether::SingleTether(const std::vector<std::string> & splitvec)
  : Polymer(splitvec)
{

  if (!flag_x0) 
    throw std::runtime_error("Need to specify x0 for single tether polymer.");


  
  rhs_of_G.resize(Nbeads+2);
  dummy_for_noise.resize(Nbeads+2);
  Gmunu.resize(Nbeads+2,Nbeads+2);


  

  rhs_of_Hhat.resize(Nbeads+2);
  tension.resize(Nbeads+2);
  Hhat.resize(Nbeads+2,Nbeads+2);


  
  dCdlambda.resize(Nbeads+2,Nbeads+2);

  


  constraint_errors.resize(Nbeads+2);

  negative_tension_change.resize(Nbeads+2);
  
  tDets.resize(Nbeads);
  bDets.resize(Nbeads+3);

}


int SingleTether::single_step(double t,double dt,
			      const std::vector<std::vector<double>> & dFdX_i,
			      int itermax, int numtries,bool throw_exception)
{

  // get here if this step has been restarted numtry times

  // get here if this step has been restarted numtry times
  if (numtries == 0) {
    if (throw_exception)
      throw std::runtime_error("could not solve constraint equation for single tethered polymer.");
    else
      return -1;
  }

  
  set_unprojected_noise(dt);
  update_G();
  update_Hhat();
  compute_noise();
  compute_effective_kappa();
  compute_uc_forces();

  for (int index = 0; index < nuc_beads.size(); index ++ ) 
    add_external_force(dFdX_i[index],nuc_beads[index]);

  compute_tension({0,0,0});
  initial_integrate(dt,3,SINGLE);
  

  update_Hhat();
  compute_effective_kappa();
  compute_uc_forces();
  for (int index = 0; index < nuc_beads.size(); index ++ )
    add_external_force(dFdX_i[index],nuc_beads[index]);

  
  compute_tension({0,0,0});

  
  int iterations = correct_tension(dt,atoms[0].R,itermax,1e-8);
  if (iterations > itermax) {
    numtries -= 1;
    std::cout << "too many iterations when correcting tension at time " << t
	      <<  ", retrying the step with new noise ( " << numtries 
	      << " attempts left). " << std::endl;
    for (int i = 0; i < get_Nbeads(); i++) 
      atoms[i].R = Rtmp[i];

    compute_tangents_and_friction();
    return single_step(t,dt,dFdX_i,itermax,numtries,throw_exception);
  }
  else {
    final_integrate(dt,3,SINGLE);
  
    compute_tangents_and_friction();
  }

  return 0;
}



int SingleTether::single_step(double t,double dt,
			      const std::vector<std::vector<double>> & dFdX_i,
			      std::function<Eigen::Vector3d (double)> X0_t,
			      std::function<Eigen::Vector3d (double)> dX0dt,
			      int itermax, int numtries,bool throw_exception)
{

  // get here if this step has been restarted numtry times
  if (numtries == 0) {
    if (throw_exception)
      throw std::runtime_error("could not solve constraint equation for single tethered polymer.");
    else
      return -1;
  }

  
  set_unprojected_noise(dt);
  update_G();
  update_Hhat();
  compute_noise();
  compute_effective_kappa();
  compute_uc_forces();

  for (int index = 0; index < nuc_beads.size(); index ++ )
    add_external_force(dFdX_i[index],nuc_beads[index]);


  compute_tension(dX0dt(t+dt/2));
  initial_integrate(dt,3,SINGLE);
  

  update_Hhat();
  compute_effective_kappa();
  compute_uc_forces();
  for (int index = 0; index < nuc_beads.size(); index ++ )
    add_external_force(dFdX_i[index],nuc_beads[index]);

  
  compute_tension(dX0dt(t+dt));

  
  int iterations = correct_tension(dt,X0_t(t+dt),itermax,1e-8);
  if (iterations > itermax) {
    numtries -= 1;
    std::cout << "too many iterations when correcting tension at time " << t
	      <<  ", retrying the step with new noise ( " << numtries 
	      << " attempts left). " << std::endl;
    for (int i = 0; i < get_Nbeads(); i++) 
      atoms[i].R = Rtmp[i];

    compute_tangents_and_friction();
    return single_step(t,dt,dFdX_i,X0_t,dX0dt,itermax,numtries,throw_exception);
  }
  else {
    final_integrate(dt,3,SINGLE);
  
    compute_tangents_and_friction();
  }

  return 0;
}




/* ---------------------------------------------------------------------------- */
/* RHS of G*eta = P. */
/* ---------------------------------------------------------------------------- */
void SingleTether::set_rhs_of_G()
{
  int offset = 3;

  rhs_of_G({0,1,2}) = atoms[0].unprojected_noise;

  Polymer::set_rhs_of_G(offset);
  
  return;
  
}


/* -------------------------------------------------------------------------- */
/* Initialse G in G*eta = P (only call once). */
/* -------------------------------------------------------------------------- */
void SingleTether::set_G()
{

  std::vector<T> coefficients = init_G_coeffsmatrix();

  Gmunu.setFromTriplets(coefficients.begin(),coefficients.end());

  Gmunu_solver.analyzePattern(Gmunu);

  return;

}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix G. */
/* -------------------------------------------------------------------------- */
std::vector<T> SingleTether::init_G_coeffsmatrix()
{
  // just set lower diagonal
  int offset = 3;
  std::vector<T> coeffs;

  
  coeffs.push_back(T(offset-3,offset-3,1));
  coeffs.push_back(T(offset-2,offset-2,1));
  coeffs.push_back(T(offset-1,offset-1,1));


  Eigen::Vector3d btmp = bonds[0].rod;
  
  coeffs.push_back(T(offset,offset-3,-btmp(0)));
  coeffs.push_back(T(offset,offset-2,-btmp(1)));
  coeffs.push_back(T(offset,offset-1,-btmp(2)));

  Polymer::init_G_coeffsmatrix(offset,coeffs);
  
  return coeffs;

  
}



/* -------------------------------------------------------------------------- */
/* Update matrix G. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void SingleTether::update_G()
{
  // update first four cols by hand

  int offset = 3;

  Gmunu.coeffRef(offset,offset-3) = -bonds[0].rod(0);
  Gmunu.coeffRef(offset,offset-2) = -bonds[0].rod(1);
  Gmunu.coeffRef(offset,offset-1) = -bonds[0].rod(2);

  Polymer::update_G(offset);
  
  return;
}

void SingleTether::compute_effective_kappa()
{

  int offset = 3;


  bDets(Nbeads-1+offset) = 1.0;
  bDets(Nbeads-2+offset) = 2 - bonds[Nbeads-2].rod.squaredNorm();
  tDets(0) = 1.0;
  tDets(1) = 2 - bonds[0].rod.squaredNorm();

  set_bdets_and_tdets(offset);
  
  bDets(-1+offset) = bDets(offset)-bonds[0].rod(2)*bonds[0].rod(2)*bDets(1+offset);
  bDets(-2+offset) = bDets(-1+offset)-bonds[0].rod(1)*bonds[0].rod(1)*bDets(1+offset);
  bDets(-3+offset) = bDets(-2+offset)-bonds[0].rod(0)*bonds[0].rod(0)*bDets(1+offset);


  double gDet = bDets(-3+offset);
  
  end_inverses({0,1,2}) = bonds[0].rod*bDets(1+offset)/(gDet*bondlength);

  for (int i = 0; i < Nbeads-2; i++) {
    k_effs(i) = (kappa - temp*bondlength*costhetas(i)*tDets(i)*bDets(i+2+offset)/gDet
		 )/(bondlength*bondlength);
  }


  
  
  return;

}




/* ---------------------------------------------------------------------------- */
/* RHS of Hhat*lambda = Q. */
/* ---------------------------------------------------------------------------- */
void SingleTether::set_rhs_of_Hhat(const Eigen::Vector3d &dXdt_at_1)
{


  int offset = 3;


  rhs_of_Hhat({0,1,2}) = atoms[0].friction*(atoms[0].Fpot+atoms[0].noise)-dXdt_at_1; 
  
  Polymer::set_rhs_of_Hhat(offset);
  
  return;
  
}






/* -------------------------------------------------------------------------- */
/* Initialise matrix H hat. */
/* -------------------------------------------------------------------------- */
void SingleTether::set_Hhat()
{

  std::vector<T> coefficients = init_Hhat_coeffsmatrix();

  Hhat.setFromTriplets(coefficients.begin(),coefficients.end());


  Hhat_solver.analyzePattern(Hhat);

  return;

}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix H hat. */
/* -------------------------------------------------------------------------- */
std::vector<T> SingleTether::init_Hhat_coeffsmatrix()
{

  // only use lower triangular part as it's a symmetric matrix.  
  std::vector<T> coeffs;

  int offset = 3;

  coeffs.push_back(T(offset-3,offset-3,Hhat_endblocks(0,0,0)));
  coeffs.push_back(T(offset-2,offset-3,Hhat_endblocks(1,0,0)));
  coeffs.push_back(T(offset-2,offset-2,Hhat_endblocks(1,1,0)));
  coeffs.push_back(T(offset-1,offset-3,Hhat_endblocks(2,0,0)));
  coeffs.push_back(T(offset-1,offset-2,Hhat_endblocks(2,1,0)));
  coeffs.push_back(T(offset-1,offset-1,Hhat_endblocks(2,2,0)));


  
  coeffs.push_back(T(offset,offset-3,Hhat_leftside(0)));
  coeffs.push_back(T(offset,offset-2,Hhat_leftside(1)));
  coeffs.push_back(T(offset,offset-1,Hhat_leftside(2)));

  
  Polymer::init_Hhat_coeffsmatrix(offset,coeffs);

  return coeffs;

  
}






/* -------------------------------------------------------------------------- */
/* Update matrix M. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void SingleTether::update_Hhat()
{

  int offset = 3;
  // update first four columns by hand

  Hhat.coeffRef(offset-3,offset-3) = Hhat_endblocks(0,0,0);
  Hhat.coeffRef(offset-2,offset-3) = Hhat_endblocks(1,0,0);
  Hhat.coeffRef(offset-2,offset-2) = Hhat_endblocks(1,1,0);
  Hhat.coeffRef(offset-1,offset-3) = Hhat_endblocks(2,0,0);
  Hhat.coeffRef(offset-1,offset-2) = Hhat_endblocks(2,1,0);
  Hhat.coeffRef(offset-1,offset-1) = Hhat_endblocks(2,2,0);


  Hhat.coeffRef(offset,offset-3) = Hhat_leftside(0);
  Hhat.coeffRef(offset,offset-2) = Hhat_leftside(1);
  Hhat.coeffRef(offset,offset-1) = Hhat_leftside(2);

  Polymer::update_Hhat(offset);

  return;
}




/* -------------------------------------------------------------------------- */
/* Initialise both the dCdlambda and mat_tmp matrices. */
/* -------------------------------------------------------------------------- */
void SingleTether::set_dCdlambda()
{

  std::vector<T> coefficients = init_dCdlambda_coeffsmatrix();

  dCdlambda.setFromTriplets(coefficients.begin(),coefficients.end());

  jacob_solver.analyzePattern(dCdlambda);

  return;

}



/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix dCdlambda. */
/* -------------------------------------------------------------------------- */
std::vector<T> SingleTether::init_dCdlambda_coeffsmatrix()
{

  // initializing the vector as the full (both upper and lower) part of Hhat, since
  // this matrix won't be symmetric.
  std::vector<T> coeffs;

  int offset = 3;

  coeffs.push_back(T(offset-3,offset-3,Hhat_endblocks(0,0,0)));
  coeffs.push_back(T(offset-3,offset-2,Hhat_endblocks(1,0,0)));
  coeffs.push_back(T(offset-3,offset-1,Hhat_endblocks(2,0,0)));
  coeffs.push_back(T(offset-3,offset,Hhat_leftside(0)));
  
  
  coeffs.push_back(T(offset-2,offset-3,Hhat_endblocks(1,0,0)));
  coeffs.push_back(T(offset-2,offset-2,Hhat_endblocks(1,1,0)));
  coeffs.push_back(T(offset-2,offset-1,Hhat_endblocks(2,1,0)));
  coeffs.push_back(T(offset-2,offset,Hhat_leftside(1)));
  
  coeffs.push_back(T(offset-1,offset-3,Hhat_endblocks(2,0,0)));
  coeffs.push_back(T(offset-1,offset-2,Hhat_endblocks(2,1,0)));
  coeffs.push_back(T(offset-1,offset-1,Hhat_endblocks(2,2,0)));
  coeffs.push_back(T(offset-1,offset,Hhat_leftside(2)));

  
  coeffs.push_back(T(offset,offset-3,Hhat_leftside(0)));
  coeffs.push_back(T(offset,offset-2,Hhat_leftside(1)));
  coeffs.push_back(T(offset,offset-1,Hhat_leftside(2)));

  
  Polymer::init_dCdlambda_coeffsmatrix(offset,coeffs);

  return coeffs;


}






void SingleTether::update_dCdlambda(double Delta_t)
{
  // must update column-wise
  int offset = 3;


  // negative two col
  dCdlambda.coeffRef(offset-3,offset-3) = -Hhat_endblocks(0,0,0)*Delta_t;
  dCdlambda.coeffRef(offset-2,offset-3) = -Hhat_endblocks(1,0,0)*Delta_t;
  dCdlambda.coeffRef(offset-1,offset-3) = -Hhat_endblocks(2,0,0)*Delta_t;

  
  dCdlambda.coeffRef(offset,offset-3) = -dCdlambda_leftside(0)*Delta_t;


  // n one col
  dCdlambda.coeffRef(offset-3,offset-2) = -Hhat_endblocks(1,0,0)*Delta_t;
  dCdlambda.coeffRef(offset-2,offset-2) = -Hhat_endblocks(1,1,0)*Delta_t;
  dCdlambda.coeffRef(offset-1,offset-2) = -Hhat_endblocks(2,1,0)*Delta_t;
  dCdlambda.coeffRef(offset,offset-2) = -dCdlambda_leftside(1)*Delta_t;


  
  // zero col
  dCdlambda.coeffRef(offset-3,offset-1) = -Hhat_endblocks(2,0,0)*Delta_t;
  dCdlambda.coeffRef(offset-2,offset-1) = -Hhat_endblocks(2,1,0)*Delta_t;
  dCdlambda.coeffRef(offset-1,offset-1) = -Hhat_endblocks(2,2,0)*Delta_t;
  dCdlambda.coeffRef(offset,offset-1) = -dCdlambda_leftside(2)*Delta_t;



  // one col
  dCdlambda.coeffRef(offset-3,offset) = -Hhat_leftside(0)*Delta_t;
  dCdlambda.coeffRef(offset-2,offset) = -Hhat_leftside(1)*Delta_t;
  dCdlambda.coeffRef(offset-1,offset) = -Hhat_leftside(2)*Delta_t;

  Polymer::update_dCdlambda(Delta_t,offset);
  

  return;

}



void SingleTether::compute_noise()
{

  int offset = 3;
  set_rhs_of_G();
  Gmunu_solver.factorize(Gmunu);

  dummy_for_noise =  Gmunu_solver.solve(rhs_of_G);

  Polymer::update_noise(offset);

  atoms[0].noise -= dummy_for_noise({0,1,2});


  return;
}


void SingleTether::compute_tension(const Eigen::Vector3d & dXdt_at_1)
{



  set_rhs_of_Hhat(dXdt_at_1);

  Hhat_solver.factorize(Hhat);
  tension =  Hhat_solver.solve(rhs_of_Hhat);
  
}

void SingleTether::test_jacob(int mu,double Delta_t,const Eigen::Vector3d & X_of_t_at_1)
{

  int offset = 3;  
  // start by initializing the vector dC at whatever tensions
  
  
  calculate_constraint_errors(X_of_t_at_1);


  update_dCdlambda(Delta_t);

  
  //  std::cout << "exact calc: " << std::endl;
  //  std::cout << dCdlambda.col(mu+offset-1) << std::endl;

  Eigen::VectorXd zeroval_change = constraint_errors;

 
  double tens_change_old = 0;
  for (double tens_change = 0.00001; tens_change > 0.000009; tens_change /= 2) {

    tension(mu+offset-1) -= tens_change_old;
    tension(mu+offset-1) += tens_change;
    
    final_integrate(Delta_t,3,SINGLE);

    calculate_constraint_errors(X_of_t_at_1);
    
    std::cout << "numerical estimate with dlamda =  " << tens_change
	      << "and lambda = " << tension(mu+offset-1) << " is = " << std::endl;
    std::cout << (constraint_errors-zeroval_change)/tens_change-dCdlambda.col(mu+offset-1)
	      << std::endl;

    tens_change_old = tens_change;

  }
  
  return;
}  

int SingleTether::correct_tension(double Delta_t,const Eigen::Vector3d & X_of_t_at_1,
				  int itermax,double tolerance)
{

  // set C_mu and dC_mu/dlambda_nu
  final_integrate(Delta_t,3,SINGLE);
  calculate_constraint_errors(X_of_t_at_1);
  update_dCdlambda(Delta_t);

  
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

    final_integrate(Delta_t,3,SINGLE);
    calculate_constraint_errors(X_of_t_at_1);
    update_dCdlambda(Delta_t);


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





void SingleTether::calculate_constraint_errors(const Eigen::Vector3d & X_of_t_at_1)
{
  int offset = 3;
  constraint_errors({0,1,2}) = atoms[0].R-X_of_t_at_1;

  Polymer::calculate_constraint_errors(offset);

  return;

}
