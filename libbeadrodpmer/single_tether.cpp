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
  : NoTether(splitvec)
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

int SingleTether::single_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			      Eigen::Ref<Eigen::Matrix3Xd> Fs,double t,double dt,
			      const std::vector<Eigen::Vector3d> & dFdX_i,
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


  Fs.setZero();
  for (int index = 0; index < nuc_beads.size(); index ++ ) 
    Fs.col(nuc_beads[index]) += -dFdX_i[index];

  first_step(xs,Fs,dt);

  
  Fs.setZero();
  for (int index = 0; index < nuc_beads.size(); index++)
    Fs.col(nuc_beads[index]) += -dFdX_i[index];  
  
  int iterations = second_step(xs,Fs,dt,itermax);
  
  if (iterations > itermax) {
    numtries -= 1;
    std::cout << "too many iterations when correcting tension at time " << t
	      <<  ", retrying the step with new noise ( " << numtries 
	      << " attempts left). " << std::endl;
    return single_step(xs,Fs,t,dt,dFdX_i,itermax,numtries,throw_exception);
  }

  return 0;
}


int SingleTether::single_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			      Eigen::Ref<Eigen::Matrix3Xd> Fs,double t,double dt,
			      const std::vector<Eigen::Vector3d> & dFdX_i,
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

  Fs.setZero();
  for (int index = 0; index < nuc_beads.size(); index++)
    Fs.col(nuc_beads[index]) += -dFdX_i[index];
  

  first_step(xs,Fs,dt,dX0dt(t+dt/2));
  
  Fs.setZero();
  
  for (int index = 0; index < nuc_beads.size(); index++)
    Fs.col(nuc_beads[index]) += -dFdX_i[index];  

  int iterations = second_step(xs,Fs,dt,itermax,X0_t(t+dt),dX0dt(t+dt/2));
  

  if (iterations > itermax) {
    numtries -= 1;
    std::cout << "too many iterations when correcting tension at time " << t
	      <<  ", retrying the step with new noise ( " << numtries 
	      << " attempts left). " << std::endl;
    return single_step(xs,Fs,t,dt,dFdX_i,X0_t,dX0dt,itermax,numtries,throw_exception);
  }

  return 0;
}



/* ============================================================================ */
/* Initial integration step. The force vector Fs needs to have been set (either
   by particle-particle interactions or set to zero). */
/* ============================================================================ */
void SingleTether::first_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			      Eigen::Ref<Eigen::Matrix3Xd> Fs,double dt)
{
  set_unprojected_noise(dt);
  update_G();
  update_Hhat();
  compute_noise();
  compute_effective_kappa();
  compute_uc_forces(Fs);
  
  compute_tension({0,0,0},Fs);
  initial_integrate(xs,Fs,dt,3,SINGLE);

}


/* ============================================================================ */
/* Initial integration step. The force vector Fs needs to have been set (either
   by particle-particle interactions or set to zero). */
/* ============================================================================ */
void SingleTether::first_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			      Eigen::Ref<Eigen::Matrix3Xd> Fs,double dt,
			      const Eigen::Ref<const Eigen::Vector3d> &dX0dt)
{
  set_unprojected_noise(dt);
  update_G();
  update_Hhat();
  compute_noise();
  compute_effective_kappa();
  compute_uc_forces(Fs);
  
  compute_tension(dX0dt,Fs);
  initial_integrate(xs,Fs,dt,3,SINGLE);

}


/* ============================================================================ */
/* Second integration step. The force vector Fs needs to have been set (either
   by particle-particle interactions or set to zero). */
/* ============================================================================ */
int SingleTether::second_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			      Eigen::Ref<Eigen::Matrix3Xd> Fs,
			      double dt,int itermax)
{
  update_Hhat();
  compute_effective_kappa();
  compute_uc_forces(Fs);
  
  compute_tension({0,0,0},Fs);

  int iterations = correct_tension(xs,Fs,dt,xs.col(0),itermax,1e-8);


  if (iterations > itermax) {
    for (int i = 0; i < get_Nbeads(); i++) 
      xs.col(i) = tmp_xs.col(i);
  } else {
    final_integrate(xs,Fs,dt,3,SINGLE);
  }
  compute_tangents_and_friction(xs);
  
  return iterations;
}



/* ============================================================================ */
/* Second integration step. The force vector Fs needs to have been set (either
   by particle-particle interactions or set to zero). */
/* ============================================================================ */
int SingleTether::second_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			      Eigen::Ref<Eigen::Matrix3Xd> Fs,
			      double dt,int itermax,
			      const Eigen::Ref<const Eigen::Vector3d> &X0_tFull,
			      const Eigen::Ref<const Eigen::Vector3d> &dX0dtHalf)
{
  update_Hhat();
  compute_effective_kappa();
  compute_uc_forces(Fs);
  
  compute_tension(dX0dtHalf,Fs);

  int iterations = correct_tension(xs,Fs,dt,X0_tFull,itermax,1e-8);


  if (iterations > itermax) {
    for (int i = 0; i < get_Nbeads(); i++) 
      xs.col(i) = tmp_xs.col(i);
  } else {
    final_integrate(xs,Fs,dt,3,SINGLE);
  }
  compute_tangents_and_friction(xs);
  
  return iterations;
}




/* ---------------------------------------------------------------------------- */
/* RHS of G*eta = P. */
/* ---------------------------------------------------------------------------- */
void SingleTether::set_rhs_of_G()
{
  int offset = 3;

  rhs_of_G({0,1,2}) = unprojected_noises.col(0);

  NoTether::set_rhs_of_G(offset);
  
  return;
  
}


/* -------------------------------------------------------------------------- */
/* Initialse G in G*eta = P (only call once). */
/* -------------------------------------------------------------------------- */
void SingleTether::set_G()
{

  std::vector<T> coefficients = init_G_coeffsmatrix();

  Gmunu.setFromTriplets(coefficients.begin(),coefficients.end());

  Gmunu_solver->analyzePattern(Gmunu);


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


  Eigen::Vector3d btmp = bonds.col(0);
  
  coeffs.push_back(T(offset,offset-3,-btmp(0)));
  coeffs.push_back(T(offset,offset-2,-btmp(1)));
  coeffs.push_back(T(offset,offset-1,-btmp(2)));

  NoTether::init_G_coeffsmatrix(offset,coeffs);
  
  return coeffs;

  
}



/* -------------------------------------------------------------------------- */
/* Update matrix G. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void SingleTether::update_G()
{
  // update first four cols by hand

  int offset = 3;

  Gmunu.coeffRef(offset,offset-3) = -bonds(0,0);
  Gmunu.coeffRef(offset,offset-2) = -bonds(1,0);
  Gmunu.coeffRef(offset,offset-1) = -bonds(2,0);

  NoTether::update_G(offset);

  return;
}


void SingleTether::set_bdets_and_tdets()
{
  int offset = 3;


  bDets(Nbeads-1+offset) = 1.0;
  bDets(Nbeads-2+offset) = 2 - bonds.col(Nbeads-2).squaredNorm();
  tDets(0) = 1.0;
  tDets(1) = 2 - bonds.col(0).squaredNorm();

  NoTether::set_bdets_and_tdets(offset);
  
  bDets(-1+offset) = bDets(offset)-bonds(2,0)*bonds(2,0)*bDets(1+offset);
  bDets(-2+offset) = bDets(-1+offset)-bonds(1,0)*bonds(1,0)*bDets(1+offset);
  bDets(-3+offset) = bDets(-2+offset)-bonds(0,0)*bonds(0,0)*bDets(1+offset);

  return;
}

void SingleTether::compute_effective_kappa()
{

  int offset = 3;

  set_bdets_and_tdets();
  
  double gDet = bDets(-3+offset);
  
  end_inverses({0,1,2}) = bonds.col(0)*bDets(1+offset)/(gDet*bondlength);

  for (int i = 0; i < Nbeads-2; i++) {
    k_effs(i) = (kappa - temp*bondlength*costhetas(i)*tDets(i)*bDets(i+2+offset)/gDet
		 )/(bondlength*bondlength);
  }

  return;

}




/* ---------------------------------------------------------------------------- */
/* RHS of Hhat*lambda = Q. */
/* ---------------------------------------------------------------------------- */
void SingleTether::set_rhs_of_Hhat(const Eigen::Vector3d &dXdt_at_1,
				   const Eigen::Ref<const Eigen::Matrix3Xd> &Fs)
{


  int offset = 3;


  rhs_of_Hhat({0,1,2}) = frictions[0]*(Fs.col(0)+noises.col(0))-dXdt_at_1; 
  
  NoTether::set_rhs_of_Hhat(offset,Fs);
  
  return;
  
}






/* -------------------------------------------------------------------------- */
/* Initialise matrix H hat. */
/* -------------------------------------------------------------------------- */
void SingleTether::set_Hhat()
{

  std::vector<T> coefficients = init_Hhat_coeffsmatrix();

  Hhat.setFromTriplets(coefficients.begin(),coefficients.end());


  Hhat_solver->analyzePattern(Hhat);

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

  
  NoTether::init_Hhat_coeffsmatrix(offset,coeffs);

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

  NoTether::update_Hhat(offset);

  return;
}




/* -------------------------------------------------------------------------- */
/* Initialise both the dCdlambda and mat_tmp matrices. */
/* -------------------------------------------------------------------------- */
void SingleTether::set_dCdlambda()
{

  std::vector<T> coefficients = init_dCdlambda_coeffsmatrix();

  dCdlambda.setFromTriplets(coefficients.begin(),coefficients.end());

  jacob_solver->analyzePattern(dCdlambda);

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

  
  NoTether::init_dCdlambda_coeffsmatrix(offset,coeffs);

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

  NoTether::update_dCdlambda(Delta_t,offset);
  

  return;

}



void SingleTether::compute_noise()
{

  int offset = 3;
  set_rhs_of_G();
  Gmunu_solver->factorize(Gmunu);

  dummy_for_noise =  Gmunu_solver->solve(rhs_of_G);

  NoTether::update_noise(offset);

  noises.col(0) -= dummy_for_noise({0,1,2});


  return;
}


void SingleTether::compute_tension(const Eigen::Vector3d & dXdt_at_1,
				   const Eigen::Ref<const Eigen::Matrix3Xd> & Fs)
{



  set_rhs_of_Hhat(dXdt_at_1,Fs);

  Hhat_solver->factorize(Hhat);
  tension =  Hhat_solver->solve(rhs_of_Hhat);
  
}
/*
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
*/
int SingleTether::correct_tension(Eigen::Ref<Eigen::Matrix3Xd> xs,
				  const Eigen::Ref<const Eigen::Matrix3Xd> & Fs,
				  double Delta_t,const Eigen::Vector3d & X_of_t_at_1,
				  int itermax,double tolerance)
{

  // set C_mu and dC_mu/dlambda_nu
  final_integrate(xs,Fs,Delta_t,3,SINGLE);
  calculate_constraint_errors(X_of_t_at_1,xs);
  update_dCdlambda(Delta_t);

  
  //and then solve

  jacob_solver->factorize(dCdlambda);

  if (jacob_solver->info() != Eigen::Success)  {
    std::cout << "Matrix factorization failed in tension correction on first step."
	      << std::endl;
    return itermax + 1;
  }

  
  negative_tension_change = jacob_solver->solve(constraint_errors);
  
  tension = tension - negative_tension_change;
  int count = 0;


  while (negative_tension_change.norm() > tolerance && count <= itermax) {

    final_integrate(xs,Fs,Delta_t,3,SINGLE);
    calculate_constraint_errors(X_of_t_at_1,xs);
    update_dCdlambda(Delta_t);


    jacob_solver->factorize(dCdlambda);


    if (jacob_solver->info() != Eigen::Success)  {
      std::cout << "Matrix factorization failed in tension correction."
		<< std::endl;
      return itermax + 1;
    }
      
    negative_tension_change = jacob_solver->solve(constraint_errors);

    
    tension = tension - negative_tension_change;

    count += 1;
  }
  return count;
  
}





void SingleTether::calculate_constraint_errors(const Eigen::Vector3d & X_of_t_at_1,
					       const Eigen::Ref<const Eigen::Matrix3Xd> &xs)
{
  int offset = 3;
  constraint_errors({0,1,2}) = xs.col(0)-X_of_t_at_1;

  NoTether::calculate_constraint_errors(offset,xs);

  return;

}
