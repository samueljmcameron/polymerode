#include "double_tether.hpp"

#include <iostream>
#include <cmath>
#include <stdexcept>

using namespace BeadRodPmer;

/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
DoubleTether::DoubleTether(const std::vector<std::string> & v_line)
  : SingleTether(v_line)
{
  
  rhs_of_G.resize(Nbeads+5);
  dummy_for_noise.resize(Nbeads+5);
  Gmunu.resize(Nbeads+5,Nbeads+5);


  rhs_of_Hhat.resize(Nbeads+5);
  tension.resize(Nbeads+5);
  Hhat.resize(Nbeads+5,Nbeads+5);


  
  dCdlambda.resize(Nbeads+5,Nbeads+5);

  


  constraint_errors.resize(Nbeads+5);

  negative_tension_change.resize(Nbeads+5);
  
  tDets.resize(Nbeads+3);
  bDets.resize(Nbeads+3);

}


/* ============================================================================ */
/* Initial integration step. The force vector Fs needs to have been set (either
   by particle-particle interactions or set to zero). */
/* ============================================================================ */
void DoubleTether::first_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			      Eigen::Ref<Eigen::Matrix3Xd> Fs,double dt)
{
  set_unprojected_noise(dt);
  update_G();
  update_Hhat();
  compute_noise();
  compute_effective_kappa();
  compute_uc_forces(Fs);
  
  compute_tension({0,0,0},{0,0,0},Fs);
  initial_integrate(xs,Fs,dt,3,DOUBLE);

}

/* ============================================================================ */
/* Initial integration step. The force vector Fs needs to have been set (either
   by particle-particle interactions or set to zero). */
/* ============================================================================ */
void DoubleTether::first_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			      Eigen::Ref<Eigen::Matrix3Xd> Fs,double dt,
			      const Eigen::Ref<const Eigen::Vector3d> &dX0dt,
			      const Eigen::Ref<const Eigen::Vector3d> & dXNdt )
{
  set_unprojected_noise(dt);
  update_G();
  update_Hhat();
  compute_noise();
  compute_effective_kappa();
  compute_uc_forces(Fs);
  
  compute_tension(dX0dt,dXNdt,Fs);
  initial_integrate(xs,Fs,dt,3,DOUBLE);

}


/* ============================================================================ */
/* Second integration step. The force vector Fs needs to have been set (either
   by particle-particle interactions or set to zero). */
/* ============================================================================ */
int DoubleTether::second_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			      Eigen::Ref<Eigen::Matrix3Xd> Fs,
			      double dt,int itermax)
{
  update_Hhat();
  compute_effective_kappa();
  compute_uc_forces(Fs);
  
  compute_tension({0,0,0},{0,0,0},Fs);

  int iterations = correct_tension(xs,Fs,dt,xs.col(0),xs.col(Nbeads-1),itermax,1e-8);

  if (iterations > itermax) {
    for (int i = 0; i < get_Nbeads(); i++) 
      xs.col(i) = tmp_xs.col(i);
  } else {
    final_integrate(xs,Fs,dt,3,DOUBLE);
  }
  compute_tangents_and_friction(xs);
  
  return iterations;
}


/* ============================================================================ */
/* Second integration step. The force vector Fs needs to have been set (either
   by particle-particle interactions or set to zero). */
/* ============================================================================ */
int DoubleTether::second_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			      Eigen::Ref<Eigen::Matrix3Xd> Fs,
			      double dt,int itermax,
			      const Eigen::Ref<const Eigen::Vector3d> &X0_tFull,
			      const Eigen::Ref<const Eigen::Vector3d> &dX0dtHalf,
			      const Eigen::Ref<const Eigen::Vector3d> & XN_tFull,
			      const Eigen::Ref<const Eigen::Vector3d> & dXNdtHalf)
{
  update_Hhat();
  compute_effective_kappa();
  compute_uc_forces(Fs);
  
  compute_tension(dX0dtHalf,dXNdtHalf,Fs);

  int iterations = correct_tension(xs,Fs,dt,X0_tFull,XN_tFull,itermax,1e-8);


  if (iterations > itermax) {
    for (int i = 0; i < get_Nbeads(); i++) 
      xs.col(i) = tmp_xs.col(i);
  } else {
    final_integrate(xs,Fs,dt,3,DOUBLE);
  }
  compute_tangents_and_friction(xs);
  
  return iterations;
}



/* ---------------------------------------------------------------------------- */
/* RHS of G*eta = P. */
/* ---------------------------------------------------------------------------- */
void DoubleTether::set_rhs_of_G()
{

  int offset = 3;
  
  SingleTether::set_rhs_of_G();

  int ind = Nbeads-1+offset;
  
  rhs_of_G({ind,ind+1,ind+2}) = unprojected_noises.col(Nbeads-1);
  
  return;
  
}




/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix G. */
/* -------------------------------------------------------------------------- */
std::vector<T> DoubleTether::init_G_coeffsmatrix()
{
  // just set lower diagonal
  int offset = 3;
  std::vector<T> coeffs;

  coeffs = SingleTether::init_G_coeffsmatrix();


  Eigen::Vector3d btmp = bonds.col(Nbeads-2);

  coeffs.push_back(T(offset+Nbeads-1,offset+Nbeads-2,btmp(0)));
  coeffs.push_back(T(offset+Nbeads-1,offset+Nbeads-1,1));
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads-2,btmp(1)));
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads,1));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads-2,btmp(2)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads+1,1));
  
  return coeffs;

  
}

/* -------------------------------------------------------------------------- */
/* Update matrix G. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void DoubleTether::update_G()
{
  // update first four cols by hand

  int offset = 3;

  SingleTether::update_G();
  
  Gmunu.coeffRef(offset+Nbeads-1,offset+Nbeads-2) = bonds(0,Nbeads-2);
  Gmunu.coeffRef(offset+Nbeads,offset+Nbeads-2) = bonds(1,Nbeads-2);
  Gmunu.coeffRef(offset+Nbeads+1,offset+Nbeads-2) = bonds(2,Nbeads-2);
  
  return;
}

void DoubleTether::compute_effective_kappa()
{

  int offset = 3;


  set_bdets_and_tdets();

  
  tDets(Nbeads) = tDets(Nbeads-1) - bonds(0,Nbeads-2)*bonds(0,Nbeads-2)*tDets(Nbeads-2);
  tDets(Nbeads+1) = tDets(Nbeads) - bonds(1,Nbeads-2)*bonds(1,Nbeads-2)*tDets(Nbeads-2);
  tDets(Nbeads+2) = tDets(Nbeads+1) - bonds(2,Nbeads-2)*bonds(2,Nbeads-2)*tDets(Nbeads-2);

  double gDet = tDets(Nbeads+2);
  end_inverses({0,1,2}) = bonds.col(0)*bDets(1+offset)/(gDet*bondlength);


  end_inverses({3,4,5}) = -bonds.col(Nbeads-2)*tDets(Nbeads-2)/(gDet*bondlength);

  for (int i = 0; i < Nbeads-2; i++) {
    k_effs(i) = (kappa + temp*bondlength*costhetas(i)*tDets(i)*bDets(i+2+offset)/gDet
		 )/(bondlength*bondlength);
  }

  return;

}


/* ---------------------------------------------------------------------------- */
/* RHS of Hhat*lambda = Q. */
/* ---------------------------------------------------------------------------- */
void DoubleTether::set_rhs_of_Hhat(const Eigen::Vector3d &dXdt_at_1,
				   const Eigen::Vector3d &dXdt_at_N,
				   const Eigen::Ref<const Eigen::Matrix3Xd>& Fs)
{


  int offset = 3;

  
  

  rhs_of_Hhat({0,1,2}) = frictions[0]*(Fs.col(0)+noises.col(0))
    - dXdt_at_1;

  NoTether::set_rhs_of_Hhat(offset,Fs);

  int ind = Nbeads-1 + offset;
  rhs_of_Hhat({ind,ind+1,ind+2})
    = frictions[Nbeads-1]*(Fs.col(Nbeads-1)+noises.col(Nbeads-1))
    - dXdt_at_N;

  
  return;
  
}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix H hat. */
/* -------------------------------------------------------------------------- */
std::vector<T> DoubleTether::init_Hhat_coeffsmatrix()
{

  // only use lower triangular part as it's a symmetric matrix.  
  std::vector<T> coeffs;

  int offset = 3;
  
  coeffs = SingleTether::init_Hhat_coeffsmatrix();

  
  coeffs.push_back(T(offset+Nbeads-1,offset+Nbeads-2,Hhat_bottomside(0)));
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads-2,Hhat_bottomside(1)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads-2,Hhat_bottomside(2)));

  coeffs.push_back(T(offset+Nbeads-1,offset+Nbeads-1,Hhat_endblocks(0,0,Nbeads-1)));
  
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads-1,Hhat_endblocks(1,0,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads,Hhat_endblocks(1,1,Nbeads-1)));

  
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads-1,Hhat_endblocks(2,0,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads,Hhat_endblocks(2,1,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads+1,Hhat_endblocks(2,2,Nbeads-1)));

  return coeffs;

  
}






/* -------------------------------------------------------------------------- */
/* Update matrix M. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void DoubleTether::update_Hhat()
{

  int offset = 3;


  SingleTether::update_Hhat();
  
  Hhat.coeffRef(offset+Nbeads-1,offset+Nbeads-2) = Hhat_bottomside(0);
  Hhat.coeffRef(offset+Nbeads,offset+Nbeads-2) = Hhat_bottomside(1);
  Hhat.coeffRef(offset+Nbeads+1,offset+Nbeads-2) = Hhat_bottomside(2);

  Hhat.coeffRef(offset+Nbeads-1,offset+Nbeads-1) = Hhat_endblocks(0,0,Nbeads-1);

  Hhat.coeffRef(offset+Nbeads,offset+Nbeads-1) = Hhat_endblocks(1,0,Nbeads-1);
  Hhat.coeffRef(offset+Nbeads,offset+Nbeads) = Hhat_endblocks(1,1,Nbeads-1);

  Hhat.coeffRef(offset+Nbeads+1,offset+Nbeads-1) = Hhat_endblocks(2,0,Nbeads-1);
  Hhat.coeffRef(offset+Nbeads+1,offset+Nbeads) = Hhat_endblocks(2,1,Nbeads-1);
  Hhat.coeffRef(offset+Nbeads+1,offset+Nbeads+1) = Hhat_endblocks(2,2,Nbeads-1);

  return;
}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix dCdlambda. */
/* -------------------------------------------------------------------------- */
std::vector<T> DoubleTether::init_dCdlambda_coeffsmatrix()
{

  // initializing the vector as the full (both upper and lower) part of Hhat, since
  // this matrix won't be symmetric.
  std::vector<T> coeffs;

  int offset = 3;

  coeffs = SingleTether::init_dCdlambda_coeffsmatrix();
  
  coeffs.push_back(T(offset+Nbeads-2,offset+Nbeads-1,Hhat_bottomside(0)));
  coeffs.push_back(T(offset+Nbeads-2,offset+Nbeads,Hhat_bottomside(1)));
  coeffs.push_back(T(offset+Nbeads-2,offset+Nbeads+1,Hhat_bottomside(2)));

  
  coeffs.push_back(T(offset+Nbeads-1,offset+Nbeads-2,Hhat_bottomside(0)));
  coeffs.push_back(T(offset+Nbeads-1,offset+Nbeads-1,Hhat_endblocks(0,0,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads-1,offset+Nbeads,Hhat_endblocks(1,0,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads-1,offset+Nbeads+1,Hhat_endblocks(2,0,Nbeads-1)));
  
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads-2,Hhat_bottomside(1)));
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads-1,Hhat_endblocks(1,0,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads,Hhat_endblocks(1,1,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads+1,Hhat_endblocks(2,1,Nbeads-1)));
  
  
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads-2,Hhat_bottomside(2)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads-1,Hhat_endblocks(2,0,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads,Hhat_endblocks(2,1,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads+1,Hhat_endblocks(2,2,Nbeads-1)));

  return coeffs;


}






void DoubleTether::update_dCdlambda(double Delta_t)
{
  // must update column-wise
  int offset = 3;


  SingleTether::update_dCdlambda(Delta_t);

  // Nbeads - 1 col

  dCdlambda.coeffRef(offset+Nbeads-1,offset+Nbeads-2) = -Hhat_bottomside(0)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads,offset+Nbeads-2) = -Hhat_bottomside(1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads+1,offset+Nbeads-2) = -Hhat_bottomside(2)*Delta_t;
  
  
  
  // Nbeads col
  dCdlambda.coeffRef(offset+Nbeads-2,offset+Nbeads-1) = -dCdlambda_bottomside(0)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads-1,offset+Nbeads-1) = -Hhat_endblocks(0,0,Nbeads-1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads,offset+Nbeads-1) = -Hhat_endblocks(1,0,Nbeads-1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads+1,offset+Nbeads-1) = -Hhat_endblocks(2,0,Nbeads-1)*Delta_t;
  
  
  // Nbeads+1 col
  dCdlambda.coeffRef(offset+Nbeads-2,offset+Nbeads) = -dCdlambda_bottomside(1)*Delta_t;

  dCdlambda.coeffRef(offset+Nbeads-1,offset+Nbeads) = -Hhat_endblocks(1,0,Nbeads-1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads,offset+Nbeads) = -Hhat_endblocks(1,1,Nbeads-1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads+1,offset+Nbeads) = -Hhat_endblocks(2,1,Nbeads-1)*Delta_t;


  // Nbeads +2 col
  dCdlambda.coeffRef(offset+Nbeads-2,offset+Nbeads+1) = -dCdlambda_bottomside(2)*Delta_t;

  dCdlambda.coeffRef(offset+Nbeads-1,offset+Nbeads+1) = -Hhat_endblocks(2,0,Nbeads-1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads,offset+Nbeads+1) = -Hhat_endblocks(2,1,Nbeads-1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads+1,offset+Nbeads+1) = -Hhat_endblocks(2,2,Nbeads-1)*Delta_t;

  return;

}



void DoubleTether::compute_noise()
{

  int offset = 3;

  SingleTether::compute_noise();

  int ind = offset + Nbeads-1;
  noises.col(Nbeads-1) -= dummy_for_noise({ind,ind+1,ind+2});

  return;
  
}


void DoubleTether::compute_tension(const Eigen::Vector3d & dXdt_at_1,
				   const Eigen::Vector3d & dXdt_at_N,
				   const Eigen::Ref<const Eigen::Matrix3Xd>& Fs)
{



  set_rhs_of_Hhat(dXdt_at_1,dXdt_at_N,Fs);

  Hhat_solver->factorize(Hhat);
  tension =  Hhat_solver->solve(rhs_of_Hhat);
  
}
/*
void DoubleTether::test_jacob(int mu,double Delta_t,const Eigen::Vector3d & X_of_t_at_1,
			      const Eigen::Vector3d & X_of_t_at_N)
{

  int offset = 3;  
  // start by initializing the vector dC at whatever tensions


  
  calculate_constraint_errors(X_of_t_at_1,X_of_t_at_N);


  update_dCdlambda(Delta_t);

  
  //  std::cout << "exact calc: " << std::endl;
  //  std::cout << dCdlambda.col(mu+offset-1) << std::endl;

  Eigen::VectorXd zeroval_change = constraint_errors;

 
  double tens_change_old = 0;
  for (double tens_change = 0.00001; tens_change > 0.000009; tens_change /= 2) {

    tension(mu+offset-1) -= tens_change_old;
    tension(mu+offset-1) += tens_change;
    
    final_integrate(Delta_t,3,DOUBLE);
    calculate_constraint_errors(X_of_t_at_1,X_of_t_at_N);

    
    std::cout << "numerical estimate with dlamda =  " << tens_change
	      << "and lambda = " << tension(mu+offset-1) << " is = " << std::endl;
    std::cout << (constraint_errors-zeroval_change)/tens_change-dCdlambda.col(mu+offset-1)
	      << std::endl;

    tens_change_old = tens_change;

  }
  
  return;
}  
*/
int DoubleTether::correct_tension(Eigen::Ref<Eigen::Matrix3Xd> xs,
				  const Eigen::Ref<const Eigen::Matrix3Xd>&Fs,
				  double Delta_t,const Eigen::Vector3d & X_of_t_at_1,
				  const Eigen::Vector3d & X_of_t_at_N,
				  int itermax,double tolerance)
{

  // set C_mu and dC_mu/dlambda_nu
  final_integrate(xs,Fs,Delta_t,3,DOUBLE);
  calculate_constraint_errors(X_of_t_at_1,X_of_t_at_N,xs);
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



    final_integrate(xs,Fs,Delta_t,3,DOUBLE);
    calculate_constraint_errors(X_of_t_at_1,X_of_t_at_N,xs);
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




void DoubleTether::calculate_constraint_errors(const Eigen::Vector3d & X_of_t_at_1,
					       const Eigen::Vector3d & X_of_t_at_N,
					       const Eigen::Ref<const Eigen::Matrix3Xd>& xs)
{
  int offset = 3;

  SingleTether::calculate_constraint_errors(X_of_t_at_1,xs);

  int ind = offset+Nbeads-1;
  constraint_errors({ind,ind+1,ind+2}) = xs.col(Nbeads-1)-X_of_t_at_N;
  
  return;

}

