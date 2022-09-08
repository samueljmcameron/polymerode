#include "double_tether.hpp"
#include "input.hpp"
#include "initialise.hpp"

#include <iostream>
#include <cmath>
#include <stdexcept>

using namespace BeadRodPmer;

/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
DoubleTether::DoubleTether(const std::vector<std::string> & splitvec)
  : Polymer(splitvec)
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

  std::cout << "initiating double tether." << std::endl;
  
  Initialise::init_atoms(splitvec,atoms,x0,xN,initspringK,initdt,inittolerance,0);
  
}

/* -------------------------------------------------------------------------- */
/* Destructor */
/* -------------------------------------------------------------------------- */
DoubleTether::~DoubleTether()
{
}



int DoubleTether::single_step(double t, double dt,
			      const std::vector<std::vector<double>> &dFdX_i,
			      int itermax, int numtries,bool throw_exception)
{

  // get here if this step has been restarted numtry times
  if (numtries == 0) {
    if (throw_exception) 
      throw std::runtime_error("could not solve constraint equation for double tethered polymer.");
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

  
  compute_tension({0,0,0},{0,0,0});
  initial_integrate(dt);
  

  
  update_Hhat();
  compute_effective_kappa();
  compute_uc_forces();

  for (int index = 0; index < nuc_beads.size(); index ++ )
    add_external_force(dFdX_i[index],nuc_beads[index]);
  
  
  compute_tension({0,0,0},{0,0,0});
  
  int iterations = correct_tension(dt,atoms[0].R,
				   atoms[get_Nbeads()-1].R,
				   itermax,1e-8);



  if (iterations > itermax) {
    numtries -= 1;
    std::cout << "too many iterations when correcting tension at time " << t
	      <<  ", retrying the step with new noise ( " << numtries 
	      << " attempts left). " << std::endl;

    for (int i = 0; i < get_Nbeads(); i++) 
      atoms[i].R = Rtmp[i];

    compute_tangents_and_friction();
    
    return single_step(t,dt,dFdX_i,itermax,numtries,throw_exception);
  }  else {
    final_integrate(dt);
  
    compute_tangents_and_friction();
  }

  return 0;
  

}


int DoubleTether::single_step(double t,double dt,
			      const std::vector<std::vector<double>> & dFdX_i,
			      std::function<Eigen::Vector3d (double)> X0_t,
			      std::function<Eigen::Vector3d (double)> XN_t,
			      std::function<Eigen::Vector3d (double)> dX0dt,
			      std::function<Eigen::Vector3d (double)> dXNdt,
			      int itermax, int numtries,bool throw_exception)
{

  // get here if this step has been restarted numtry times
  if (numtries == 0) {
    if (throw_exception)
      throw std::runtime_error("could not solve constraint equation for double tethered polymer.");
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

  
  compute_tension(dX0dt(t+dt/2),dXNdt(t+dt/2));
  initial_integrate(dt);
  

  
  update_Hhat();
  compute_effective_kappa();
  compute_uc_forces();

  for (int index = 0; index < nuc_beads.size(); index ++ )
    add_external_force(dFdX_i[index],nuc_beads[index]);
  
  
  compute_tension(dX0dt(t+dt),dXNdt(t+dt));
  
  int iterations = correct_tension(dt,X0_t(t+dt),XN_t(t+dt),
				   itermax,1e-8);



  if (iterations > itermax) {
    numtries -= 1;
    std::cout << "too many iterations when correcting tension at time " << t
	      <<  ", retrying the step with new noise ( " << numtries 
	      << " attempts left). " << std::endl;
    for (int i = 0; i < get_Nbeads(); i++) 
      atoms[i].R = Rtmp[i];

    compute_tangents_and_friction();
    return single_step(t,dt,dFdX_i,X0_t,XN_t,dX0dt,dXNdt,itermax,numtries,throw_exception);
  }  else {
    final_integrate(dt);
  
    compute_tangents_and_friction();
  }

  return 0;
}


/* ---------------------------------------------------------------------------- */
/* RHS of G*eta = P. */
/* ---------------------------------------------------------------------------- */
void DoubleTether::set_rhs_of_G()
{
  int offset = 2;
  rhs_of_G(-2+offset) = atoms[0].unprojected_noise(0);
  rhs_of_G(-1+offset) = atoms[0].unprojected_noise(1);
  rhs_of_G(0+offset) = atoms[0].unprojected_noise(2);
  for (int i = 0; i < Nbeads-1; i++) {
    rhs_of_G(i+1+offset) = bonds[i].rod.dot(atoms[i+1].unprojected_noise
					  -atoms[i].unprojected_noise);

  }

  rhs_of_G(Nbeads+offset) = atoms[Nbeads-1].unprojected_noise(0);
  rhs_of_G(Nbeads+1+offset) = atoms[Nbeads-1].unprojected_noise(1);
  rhs_of_G(Nbeads+2+offset) = atoms[Nbeads-1].unprojected_noise(2);
  
  return;
  
}


/* -------------------------------------------------------------------------- */
/* Initialse G in G*eta = P (only call once). */
/* -------------------------------------------------------------------------- */
void DoubleTether::set_G()
{

  std::vector<T> coefficients = init_G_coeffsmatrix();

  Gmunu.setFromTriplets(coefficients.begin(),coefficients.end());

  Gmunu_solver.analyzePattern(Gmunu);

  return;

}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix G. */
/* -------------------------------------------------------------------------- */
std::vector<T> DoubleTether::init_G_coeffsmatrix()
{
  // just set lower diagonal
  int offset = 2;
  std::vector<T> coeffs;
  coeffs.push_back(T(offset-2,offset-2,1));

  coeffs.push_back(T(offset-1,offset-1,1));


  coeffs.push_back(T(offset,offset,1));
  
  coeffs.push_back(T(offset+1,offset-2,-bonds[0].rod(0)));
  coeffs.push_back(T(offset+1,offset-1,-bonds[0].rod(1)));
  coeffs.push_back(T(offset+1,offset,-bonds[0].rod(2)));

  coeffs.push_back(T(offset+1,offset+1,2));
  
  for (int i = 1; i < Nbeads-1; i++) {
    coeffs.push_back(T(i+offset+1,i+offset+1,2));
    coeffs.push_back(T(i+offset+1,i+offset,-costhetas(i-1)));

  }


  coeffs.push_back(T(offset+Nbeads,offset+Nbeads-1,bonds[Nbeads-2].rod(0)));
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads,1));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads-1,bonds[Nbeads-2].rod(1)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads+1,1));
  coeffs.push_back(T(offset+Nbeads+2,offset+Nbeads-1,bonds[Nbeads-2].rod(2)));
  coeffs.push_back(T(offset+Nbeads+2,offset+Nbeads+2,1));
  
  return coeffs;

  
}

/* -------------------------------------------------------------------------- */
/* Update matrix G. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void DoubleTether::update_G()
{
  // update first four cols by hand

  int offset = 2;

  Gmunu.coeffRef(offset+1,offset-2) = -bonds[0].rod(0);
  Gmunu.coeffRef(offset+1,offset-1) = -bonds[0].rod(1);
  Gmunu.coeffRef(offset+1,offset) = -bonds[0].rod(2);
  Gmunu.coeffRef(offset+2,offset+1) = -costhetas(0);


  // update middle columns
  for (int k = 1; k < Nbeads-2; k++) {


    int count = 0;
    for (SpMat::InnerIterator it(Gmunu,k+offset+1); it; ++it) {
      if (count == 1) {
	it.valueRef() = -costhetas(k);

      }
      count += 1;
    }
  }

  Gmunu.coeffRef(offset+Nbeads,offset+Nbeads-1) = bonds[Nbeads-2].rod(0);
  Gmunu.coeffRef(offset+Nbeads+1,offset+Nbeads-1) = bonds[Nbeads-2].rod(1);
  Gmunu.coeffRef(offset+Nbeads+2,offset+Nbeads-1) = bonds[Nbeads-2].rod(2);
  
  return;
}

void DoubleTether::compute_effective_kappa()
{

  int offset = 2;


  bDets(Nbeads+offset) = 1.0;
  bDets(Nbeads-1+offset) = 2 - bonds[Nbeads-2].rod.squaredNorm();
  tDets(0) = 1.0;
  tDets(1) = 2 - bonds[0].rod.squaredNorm();

  int mu;
  for (int i = 0; i < Nbeads - 2; i++ ) {

    mu = Nbeads-2-i;
    bDets(mu+offset) = 2*bDets(mu+1+offset)
      - costhetas(mu+1-2)*costhetas(mu+1-2)*bDets(mu+2+offset);

    mu = i+2;

    tDets(mu) = 2*tDets(mu-1)-costhetas(mu-2)*costhetas(mu-2)*tDets(mu-2);


  }

  bDets(offset) = bDets(1+offset)-bonds[0].rod(2)*bonds[0].rod(2)*bDets(2+offset);
  bDets(-1+offset) = bDets(offset)-bonds[0].rod(1)*bonds[0].rod(1)*bDets(2+offset);
  bDets(-2+offset) = bDets(-1+offset)-bonds[0].rod(0)*bonds[0].rod(0)*bDets(2+offset);

  tDets(Nbeads) = tDets(Nbeads-1) - bonds[Nbeads-2].rod(0)*bonds[Nbeads-2].rod(0)*tDets(Nbeads-2);
  tDets(Nbeads+1) = tDets(Nbeads) - bonds[Nbeads-2].rod(1)*bonds[Nbeads-2].rod(1)*tDets(Nbeads-2);
  tDets(Nbeads+2) = tDets(Nbeads+1) - bonds[Nbeads-2].rod(2)*bonds[Nbeads-2].rod(2)*tDets(Nbeads-2);

  double gDet = tDets(Nbeads+2);
  end_inverses(0) = bonds[0].rod(0)*bDets(2+offset)/(gDet*bondlength);
  end_inverses(1) = bonds[0].rod(1)*bDets(2+offset)/(gDet*bondlength);
  end_inverses(2) = bonds[0].rod(2)*bDets(2+offset)/(gDet*bondlength);


  end_inverses(3) = -bonds[Nbeads-2].rod(0)*tDets(Nbeads-2)/(gDet*bondlength);
  end_inverses(4) = -bonds[Nbeads-2].rod(1)*tDets(Nbeads-2)/(gDet*bondlength);
  end_inverses(5) = -bonds[Nbeads-2].rod(2)*tDets(Nbeads-2)/(gDet*bondlength);

  

  
  for (int i = 0; i < Nbeads-2; i++) {
    mu = i+2;
    k_effs(i) = (kappa - temp*bondlength*costhetas(i)*tDets(mu-2)*bDets(mu+offset+1)/gDet
		 )/(bondlength*bondlength);
  }


  
  
  return;

}


/* ---------------------------------------------------------------------------- */
/* RHS of Hhat*lambda = Q. */
/* ---------------------------------------------------------------------------- */
void DoubleTether::set_rhs_of_Hhat(const Eigen::Vector3d &dXdt_at_1,
				   const Eigen::Vector3d &dXdt_at_N)
{


  int offset = 2;

  
  Eigen::Vector3d tmp = atoms[0].friction*(atoms[0].Fpot+atoms[0].noise) - dXdt_at_1;

  rhs_of_Hhat(-2+offset) = tmp(0);
  rhs_of_Hhat(-1+offset) = tmp(1);
  rhs_of_Hhat(0+offset) = tmp(2);



  
  for (int i = 0; i< Nbeads-1; i++) {


    rhs_of_Hhat(i+offset+1) = bonds[i].rod.dot(atoms[i+1].friction*(atoms[i+1].Fpot
								    +atoms[i+1].noise)
					       -atoms[i].friction*(atoms[i].Fpot
								   +atoms[i].noise));

  }

  tmp = atoms[Nbeads-1].friction*(atoms[Nbeads-1].Fpot+atoms[Nbeads-1].noise) - dXdt_at_N;

  rhs_of_Hhat(offset+Nbeads) = tmp(0);
  rhs_of_Hhat(offset+Nbeads+1) = tmp(1);
  rhs_of_Hhat(offset+Nbeads+2) = tmp(2);


  
  return;
  
}






/* -------------------------------------------------------------------------- */
/* Initialise matrix H hat. */
/* -------------------------------------------------------------------------- */
void DoubleTether::set_Hhat()
{

  std::vector<T> coefficients = init_Hhat_coeffsmatrix();

  Hhat.setFromTriplets(coefficients.begin(),coefficients.end());


  Hhat_solver.analyzePattern(Hhat);

  return;

}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix H hat. */
/* -------------------------------------------------------------------------- */
std::vector<T> DoubleTether::init_Hhat_coeffsmatrix()
{

  // only use lower triangular part as it's a symmetric matrix.  
  std::vector<T> coeffs;

  int offset = 2;

  coeffs.push_back(T(offset-2,offset-2,Hhat_endblocks(0,0,0)));
  coeffs.push_back(T(offset-1,offset-2,Hhat_endblocks(1,0,0)));
  coeffs.push_back(T(offset-1,offset-1,Hhat_endblocks(1,1,0)));
  coeffs.push_back(T(offset,offset-2,Hhat_endblocks(2,0,0)));
  coeffs.push_back(T(offset,offset-1,Hhat_endblocks(2,1,0)));
  coeffs.push_back(T(offset,offset,Hhat_endblocks(2,2,0)));


  
  coeffs.push_back(T(offset+1,offset-2,Hhat_leftside(0)));
  coeffs.push_back(T(offset+1,offset-1,Hhat_leftside(1)));
  coeffs.push_back(T(offset+1,offset,Hhat_leftside(2)));

  
  coeffs.push_back(T(offset+1,offset+1,Hhat_diag_val(0)));
  
  
  for (int i = 1; i < Nbeads-1; i++) {
    coeffs.push_back(T(i+offset+1,i+offset+1,Hhat_diag_val(i)));
    coeffs.push_back(T(i+offset+1,i+offset,Hhat_loweroff_val(i)));

  }

  coeffs.push_back(T(offset+Nbeads,offset+Nbeads-1,Hhat_bottomside(0)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads-1,Hhat_bottomside(1)));
  coeffs.push_back(T(offset+Nbeads+2,offset+Nbeads-1,Hhat_bottomside(2)));

  coeffs.push_back(T(offset+Nbeads,offset+Nbeads,Hhat_endblocks(0,0,Nbeads-1)));
  
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads,Hhat_endblocks(1,0,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads+1,Hhat_endblocks(1,1,Nbeads-1)));

  
  coeffs.push_back(T(offset+Nbeads+2,offset+Nbeads,Hhat_endblocks(2,0,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads+2,offset+Nbeads+1,Hhat_endblocks(2,1,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads+2,offset+Nbeads+2,Hhat_endblocks(2,2,Nbeads-1)));

  return coeffs;

  
}






/* -------------------------------------------------------------------------- */
/* Update matrix M. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void DoubleTether::update_Hhat()
{

  int offset = 2;
  // update first four columns by hand

  Hhat.coeffRef(offset-2,offset-2) = Hhat_endblocks(0,0,0);
  Hhat.coeffRef(offset-1,offset-2) = Hhat_endblocks(1,0,0);
  Hhat.coeffRef(offset-1,offset-1) = Hhat_endblocks(1,1,0);
  Hhat.coeffRef(offset,offset-2) = Hhat_endblocks(2,0,0);
  Hhat.coeffRef(offset,offset-1) = Hhat_endblocks(2,1,0);
  Hhat.coeffRef(offset,offset) = Hhat_endblocks(2,2,0);


  Hhat.coeffRef(offset+1,offset-2) = Hhat_leftside(0);
  Hhat.coeffRef(offset+1,offset-1) = Hhat_leftside(1);
  Hhat.coeffRef(offset+1,offset) = Hhat_leftside(2);
    
  Hhat.coeffRef(offset+1,offset+1) = Hhat_diag_val(0);

  Hhat.coeffRef(offset+2,offset+1) = Hhat_loweroff_val(1);

  // update middle columns
  for (int k = 1; k < Nbeads-2; k++) {

    int count = 0;
    for (SpMat::InnerIterator it(Hhat,k+offset+1); it; ++it) {
      if (count == 0) {
	it.valueRef() = Hhat_diag_val(k);

      } else {
	it.valueRef() = Hhat_loweroff_val(k+1);
      }

      count += 1;
    }
  }

  Hhat.coeffRef(offset+Nbeads-1,offset+Nbeads-1) = Hhat_diag_val(Nbeads-2);

  Hhat.coeffRef(offset+Nbeads,offset+Nbeads-1) = Hhat_bottomside(0);
  Hhat.coeffRef(offset+Nbeads+1,offset+Nbeads-1) = Hhat_bottomside(1);
  Hhat.coeffRef(offset+Nbeads+2,offset+Nbeads-1) = Hhat_bottomside(2);

  Hhat.coeffRef(offset+Nbeads,offset+Nbeads) = Hhat_endblocks(0,0,Nbeads-1);

  Hhat.coeffRef(offset+Nbeads+1,offset+Nbeads) = Hhat_endblocks(1,0,Nbeads-1);
  Hhat.coeffRef(offset+Nbeads+1,offset+Nbeads+1) = Hhat_endblocks(1,1,Nbeads-1);

  Hhat.coeffRef(offset+Nbeads+2,offset+Nbeads) = Hhat_endblocks(2,0,Nbeads-1);
  Hhat.coeffRef(offset+Nbeads+2,offset+Nbeads+1) = Hhat_endblocks(2,1,Nbeads-1);
  Hhat.coeffRef(offset+Nbeads+2,offset+Nbeads+2) = Hhat_endblocks(2,2,Nbeads-1);

  return;
}




/* -------------------------------------------------------------------------- */
/* Initialise the dCdlambda matrix. */
/* -------------------------------------------------------------------------- */
void DoubleTether::set_dCdlambda()
{

  std::vector<T> coefficients = init_dCdlambda_coeffsmatrix();

  dCdlambda.setFromTriplets(coefficients.begin(),coefficients.end());

  jacob_solver.analyzePattern(dCdlambda);

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

  int offset = 2;

  coeffs.push_back(T(offset-2,offset-2,Hhat_endblocks(0,0,0)));
  coeffs.push_back(T(offset-2,offset-1,Hhat_endblocks(1,0,0)));
  coeffs.push_back(T(offset-2,offset,Hhat_endblocks(2,0,0)));
  coeffs.push_back(T(offset-2,offset+1,Hhat_leftside(0)));
  
  
  coeffs.push_back(T(offset-1,offset-2,Hhat_endblocks(1,0,0)));
  coeffs.push_back(T(offset-1,offset-1,Hhat_endblocks(1,1,0)));
  coeffs.push_back(T(offset-1,offset,Hhat_endblocks(2,1,0)));
  coeffs.push_back(T(offset-1,offset+1,Hhat_leftside(1)));
  
  coeffs.push_back(T(offset,offset-2,Hhat_endblocks(2,0,0)));
  coeffs.push_back(T(offset,offset-1,Hhat_endblocks(2,1,0)));
  coeffs.push_back(T(offset,offset,Hhat_endblocks(2,2,0)));
  coeffs.push_back(T(offset,offset+1,Hhat_leftside(2)));

  
  coeffs.push_back(T(offset+1,offset-2,Hhat_leftside(0)));
  coeffs.push_back(T(offset+1,offset-1,Hhat_leftside(1)));
  coeffs.push_back(T(offset+1,offset,Hhat_leftside(2)));
  coeffs.push_back(T(offset+1,offset+1,Hhat_diag_val(0)));
  
  
  for (int i = 1; i < Nbeads-1; i++) {
    coeffs.push_back(T(i+offset,i+offset+1,Hhat_loweroff_val(i)));
    coeffs.push_back(T(i+offset+1,i+offset+1,Hhat_diag_val(i)));
    coeffs.push_back(T(i+offset+1,i+offset,Hhat_loweroff_val(i)));


  }

  
  coeffs.push_back(T(offset+Nbeads-1,offset+Nbeads,Hhat_bottomside(0)));
  coeffs.push_back(T(offset+Nbeads-1,offset+Nbeads+1,Hhat_bottomside(1)));
  coeffs.push_back(T(offset+Nbeads-1,offset+Nbeads+2,Hhat_bottomside(2)));

  
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads-1,Hhat_bottomside(0)));
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads,Hhat_endblocks(0,0,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads+1,Hhat_endblocks(1,0,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads,offset+Nbeads+2,Hhat_endblocks(2,0,Nbeads-1)));
  
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads-1,Hhat_bottomside(1)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads,Hhat_endblocks(1,0,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads+1,Hhat_endblocks(1,1,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads+1,offset+Nbeads+2,Hhat_endblocks(2,1,Nbeads-1)));
  
  
  coeffs.push_back(T(offset+Nbeads+2,offset+Nbeads-1,Hhat_bottomside(2)));
  coeffs.push_back(T(offset+Nbeads+2,offset+Nbeads,Hhat_endblocks(2,0,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads+2,offset+Nbeads+1,Hhat_endblocks(2,1,Nbeads-1)));
  coeffs.push_back(T(offset+Nbeads+2,offset+Nbeads+2,Hhat_endblocks(2,2,Nbeads-1)));

  return coeffs;


}






void DoubleTether::update_dCdlambda(double Delta_t)
{
  // must update column-wise
  int offset = 2;


  // negative two col
  dCdlambda.coeffRef(offset-2,offset-2) = -Hhat_endblocks(0,0,0)*Delta_t;
  dCdlambda.coeffRef(offset-1,offset-2) = -Hhat_endblocks(1,0,0)*Delta_t;
  dCdlambda.coeffRef(offset,offset-2) = -Hhat_endblocks(2,0,0)*Delta_t;

  
  dCdlambda.coeffRef(offset+1,offset-2) = -dCdlambda_leftside(0)*Delta_t;


  // n one col
  dCdlambda.coeffRef(offset-2,offset-1) = -Hhat_endblocks(1,0,0)*Delta_t;
  dCdlambda.coeffRef(offset-1,offset-1) = -Hhat_endblocks(1,1,0)*Delta_t;
  dCdlambda.coeffRef(offset,offset-1) = -Hhat_endblocks(2,1,0)*Delta_t;
  dCdlambda.coeffRef(offset+1,offset-1) = -dCdlambda_leftside(1)*Delta_t;


  
  // zero col
  dCdlambda.coeffRef(offset-2,offset) = -Hhat_endblocks(2,0,0)*Delta_t;
  dCdlambda.coeffRef(offset-1,offset) = -Hhat_endblocks(2,1,0)*Delta_t;
  dCdlambda.coeffRef(offset,offset) = -Hhat_endblocks(2,2,0)*Delta_t;
  dCdlambda.coeffRef(offset+1,offset) = -dCdlambda_leftside(2)*Delta_t;



  // one col
  dCdlambda.coeffRef(offset-2,offset+1) = -Hhat_leftside(0)*Delta_t;
  dCdlambda.coeffRef(offset-1,offset+1) = -Hhat_leftside(1)*Delta_t;
  dCdlambda.coeffRef(offset,offset+1) = -Hhat_leftside(2)*Delta_t;
  dCdlambda.coeffRef(offset+1,offset+1) = -dCdlambda_diag_val(0)*Delta_t;
  dCdlambda.coeffRef(offset+2,offset+1)  = -dCdlambda_loweroff_val(1)*Delta_t;
  

  // update middle columns
  for (int k = 1; k < Nbeads-2; k++) {


    int count = 0;
    for (SpMat::InnerIterator it(dCdlambda,k+offset+1); it; ++it) {

      if (count == 0) {
    	it.valueRef() = -dCdlambda_upperoff_val(k)*Delta_t;
	
      } else if (count == 1) {
	
	it.valueRef() = -dCdlambda_diag_val(k)*Delta_t;
	
      } else {

	it.valueRef() = -dCdlambda_loweroff_val(k+1)*Delta_t;
      }

      count += 1;
    }
  }


  // Nbeads - 1 col
  dCdlambda.coeffRef(offset+Nbeads-2, offset+Nbeads-1) = - dCdlambda_upperoff_val(Nbeads-2)*Delta_t;
  							 
  dCdlambda.coeffRef(offset+Nbeads-1, offset+Nbeads-1) = - dCdlambda_diag_val(Nbeads-2)*Delta_t;

  dCdlambda.coeffRef(offset+Nbeads,offset+Nbeads-1) = -Hhat_bottomside(0)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads+1,offset+Nbeads-1) = -Hhat_bottomside(1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads+2,offset+Nbeads-1) = -Hhat_bottomside(2)*Delta_t;
  
  
  
  // Nbeads col
  dCdlambda.coeffRef(offset+Nbeads-1,offset+Nbeads) = -dCdlambda_bottomside(0)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads,offset+Nbeads) = -Hhat_endblocks(0,0,Nbeads-1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads+1,offset+Nbeads) = -Hhat_endblocks(1,0,Nbeads-1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads+2,offset+Nbeads) = -Hhat_endblocks(2,0,Nbeads-1)*Delta_t;
  
  
  // Nbeads+1 col
  dCdlambda.coeffRef(offset+Nbeads-1,offset+Nbeads+1) = -dCdlambda_bottomside(1)*Delta_t;

  dCdlambda.coeffRef(offset+Nbeads,offset+Nbeads+1) = -Hhat_endblocks(1,0,Nbeads-1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads+1,offset+Nbeads+1) = -Hhat_endblocks(1,1,Nbeads-1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads+2,offset+Nbeads+1) = -Hhat_endblocks(2,1,Nbeads-1)*Delta_t;


  // Nbeads +2 col
  dCdlambda.coeffRef(offset+Nbeads-1,offset+Nbeads+2) = -dCdlambda_bottomside(2)*Delta_t;

  dCdlambda.coeffRef(offset+Nbeads,offset+Nbeads+2) = -Hhat_endblocks(2,0,Nbeads-1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads+1,offset+Nbeads+2) = -Hhat_endblocks(2,1,Nbeads-1)*Delta_t;
  dCdlambda.coeffRef(offset+Nbeads+2,offset+Nbeads+2) = -Hhat_endblocks(2,2,Nbeads-1)*Delta_t;

  return;

}



void DoubleTether::compute_noise()
{

  int offset = 2;
  set_rhs_of_G();
  Gmunu_solver.factorize(Gmunu);

  dummy_for_noise =  Gmunu_solver.solve(rhs_of_G);

  int i = 0;
  atoms[i].noise = (atoms[i].unprojected_noise
		    + dummy_for_noise(i+offset+1)*bonds[i].rod);

  atoms[i].noise(0) -= dummy_for_noise(-2+offset);
  atoms[i].noise(1) -= dummy_for_noise(-1+offset);
  atoms[i].noise(2) -= dummy_for_noise(offset);
  
  
  for (i = 1; i < Nbeads-1; i++) {
    atoms[i].noise = (atoms[i].unprojected_noise
		      + dummy_for_noise(i+offset+1)*bonds[i].rod
		      - dummy_for_noise(i+offset)*bonds[i-1].rod);
  }


  i = Nbeads-1;

  atoms[i].noise = (atoms[i].unprojected_noise
		    - dummy_for_noise(i+offset)*bonds[i-1].rod);

  atoms[i].noise(0) -= dummy_for_noise(offset+1+Nbeads-1);
  atoms[i].noise(1) -= dummy_for_noise(offset+1+Nbeads);
  atoms[i].noise(2) -= dummy_for_noise(offset+1+Nbeads+1);
  
}


void DoubleTether::compute_tension(const Eigen::Vector3d & dXdt_at_1,
				   const Eigen::Vector3d & dXdt_at_N)
{



  set_rhs_of_Hhat(dXdt_at_1,dXdt_at_N);

  Hhat_solver.factorize(Hhat);
  tension =  Hhat_solver.solve(rhs_of_Hhat);
  
}

void DoubleTether::test_jacob(int mu,double Delta_t,const Eigen::Vector3d & X_of_t_at_1,
			      const Eigen::Vector3d & X_of_t_at_N)
{

  int offset = 2;  
  // start by initializing the vector dC at whatever tensions


  
  calculate_constraint_errors(X_of_t_at_1,X_of_t_at_N);


  update_dCdlambda(Delta_t);

  
  //  std::cout << "exact calc: " << std::endl;
  //  std::cout << dCdlambda.col(mu+offset) << std::endl;

  Eigen::VectorXd zeroval_change = constraint_errors;

 
  double tens_change_old = 0;
  for (double tens_change = 0.00001; tens_change > 0.000009; tens_change /= 2) {

    tension(mu+offset) -= tens_change_old;
    tension(mu+offset) += tens_change;
    
    final_integrate(Delta_t);
    calculate_constraint_errors(X_of_t_at_1,X_of_t_at_N);

    
    std::cout << "numerical estimate with dlamda =  " << tens_change
	      << "and lambda = " << tension(mu+offset) << " is = " << std::endl;
    std::cout << (constraint_errors-zeroval_change)/tens_change-dCdlambda.col(mu+offset)
	      << std::endl;

    tens_change_old = tens_change;

  }
  
  return;
}  

int DoubleTether::correct_tension(double Delta_t,const Eigen::Vector3d & X_of_t_at_1,
				  const Eigen::Vector3d & X_of_t_at_N,
				  int itermax,double tolerance)
{

  // set C_mu and dC_mu/dlambda_nu
  final_integrate(Delta_t);
  calculate_constraint_errors(X_of_t_at_1,X_of_t_at_N);
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



    final_integrate(Delta_t);
    calculate_constraint_errors(X_of_t_at_1,X_of_t_at_N);
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



/*----------------------------------------------------------------------------*/
/* Computes the bead positions at the midstep, along with updating the
   un-normalised bond tangents, the (normalised) bead tangents, friction
   tensors, and cos(theta_i)s. */
/*----------------------------------------------------------------------------*/
void DoubleTether::initial_integrate(double Delta_t)
{
  
  double tmp = Delta_t/2.0;
  int offset = 2;

  double tangentnorm;

  
  int i = 0;
  Rtmp[0] = atoms[0].R;

  atoms[0].t_force = tension(offset+1)*bonds[0].rod;

  atoms[0].t_force(0) -= tension(offset-2);
  atoms[0].t_force(1) -= tension(offset-1);
  atoms[0].t_force(2) -= tension(offset);

  atoms[0].R += tmp*atoms[0].friction*(atoms[0].Fpot+atoms[0].noise+atoms[0].t_force);

  i = 1;
  Rtmp[i] = atoms[i].R;


  atoms[i].t_force = tension(offset+1+i)*bonds[i].rod-tension(offset+i)*bonds[i-1].rod;
  
  
  atoms[i].R += tmp*atoms[i].friction*(atoms[i].Fpot+atoms[i].noise+atoms[i].t_force);

  bonds[i-1].rod = (atoms[i].R-atoms[i-1].R)/bondlength;
  tmpbonds[i-1].rod = bonds[i-1].rod;
  
  tangentnorm = sqrt(bonds[i-1].rod.dot(bonds[i-1].rod));
  atoms[i-1].tangent = bonds[i-1].rod/tangentnorm;


  // compute friction with new values of R on bead i-1
  single_inv_friction(i-1);
  
  for (i = 2; i < Nbeads-1; i++) {

    Rtmp[i] = atoms[i].R;

    atoms[i].t_force = tension(offset+1+i)*bonds[i].rod-tension(offset+i)*bonds[i-1].rod;


    atoms[i].R += tmp*atoms[i].friction*(atoms[i].Fpot+atoms[i].noise+atoms[i].t_force);

    bonds[i-1].rod = (atoms[i].R-atoms[i-1].R)/bondlength;
    tmpbonds[i-1].rod = bonds[i-1].rod;

    costhetas(i-2) = bonds[i-1].rod.dot(bonds[i-2].rod)
      /(bonds[i-1].rod.norm()*bonds[i-2].rod.norm());


    tangentnorm = sqrt((bonds[i-1].rod+bonds[i-2].rod
			).dot(bonds[i-1].rod+bonds[i-2].rod));


    atoms[i-1].tangent = (bonds[i-1].rod+bonds[i-2].rod)/tangentnorm;

    single_inv_friction(i-1);

  }


  i = Nbeads -1;

  Rtmp[i] = atoms[i].R;

  atoms[i].t_force = -tension(offset+i)*bonds[i-1].rod;

  atoms[i].t_force(0) -= tension(offset+Nbeads);
  atoms[i].t_force(1) -= tension(offset+Nbeads+1);
  atoms[i].t_force(2) -= tension(offset+Nbeads+2);


  atoms[i].R += tmp*atoms[i].friction*(atoms[i].Fpot+atoms[i].noise+atoms[i].t_force);

  bonds[i-1].rod = (atoms[i].R-atoms[i-1].R)/bondlength;
  tmpbonds[i-1].rod = bonds[i-1].rod;


  costhetas(i-2) = bonds[i-1].rod.dot(bonds[i-2].rod)
    /(bonds[i-1].rod.norm()*bonds[i-2].rod.norm());

  
  tangentnorm = sqrt((bonds[i-1].rod+bonds[i-2].rod
		      ).dot(bonds[i-1].rod+bonds[i-2].rod));
  
  
  atoms[i-1].tangent = (bonds[i-1].rod+bonds[i-2].rod)/tangentnorm;
  
  single_inv_friction(i-1);
  
  tangentnorm = sqrt(bonds[i-1].rod.dot(bonds[i-1].rod));
  atoms[i].tangent = bonds[i-1].rod/tangentnorm;

  single_inv_friction(i);

  return;
}

/*----------------------------------------------------------------------------*/
/* Computes the final bead positions and unnormalised bond tangents. */
/*----------------------------------------------------------------------------*/
void DoubleTether::final_integrate(double Delta_t)
{
  
  double tmp = Delta_t;

  int offset=2;


  atoms[0].t_force= tension(offset+1)*bonds[0].rod;

  atoms[0].t_force(0) -= tension(offset-2);
  atoms[0].t_force(1) -= tension(offset-1);
  atoms[0].t_force(2) -= tension(offset);
  

  atoms[0].R = Rtmp[0] + tmp*atoms[0].friction*(atoms[0].Fpot+atoms[0].noise+atoms[0].t_force);


    
  for (int i = 1; i < Nbeads-1; i++) {

    atoms[i].t_force = tension(offset+1+i)*bonds[i].rod-tension(offset+i)*bonds[i-1].rod;

    atoms[i].R = Rtmp[i] + tmp*atoms[i].friction*(atoms[i].Fpot+atoms[i].noise+atoms[i].t_force);
    tmpbonds[i-1].rod = (atoms[i].R-atoms[i-1].R)/bondlength;
    
  }

  int i = Nbeads -1;

  atoms[i].t_force = -tension(offset+i)*bonds[i-1].rod;

  atoms[i].t_force(0) -= tension(offset+Nbeads);
  atoms[i].t_force(1) -= tension(offset+Nbeads+1);
  atoms[i].t_force(2) -= tension(offset+Nbeads+2);

  atoms[i].R = Rtmp[i] + tmp*atoms[i].friction*(atoms[i].Fpot+atoms[i].noise+atoms[i].t_force);
  tmpbonds[i-1].rod = (atoms[i].R-atoms[i-1].R)/bondlength;
  return;
}



void DoubleTether::calculate_constraint_errors(const Eigen::Vector3d & X_of_t_at_1,
					       const Eigen::Vector3d & X_of_t_at_N)
{
  int offset = 2;
  constraint_errors(offset-2) = atoms[0].R(0)-X_of_t_at_1(0);
  constraint_errors(offset-1) = atoms[0].R(1)-X_of_t_at_1(1);
  constraint_errors(offset) = atoms[0].R(2)-X_of_t_at_1(2);
  for (int mu = 1; mu < Nbeads; mu++) {
    constraint_errors(offset+mu) = (atoms[mu].R-atoms[mu-1].R).norm()-bondlength;
  }

  constraint_errors(offset+Nbeads) = atoms[Nbeads-1].R(0)-X_of_t_at_N(0);
  constraint_errors(offset+Nbeads+1) = atoms[Nbeads-1].R(1)-X_of_t_at_N(1);
  constraint_errors(offset+Nbeads+2) = atoms[Nbeads-1].R(2)-X_of_t_at_N(2);
  return;

}

