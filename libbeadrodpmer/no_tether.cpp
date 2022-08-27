#include "no_tether.hpp"
#include "input.hpp"

#include <iostream>
#include <cmath>
#include <stdexcept>

#define SMALL 1e-14

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

  end_inverses.setZero();

  init_atoms();

}

/* -------------------------------------------------------------------------- */
/* Destructor */
/* -------------------------------------------------------------------------- */
NoTether::~NoTether()
{
}





void NoTether::single_step(double t, double dt,
			   const std::vector<std::vector<double>> & dFdX_i,
			   int itermax, int numtries)
{

  // get here if this step has been restarted numtry times
  if (numtries == 0)
    throw std::runtime_error("could not solve constraint equation for no tether polymer.");

  
  set_unprojected_noise(dt);
  update_G();
  update_Hhat();
  compute_noise();
  compute_effective_kappa();
  compute_uc_forces();

  for (int index = 0; index < nuc_beads.size(); index ++ )
    add_external_force(dFdX_i[index],nuc_beads[index]);
  
  compute_tension();
  initial_integrate(dt);
  

  
  update_Hhat();
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
      atoms[i].R = Rtmp[i];


    single_step(t,dt,dFdX_i,itermax,numtries);
  } else {
    final_integrate(dt);
  
    compute_tangents_and_friction();
  }

  return;
}




/* ---------------------------------------------------------------------------- */
/* RHS of G*eta = P. */
/* ---------------------------------------------------------------------------- */
void NoTether::set_rhs_of_G()
{


  
  for (int i = 0; i < Nbeads-1; i++) {
    rhs_of_G(i) = bonds[i].rod.dot(atoms[i+1].unprojected_noise
				   -atoms[i].unprojected_noise);

  }

  
  return;
  
}


/* -------------------------------------------------------------------------- */
/* Initialse G in G*eta = P (only call once). */
/* -------------------------------------------------------------------------- */
void NoTether::set_G()
{

  std::vector<T> coefficients = init_G_coeffsmatrix();

  Gmunu.setFromTriplets(coefficients.begin(),coefficients.end());

  Gmunu_solver.analyzePattern(Gmunu);

  return;

}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix G. */
/* -------------------------------------------------------------------------- */
std::vector<T> NoTether::init_G_coeffsmatrix()
{
  // just set lower diagonal

  std::vector<T> coeffs;

  coeffs.push_back(T(0,0,2));
  
  for (int i = 1; i < Nbeads-1; i++) {
    coeffs.push_back(T(i,i,2));
    coeffs.push_back(T(i,i-1,-costhetas(i-1)));

  }
  
  return coeffs;

  
}

/* -------------------------------------------------------------------------- */
/* Update matrix G. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void NoTether::update_G()
{
  // update first column by hand



  Gmunu.coeffRef(1,0) = -costhetas(0);


  // update middle columns
  for (int k = 1; k < Nbeads-2; k++) {


    int count = 0;
    for (SpMat::InnerIterator it(Gmunu,k); it; ++it) {
      if (count == 1) {
	it.valueRef() = -costhetas(k);

      }
      count += 1;
    }
  }

  
  return;
}

void NoTether::compute_effective_kappa()
{


  bDets(Nbeads-1) = 1.0;
  tDets(0) = 1.0;

  bDets(Nbeads-2) = 2;
  tDets(1) = 2;

  int mu;
  for (int i = 0; i < Nbeads - 2; i++ ) {

    mu = Nbeads-2-i;
    bDets(mu-1) = 2*bDets(mu)
      - costhetas(mu-1)*costhetas(mu-1)*bDets(mu+1);

    mu = i+2;

    tDets(mu) = 2*tDets(mu-1)-costhetas(mu-2)*costhetas(mu-2)*tDets(mu-2);


  }



  double gDet = tDets(Nbeads-1);
  
  for (int i = 0; i < Nbeads-2; i++) {

    k_effs(i) = (kappa - temp*bondlength*costhetas(i)*tDets(i)*bDets(i+2)/gDet
		 )/(bondlength*bondlength);
  }


  
  
  return;

}




/* ---------------------------------------------------------------------------- */
/* RHS of Hhat*lambda = Q. */
/* ---------------------------------------------------------------------------- */
void NoTether::set_rhs_of_Hhat()
{

  
  
  for (int i = 0; i< Nbeads-1; i++) {


    rhs_of_Hhat(i) = bonds[i].rod.dot(atoms[i+1].friction*(atoms[i+1].Fpot
							   +atoms[i+1].noise)
				      -atoms[i].friction*(atoms[i].Fpot
							  +atoms[i].noise));
    
  }

  return;
  
}






/* -------------------------------------------------------------------------- */
/* Initialise matrix H hat. */
/* -------------------------------------------------------------------------- */
void NoTether::set_Hhat()
{

  std::vector<T> coefficients = init_Hhat_coeffsmatrix();

  Hhat.setFromTriplets(coefficients.begin(),coefficients.end());


  Hhat_solver.analyzePattern(Hhat);

  return;

}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix H hat. */
/* -------------------------------------------------------------------------- */
std::vector<T> NoTether::init_Hhat_coeffsmatrix()
{

  // only use lower triangular part as it's a symmetric matrix.  
  std::vector<T> coeffs;



  
  coeffs.push_back(T(0,0,Hhat_diag_val(0)));
  
  
  for (int i = 1; i < Nbeads-1; i++) {
    coeffs.push_back(T(i,i,Hhat_diag_val(i)));
    coeffs.push_back(T(i,i-1,Hhat_loweroff_val(i)));

  }


  return coeffs;

  
}






/* -------------------------------------------------------------------------- */
/* Update matrix M. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void NoTether::update_Hhat()
{


  // update first column by hand

    
  Hhat.coeffRef(0,0) = Hhat_diag_val(0);

  Hhat.coeffRef(1,0) = Hhat_loweroff_val(1);

  // update middle columns
  for (int k = 1; k < Nbeads-2; k++) {

    int count = 0;
    for (SpMat::InnerIterator it(Hhat,k); it; ++it) {
      if (count == 0) {
	it.valueRef() = Hhat_diag_val(k);

      } else {
	it.valueRef() = Hhat_loweroff_val(k+1);
      }

      count += 1;
    }
  }
  Hhat.coeffRef(Nbeads-2,Nbeads-2) = Hhat_diag_val(Nbeads-2);
  return;
}




/* -------------------------------------------------------------------------- */
/* Initialise both the dCdlambda and mat_tmp matrices. */
/* -------------------------------------------------------------------------- */
void NoTether::set_dCdlambda()
{

  std::vector<T> coefficients = init_dCdlambda_coeffsmatrix();

  dCdlambda.setFromTriplets(coefficients.begin(),coefficients.end());

  jacob_solver.analyzePattern(dCdlambda);

  return;

}



/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix dCdlambda. */
/* -------------------------------------------------------------------------- */
std::vector<T> NoTether::init_dCdlambda_coeffsmatrix()
{

  // initializing the vector as the full (both upper and lower) part of Hhat, since
  // this matrix won't be symmetric.
  std::vector<T> coeffs;

  coeffs.push_back(T(0,0,Hhat_diag_val(0)));
  
  
  for (int i = 1; i < Nbeads-1; i++) {
    coeffs.push_back(T(i-1,i,Hhat_loweroff_val(i)));
    coeffs.push_back(T(i,i,Hhat_diag_val(i)));
    coeffs.push_back(T(i,i-1,Hhat_loweroff_val(i)));


  }

  

  return coeffs;


}






void NoTether::update_dCdlambda(double Delta_t)
{



  dCdlambda.coeffRef(0,0) = -dCdlambda_diag_val(0)*Delta_t;
  dCdlambda.coeffRef(1,0)  = -dCdlambda_loweroff_val(1)*Delta_t;
  

  // update middle columns
  for (int k = 1; k < Nbeads-2; k++) {


    int count = 0;
    for (SpMat::InnerIterator it(dCdlambda,k); it; ++it) {

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
  dCdlambda.coeffRef(Nbeads-3, Nbeads-2) = - dCdlambda_upperoff_val(Nbeads-2)*Delta_t;
  							 
  dCdlambda.coeffRef(Nbeads-2, Nbeads-2) = - dCdlambda_diag_val(Nbeads-2)*Delta_t;


  return;

}



void NoTether::compute_noise()
{


  set_rhs_of_G();
  Gmunu_solver.factorize(Gmunu);

  dummy_for_noise =  Gmunu_solver.solve(rhs_of_G);

  int i = 0;
  atoms[i].noise = (atoms[i].unprojected_noise
		    + dummy_for_noise(i)*bonds[i].rod);

  
  
  for (i = 1; i < Nbeads-1; i++) {
    atoms[i].noise = (atoms[i].unprojected_noise
		      + dummy_for_noise(i)*bonds[i].rod
		      - dummy_for_noise(i-1)*bonds[i-1].rod);
  }


  i = Nbeads-1;

  atoms[i].noise = (atoms[i].unprojected_noise
		    - dummy_for_noise(i-1)*bonds[i-1].rod);

  
}


void NoTether::compute_tension()
{



  set_rhs_of_Hhat();

  Hhat_solver.factorize(Hhat);
  tension =  Hhat_solver.solve(rhs_of_Hhat);
  
}

void NoTether::test_jacob(int mu,double Delta_t)
{


  // start by initializing the vector dC at whatever tensions
  
  
  calculate_constraint_errors();


  update_dCdlambda(Delta_t);

  
  //  std::cout << "exact calc: " << std::endl;
  //  std::cout << dCdlambda.col(mu) << std::endl;

  Eigen::VectorXd zeroval_change = constraint_errors;

 
  double tens_change_old = 0;
  for (double tens_change = 0.00001; tens_change > 0.000009; tens_change /= 2) {

    tension(mu) -= tens_change_old;
    tension(mu) += tens_change;
    
    final_integrate(Delta_t);

    calculate_constraint_errors();
    
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
  final_integrate(Delta_t);
  calculate_constraint_errors();
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
    calculate_constraint_errors();
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
void NoTether::initial_integrate(double Delta_t)
{
  
  double tmp = Delta_t/2.0;

  double tangentnorm;

  
  int i = 0;
  Rtmp[0] = atoms[0].R;

  atoms[0].t_force = tension(0)*bonds[0].rod;


  atoms[0].R += tmp*atoms[0].friction*(atoms[0].Fpot+atoms[0].noise+atoms[0].t_force);

  i = 1;
  Rtmp[i] = atoms[i].R;


  atoms[i].t_force = tension(i)*bonds[i].rod-tension(i-1)*bonds[i-1].rod;
  
  
  atoms[i].R += tmp*atoms[i].friction*(atoms[i].Fpot+atoms[i].noise+atoms[i].t_force);

  bonds[i-1].rod = (atoms[i].R-atoms[i-1].R)/bondlength;
  tmpbonds[i-1].rod = bonds[i-1].rod;
  
  tangentnorm = sqrt(bonds[i-1].rod.dot(bonds[i-1].rod));
  atoms[i-1].tangent = bonds[i-1].rod/tangentnorm;


  // compute friction with new values of R on bead i-1
  single_inv_friction(i-1);
  
  for (i = 2; i < Nbeads-1; i++) {

    Rtmp[i] = atoms[i].R;

    atoms[i].t_force = tension(i)*bonds[i].rod-tension(i-1)*bonds[i-1].rod;


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

  atoms[i].t_force = -tension(i-1)*bonds[i-1].rod;



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
void NoTether::final_integrate(double Delta_t)
{
  
  double tmp = Delta_t;




  atoms[0].t_force= tension(0)*bonds[0].rod;

  

  atoms[0].R = Rtmp[0] + tmp*atoms[0].friction*(atoms[0].Fpot+atoms[0].noise+atoms[0].t_force);


    
  for (int i = 1; i < Nbeads-1; i++) {

    atoms[i].t_force = tension(i)*bonds[i].rod-tension(i-1)*bonds[i-1].rod;

    atoms[i].R = Rtmp[i] + tmp*atoms[i].friction*(atoms[i].Fpot+atoms[i].noise+atoms[i].t_force);
    tmpbonds[i-1].rod = (atoms[i].R-atoms[i-1].R)/bondlength;
    
  }

  int i = Nbeads -1;

  atoms[i].t_force = -tension(i-1)*bonds[i-1].rod;


  atoms[i].R = Rtmp[i] + tmp*atoms[i].friction*(atoms[i].Fpot+atoms[i].noise+atoms[i].t_force);
  tmpbonds[i-1].rod = (atoms[i].R-atoms[i-1].R)/bondlength;
  return;
}




void NoTether::calculate_constraint_errors()
{
  for (int mu = 1; mu < Nbeads; mu++) {
    constraint_errors(mu-1) = (atoms[mu].R-atoms[mu-1].R).norm()-bondlength;
  }

  return;

}



namespace gramschmidt {
void gramschmidt(Eigen::Vector3d& x1,Eigen::Vector3d& ey,
		 Eigen::Vector3d& ez)
{
  
  
  
  if (std::abs(std::abs(ey.dot(ez))-std::abs(ey.dot(ey))) < SMALL) {
    ey(0) += 0.1;
    ey(1) += -0.2;
    ey(3) += 0.3;
    ey = ey/sqrt(ey.dot(ey));
  }
  
  if (std::abs(std::abs(x1.dot(ey))-std::abs(x1.dot(x1))) < SMALL) {
    ey(0) += 0.1;
    ey(1) += -0.2;
    ey(3) += 0.3;
    ey = ey/sqrt(ey.dot(ey));
  } else if (std::abs(std::abs(x1.dot(ez))-std::abs(x1.dot(x1))) < SMALL) {
    ez(0) += -0.3;
    ez(1) += 0.12;
    ez(3) += -0.5;
    ez = ez/sqrt(ez.dot(ez));
  }
  
  // Proceed to orthonormalise
  
  ey = ey - ey.dot(x1)*x1;
  ey = ey/sqrt(ey.dot(ey));
  
  ez = ez - ez.dot(ey)*ey - ez.dot(x1)*x1;
  ez = ez/sqrt(ez.dot(ez));
  
  return;
  
}

  
}

  

/* -------------------------------------------------------------------------- */
/* Call to set random configuration. */
/* -------------------------------------------------------------------------- */
void NoTether::init_atoms()
{

  Eigen::Vector3d u(1,0,0);
  Eigen::Vector3d v(0,1,0);
  Eigen::Vector3d w(0,0,1);

  Eigen::Vector3d diffbasis(0,0,0);
  
  double tmp,swich;
  double dum;

  for (int mu = 0; mu < Nbeads-1; mu++) {

    if (mu > 0) {
      gramschmidt::gramschmidt(u,v,w); // redundant in first iteration
    }

    diffbasis(0) = 2*dist(gen);
    
    tmp = 1-diffbasis(0)*diffbasis(0);
    if (tmp > SMALL) {
      dum = dist(gen);
      diffbasis(1) = 2*dum*sqrt(tmp);

      swich = dist(gen);
      
      dum = tmp - diffbasis(1)*diffbasis(1);


      
      if (swich > 0) {
	
	diffbasis(2) = sqrt(tmp - diffbasis(1)*diffbasis(1));
      } else {
	diffbasis(2) = -1*sqrt(tmp - diffbasis(1)*diffbasis(1));
      }
      
      
    } else {
      diffbasis(1) = diffbasis(2) = 0;
    }
    
    bonds[mu].rod(0) = diffbasis(0)*u(0) + diffbasis(1)*v(0) + diffbasis(2)*w(0);
    bonds[mu].rod(1) = diffbasis(0)*u(1) + diffbasis(1)*v(1) + diffbasis(2)*w(1);
    bonds[mu].rod(2) = diffbasis(0)*u(2) + diffbasis(1)*v(2) + diffbasis(2)*w(2);
    
    u = bonds[mu].rod;

  
  }

  
  atoms[0].R(0) = x0(0);
  atoms[0].R(1) = x0(1);

  atoms[0].R(2) = x0(2);


  for (int mu = 0; mu< Nbeads-1; mu++) {
    atoms[mu+1].R = atoms[mu].R + bondlength*bonds[mu].rod;
  }

  // shift so that centre of mass is halfway between the two end points specified

  Eigen::Vector3d com = {0,0,0};

  for (int i = 0; i < Nbeads; i++) {
    com += atoms[i].R;
  }

  com /= Nbeads;

  for (int i = 0; i < Nbeads; i++) {
    atoms[i].R = atoms[i].R + (x0 + xN)/2 - com;
  }  
    
  return;

}


}
