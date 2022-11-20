#include "no_tether.hpp"
#include "input.hpp"
#include "initialise.hpp"

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

  

  // set ^ configuration, with the bottoms of the ^ being x0 and xN
  init_atoms_caret();

}





int NoTether::single_step(double t, double dt,
			  const std::vector<std::vector<double>> & dFdX_i,
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
      atoms[i].R = Rtmp[i];

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
void NoTether::init_atoms_rand()
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

/* -------------------------------------------------------------------------- */
/* Call to set line configuration so that polymer intersects initial points
   x0 and xN with same centre of mass as the initial points. */
/* -------------------------------------------------------------------------- */
void NoTether::init_atoms_line()
{

  double length_init = (xN-x0).norm();
  
  double t1 = 0.5*(1-(Nbeads-1)*bondlength/length_init);
  for (int i = 0; i < Nbeads; i++) 
    atoms[i].R = x0 + (xN-x0)*(t1+i*bondlength/length_init);
    
  return;

}

void NoTether::init_atoms_caret()
{


  
  Eigen::Vector3d nhat;

  nhat(0) = 2*dist(gen);
  nhat(1) = 2*dist(gen);
  nhat(2) = 2*dist(gen);
  
  Eigen::Vector3d dd = xN-x0;

  if (dd.norm() > bondlength*(Nbeads-1)) {
    throw std::runtime_error("|x0-xN| is longer than the polymer.");
  }


  if (dd.norm() < SMALL) {
    throw std::runtime_error("|x0-xN| is longer than the polymer.");
  }

  if (std::abs(dd(2)) < SMALL) {
    if (std::abs(dd(1)) < SMALL) {
      nhat(0) = 0;
      // loop here is to ensure that not all components of normal vector are zero
      while (std::abs(nhat(1)) < SMALL) {
	nhat(1) = 2*dist(gen);
      }
    } else {
      nhat(1) = -nhat(0)*dd(0)/dd(1);
    }
  } else {
    nhat(2) = -(nhat(1)*dd(1) + nhat(0)*dd(0))/dd(2);
  }

  nhat = nhat/nhat.norm();


  double fakelength;


  // since the parameterisation only works for an odd number of beads, we overshoot
  // the point x0 for an even number of beads by pretending the system only has
  // Nbeads-1 beads, so the final bead gets put past the point x0
  if (Nbeads % 2 == 0) {
    fakelength = (Nbeads-2)*bondlength;
  } else {
    fakelength = (Nbeads-1)*bondlength;
  }
  

  double gamma = 0.5*sqrt(fakelength*fakelength-dd.squaredNorm());

  Eigen::Vector3d alpha = 0.5*dd+gamma*nhat;
  Eigen::Vector3d beta = -0.5*dd+gamma*nhat;


  double tp = 0;
  double dtp = 1.0*bondlength/fakelength;

  
  for (int i = 0; i < Nbeads; i++) {
    
    tp = i*dtp;

    // integer division rounds toward zero
    if (i < Nbeads/2) {
      atoms[i].R = x0 + 2*alpha*tp;
    } else {
      atoms[i].R = xN + 2*beta*(1-tp);
    }
  }
  return;

}

}
