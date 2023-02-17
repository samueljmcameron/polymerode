#include "no_tether.hpp"

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <Eigen/Dense>

namespace BeadRodPmer {

/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
NoTether::NoTether(const std::vector<std::string> & v_line)
  : Polymer(v_line),Gmunu_solver(new Eigen::SimplicialLDLT< SpMat, Eigen::Lower >),
    Hhat_solver (new Eigen::SimplicialLDLT< SpMat, Eigen::Lower >),
    jacob_solver( new Eigen::SparseLU< SpMat >)
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

  tmp_bonds.resize(Eigen::NoChange,Nbeads-1);
  tmp_xs.resize(Eigen::NoChange,Nbeads);

  k_effs.resize(Nbeads-2);
  end_inverses.resize(6);
  end_inverses.setZero();  
	      
}



void NoTether::setup(const Eigen::Ref<const Eigen::Matrix3Xd> &xs) {
  compute_tangents_and_friction(xs);
  set_G();
  set_Hhat();
  set_dCdlambda();
}
  /* This function is outdated, only works when no force is applied. Leaving here for now
     as reminder of what a full step would look like. 

int NoTether::single_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			  Eigen::Ref<Eigen::Matrix3Xd> Fs,double t, double dt,
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

  Fs.setZero();
  for (int index = 0; index < nuc_beads.size(); index++)
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

  */

/* ============================================================================ */
/* Initial integration step. The force vector Fs needs to have been set (either
   by particle-particle interactions or set to zero). */
/* ============================================================================ */
void NoTether::first_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			  Eigen::Ref<Eigen::Matrix3Xd> Fs,double dt)
{
  set_unprojected_noise(dt);
  update_G(0);
  update_Hhat(0);
  compute_noise();
  compute_effective_kappa();
  compute_uc_forces(Fs);
  
  compute_tension(Fs);
  initial_integrate(xs,Fs,dt,0,NONE);
}


/* ============================================================================ */
/* Second integration step. The force vector Fs needs to have been set (either
   by particle-particle interactions or set to zero). */
/* ============================================================================ */
int NoTether::second_step(Eigen::Ref<Eigen::Matrix3Xd> xs,
			  Eigen::Ref<Eigen::Matrix3Xd> Fs,
			  double dt,int itermax)
{
  update_Hhat(0);
  compute_effective_kappa();
  compute_uc_forces(Fs);
  
  compute_tension(Fs);

  int iterations = correct_tension(xs,Fs,dt,itermax,1e-8);


  if (iterations > itermax) {
    for (int i = 0; i < get_Nbeads(); i++) 
      xs.col(i) = tmp_xs.col(i);
  } else {
    final_integrate(xs,Fs,dt,0,NONE);
  }
  compute_tangents_and_friction(xs);
  
  return iterations;
}


  

/* -------------------------------------------------------------------------- */
/* Compute forces on the particles. */
/* -------------------------------------------------------------------------- */
void NoTether::compute_uc_forces(Eigen::Ref<Eigen::Matrix3Xd> Fs)
{

  Fs.col(0) += (end_inverses({0,1,2}).dot(bonds.col(0))*bonds.col(0)
		-1*k_effs(0)*(bonds.col(1)-costhetas(0)*bonds.col(0)));

  Fs.col(0) -= end_inverses({0,1,2});

  Fs.col(1) += (-end_inverses({0,1,2}).dot(bonds.col(0))*bonds.col(0)
		+k_effs(0)*(bonds.col(1)-costhetas(0)*bonds.col(0)
			    - (bonds.col(0)-costhetas(0)*bonds.col(1)))
		- k_effs(1)*(bonds.col(2)-costhetas(1)*bonds.col(1)));


  Fs.col(1) += end_inverses({0,1,2});

  			
  for (int k = 2; k < Nbeads-2; k++) {

    Fs.col(k) += (k_effs(k-1)*(bonds.col(k)-costhetas(k-1)*bonds.col(k-1)
			       -(bonds.col(k-1)-costhetas(k-1)*bonds.col(k)))
		  -k_effs(k)*(bonds.col(k+1)-costhetas(k)*bonds.col(k))
		  +k_effs(k-2)*(bonds.col(k-2)-costhetas(k-2)*bonds.col(k-1)));
    
  }

  int k = Nbeads-2;
  
  Fs.col(k) += (-end_inverses({3,4,5}).dot(bonds.col(k))*bonds.col(k)
		+k_effs(k-1)*(bonds.col(k)-costhetas(k-1)*bonds.col(k-1)
			      -(bonds.col(k-1)-costhetas(k-1)*bonds.col(k)))
		+k_effs(k-2)*(bonds.col(k-2)-costhetas(k-2)*bonds.col(k-1)));
  
  Fs.col(k) += end_inverses({3,4,5});

  k = Nbeads-1;

  Fs.col(k) += (end_inverses({3,4,5}).dot(bonds.col(k-1))*bonds.col(k-1)
		+k_effs(k-2)*(bonds.col(k-2)-costhetas(k-2)*bonds.col(k-1)));


  Fs.col(k) -= end_inverses({3,4,5});

  return;
  
}


double NoTether::Hhat_diag_val(int mu)
{


  double tmp1 = tangents.col(mu).dot(bonds.col(mu));
  tmp1 = tmp1*tmp1;


  double tmp2 =  tangents.col(mu+1).dot(bonds.col(mu));
  tmp2 = tmp2*tmp2;


  return 2./zperp + (1./zpara-1./zperp)*(tmp1+tmp2);

}



double NoTether::dCdlambda_diag_val(int mu)
{


  double tmp1 = tangents.col(mu).dot(tmp_bonds.col(mu));
  tmp1 = tmp1*tangents.col(mu).dot(bonds.col(mu));


  double tmp2 =  tangents.col(mu+1).dot(tmp_bonds.col(mu));
  tmp2 = tmp2*tangents.col(mu+1).dot(bonds.col(mu));

  double tmp3 = tmp_bonds.col(mu).dot(bonds.col(mu));


  return (2*tmp3/zperp + (1./zpara-1./zperp)*(tmp1+tmp2))/tmp_bonds.col(mu).norm();

}


/* NOte that there is a typo in the paper I'm basing this off of
   (Montesi et al, 2005) , which is why there is a difference between
   my off-diagonal H values and theirs. In their equation 34, anywhere
   there is a mu+1, it should be a mu-1. */
double NoTether::Hhat_loweroff_val(int mu)
{

  double tmp1 = bonds.col(mu-1).dot(bonds.col(mu));
  double tmp2 = bonds.col(mu-1).dot(tangents.col(mu));
  double tmp3 = tangents.col(mu).dot(bonds.col(mu));

  return -1./zperp*tmp1-(1./zpara-1./zperp)*tmp2*tmp3;

}

double NoTether::dCdlambda_loweroff_val(int mu)
{

  double tmp1 = bonds.col(mu-1).dot(tmp_bonds.col(mu));
  double tmp2 = bonds.col(mu-1).dot(tangents.col(mu));
  double tmp3 = tangents.col(mu).dot(tmp_bonds.col(mu));

  return (-1./zperp*tmp1-(1./zpara-1./zperp)*tmp2*tmp3)/tmp_bonds.col(mu).norm();

}


double NoTether::dCdlambda_upperoff_val(int mu)
{

  double tmp1 = bonds.col(mu).dot(tmp_bonds.col(mu-1));
  double tmp2 = tmp_bonds.col(mu-1).dot(tangents.col(mu));
  double tmp3 = tangents.col(mu).dot(bonds.col(mu));

  return (-1./zperp*tmp1-(1./zpara-1./zperp)*tmp2*tmp3)/tmp_bonds.col(mu-1).norm();

}



double NoTether::Hhat_endblocks(int first, int second,int mu)
{
  double extra = 0;

  if (first == second) extra = 1;
  
  return (extra/zperp + (1./zpara-1./zperp)*tangents(first,mu)
	  *tangents(second,mu));

}


double NoTether::Hhat_leftside(int first)
{

  double tmp =  tangents.col(0).dot(bonds.col(0))*tangents(first,0);

  
  return -1*(1./zpara-1./zperp)*tmp - bonds(first,0)/zperp;

}


double NoTether::dCdlambda_leftside(int first)
{

  double tmp =  tangents.col(0).dot(tmp_bonds.col(0))*tangents(first,0);

  
  return (-1*(1./zpara-1./zperp)*tmp - tmp_bonds(first,0)/zperp)/tmp_bonds.col(0).norm();

}



double NoTether::Hhat_bottomside(int first)
{
  
  double tmp =  tangents.col(Nbeads-1).dot(bonds.col(Nbeads-2))*tangents(first,Nbeads-1);
  
  
  return (1./zpara-1./zperp)*tmp + bonds(first,Nbeads-2)/zperp;
  
}


double NoTether::dCdlambda_bottomside(int first)
{
  
  double tmp =  tangents.col(Nbeads-1).dot(tmp_bonds.col(Nbeads-2))
    *tangents(first,Nbeads-1);
  
  
  return ((1./zpara-1./zperp)*tmp + tmp_bonds(first,Nbeads-2)/zperp)
    /tmp_bonds.col(Nbeads-2).norm();
  
}
  
  
/* ---------------------------------------------------------------------------- */
/* RHS of G*eta = P. To be used in derived classes. */
/* ---------------------------------------------------------------------------- */
void NoTether::set_rhs_of_G(int offset)
{
  
  
  
  for (int i = 0; i < Nbeads-1; i++) {
    rhs_of_G(i+offset) = bonds.col(i).dot(unprojected_noises.col(i+1)
					  -unprojected_noises.col(i));
    
  }
  
  
  return;
  
}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix G. */
/* -------------------------------------------------------------------------- */
void NoTether::init_G_coeffsmatrix(int offset,std::vector<T> &coeffs)
{
  // just set lower diagonal

  coeffs.push_back(T(offset,offset,2));
  
  for (int i = 1; i < Nbeads-1; i++) {
    coeffs.push_back(T(i+offset,i+offset,2));
    coeffs.push_back(T(i+offset,i+offset-1,-costhetas(i-1)));

  }
  
  return;

  
}


/* -------------------------------------------------------------------------- */
/* Update matrix G. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void NoTether::update_G(int offset)
{
  // update first column by hand



  Gmunu.coeffRef(offset+1,offset) = -costhetas(0);


  // update middle columns
  for (int k = 1; k < Nbeads-2; k++) {


    int count = 0;
    for (SpMat::InnerIterator it(Gmunu,k+offset); it; ++it) {
      if (count == 1) {
	it.valueRef() = -costhetas(k);

      }
      count += 1;
    }
  }


  return;
}


void NoTether::set_bdets_and_tdets(int offset)
{

  int mu;
  for (int i = 0; i < Nbeads - 2; i++ ) {

    mu = Nbeads-2-i;
    bDets(mu-1+offset) = 2*bDets(mu+offset)
      - costhetas(mu-1)*costhetas(mu-1)*bDets(mu+1+offset);

    mu = i+2;

    tDets(mu) = 2*tDets(mu-1)-costhetas(mu-2)*costhetas(mu-2)*tDets(mu-2);


  }
  
}


/* ---------------------------------------------------------------------------- */
/* RHS of Hhat*lambda = Q. */
/* ---------------------------------------------------------------------------- */
void NoTether::set_rhs_of_Hhat(int offset,const Eigen::Ref<const Eigen::Matrix3Xd> &Fs)
{

  
  
  for (int i = 0; i< Nbeads-1; i++) {


    rhs_of_Hhat(i+offset) = bonds.col(i).dot(frictions[i+1]*(Fs.col(i+1)
								  +noises.col(i+1))
					     -frictions[i]*(Fs.col(i)
								 +noises.col(i)));
    
  }

  return;
  
}

/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix H hat. */
/* -------------------------------------------------------------------------- */
void NoTether::init_Hhat_coeffsmatrix(int offset,std::vector<T> & coeffs)
{



  
  coeffs.push_back(T(offset,offset,Hhat_diag_val(0)));
  
  
  for (int i = 1; i < Nbeads-1; i++) {
    coeffs.push_back(T(i+offset,i+offset,Hhat_diag_val(i)));
    coeffs.push_back(T(i+offset,i+offset-1,Hhat_loweroff_val(i)));

  }


  return;

  
}


/* -------------------------------------------------------------------------- */
/* Update matrix M. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void NoTether::update_Hhat(int offset)
{

    
  Hhat.coeffRef(offset,offset) = Hhat_diag_val(0);

  Hhat.coeffRef(offset+1,offset) = Hhat_loweroff_val(1);

  // update middle columns
  for (int k = 1; k < Nbeads-2; k++) {

    int count = 0;
    for (SpMat::InnerIterator it(Hhat,k+offset); it; ++it) {
      if (count == 0) {
	it.valueRef() = Hhat_diag_val(k);

      } else {
	it.valueRef() = Hhat_loweroff_val(k+1);
      }

      count += 1;
    }
  }

  Hhat.coeffRef(offset+Nbeads-2,offset+Nbeads-2) = Hhat_diag_val(Nbeads-2);
  return;
}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix dCdlambda. */
/* -------------------------------------------------------------------------- */
void NoTether::init_dCdlambda_coeffsmatrix(int offset,std::vector<T> &coeffs)
{

  // initializing the vector as the full (both upper and lower) part of Hhat, since
  // this matrix won't be symmetric.

  coeffs.push_back(T(offset,offset,Hhat_diag_val(0)));
  
  for (int i = 1; i < Nbeads-1; i++) {
    coeffs.push_back(T(i-1+offset,i+offset,Hhat_loweroff_val(i)));
    coeffs.push_back(T(i+offset,i+offset,Hhat_diag_val(i)));
    coeffs.push_back(T(i+offset,i-1+offset,Hhat_loweroff_val(i)));


  }

  return ;

}

void NoTether::update_dCdlambda(double Delta_t,int offset)
{



  dCdlambda.coeffRef(offset,offset) = -dCdlambda_diag_val(0)*Delta_t;
  dCdlambda.coeffRef(offset+1,offset)  = -dCdlambda_loweroff_val(1)*Delta_t;
  

  // update middle columns
  for (int k = 1; k < Nbeads-2; k++) {


    int count = 0;
    for (SpMat::InnerIterator it(dCdlambda,k+offset); it; ++it) {

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
  dCdlambda.coeffRef(offset+Nbeads-3, offset+Nbeads-2) = - dCdlambda_upperoff_val(Nbeads-2)*Delta_t;			 
  dCdlambda.coeffRef(offset+Nbeads-2, offset+Nbeads-2) = - dCdlambda_diag_val(Nbeads-2)*Delta_t;


  return;

}

void NoTether::update_noise(int offset)
{
  int i = 0;
  noises.col(i) = (unprojected_noises.col(i)
			 + dummy_for_noise(i+offset)*bonds.col(i));
  
  
  
  for (i = 1; i < Nbeads-1; i++) {
    noises.col(i) = (unprojected_noises.col(i)
			   + dummy_for_noise(i+offset)*bonds.col(i)
			   - dummy_for_noise(i+offset-1)*bonds.col(i-1));
  }


  i = Nbeads-1;

  noises.col(i) = (unprojected_noises.col(i)
			 - dummy_for_noise(i+offset-1)*bonds.col(i-1));
  return;

}





void NoTether::initial_integrate(Eigen::Ref<Eigen::Matrix3Xd> xs,
				 const Eigen::Ref<const Eigen::Matrix3Xd> &Fs,
				 double Delta_t,int offset, PTYPE cflag)
{
  
  double tmp = Delta_t/2.0;

  double tangentnorm;

  
  int i = 0;
  tmp_xs.col(0) = xs.col(0);

  t_forces.col(0) = tension(offset)*bonds.col(0);

  if (cflag == SINGLE || cflag == DOUBLE) t_forces.col(0) -= tension({0,1,2});

  xs.col(0) += tmp*frictions[0]*(Fs.col(0)+noises.col(0)+t_forces.col(0));

  i = 1;
  tmp_xs.col(i) = xs.col(i);


  t_forces.col(i) = tension(offset+i)*bonds.col(i)-tension(offset+i-1)*bonds.col(i-1);
  
  
  xs.col(i) += tmp*frictions[i]*(Fs.col(i)+noises.col(i)+t_forces.col(i));

  bonds.col(i-1) = (xs.col(i)-xs.col(i-1))/bondlength;
  tmp_bonds.col(i-1) = bonds.col(i-1);
  
  tangentnorm = sqrt(bonds.col(i-1).dot(bonds.col(i-1)));
  tangents.col(i-1) = bonds.col(i-1)/tangentnorm;


  // compute friction with new values of R on bead i-1
  single_inv_friction(i-1);
  
  for (i = 2; i < Nbeads-1; i++) {

    tmp_xs.col(i) = xs.col(i);

    t_forces.col(i) = tension(offset+i)*bonds.col(i)-tension(offset-1+i)*bonds.col(i-1);


    xs.col(i) += tmp*frictions[i]*(Fs.col(i)+noises.col(i)+t_forces.col(i));

    bonds.col(i-1) = (xs.col(i)-xs.col(i-1))/bondlength;
    tmp_bonds.col(i-1) = bonds.col(i-1);

    costhetas(i-2) = bonds.col(i-1).dot(bonds.col(i-2))
      /(bonds.col(i-1).norm()*bonds.col(i-2).norm());


    tangentnorm = sqrt((bonds.col(i-1)+bonds.col(i-2)
			).dot(bonds.col(i-1)+bonds.col(i-2)));


    tangents.col(i-1) = (bonds.col(i-1)+bonds.col(i-2))/tangentnorm;

    single_inv_friction(i-1);

  }


  i = Nbeads -1;

  tmp_xs.col(i) = xs.col(i);

  t_forces.col(i) = -tension(offset-1+i)*bonds.col(i-1);

  if (cflag == DOUBLE) {
    int ind = offset+Nbeads-1;
    t_forces.col(i) -= tension({ind,ind+1,ind+2});
  }


  xs.col(i) += tmp*frictions[i]*(Fs.col(i)+noises.col(i)+t_forces.col(i));

  bonds.col(i-1) = (xs.col(i)-xs.col(i-1))/bondlength;
  tmp_bonds.col(i-1) = bonds.col(i-1);


  costhetas(i-2) = bonds.col(i-1).dot(bonds.col(i-2))
    /(bonds.col(i-1).norm()*bonds.col(i-2).norm());

  
  tangentnorm = sqrt((bonds.col(i-1)+bonds.col(i-2)
		      ).dot(bonds.col(i-1)+bonds.col(i-2)));
  
  
  tangents.col(i-1) = (bonds.col(i-1)+bonds.col(i-2))/tangentnorm;
  
  single_inv_friction(i-1);
  
  tangentnorm = sqrt(bonds.col(i-1).dot(bonds.col(i-1)));
  tangents.col(i) = bonds.col(i-1)/tangentnorm;

  single_inv_friction(i);

  return;
}

/*----------------------------------------------------------------------------*/
/* Computes the final bead positions and unnormalised bond tangents. */
/*----------------------------------------------------------------------------*/
void NoTether::final_integrate(Eigen::Ref<Eigen::Matrix3Xd> xs,
			       const Eigen::Ref<const Eigen::Matrix3Xd> &Fs,
			       double Delta_t, int offset,PTYPE cflag)
{
  
  double tmp = Delta_t;


  t_forces.col(0)= tension(offset)*bonds.col(0);

  if (cflag == SINGLE || cflag == DOUBLE) {
    t_forces.col(0) -= tension({0,1,2});
  }

  xs.col(0) = tmp_xs.col(0) + tmp*frictions[0]*(Fs.col(0)+noises.col(0)+t_forces.col(0));


    
  for (int i = 1; i < Nbeads-1; i++) {

    t_forces.col(i) = tension(offset+i)*bonds.col(i)-tension(offset-1+i)*bonds.col(i-1);

    xs.col(i) = tmp_xs.col(i) + tmp*frictions[i]*(Fs.col(i)+noises.col(i)+t_forces.col(i));
    tmp_bonds.col(i-1) = (xs.col(i)-xs.col(i-1))/bondlength;
    
  }

  int i = Nbeads -1;

  t_forces.col(i) = -tension(offset-1+i)*bonds.col(i-1);

  if (cflag == DOUBLE) {
    int ind = offset+Nbeads-1;
    t_forces.col(i) -= tension({ind,ind+1,ind+2});
  }

  xs.col(i) = tmp_xs.col(i) + tmp*frictions[i]*(Fs.col(i)+noises.col(i)+t_forces.col(i));
  tmp_bonds.col(i-1) = (xs.col(i)-xs.col(i-1))/bondlength;

  
  return;
}


  
  
  
void NoTether::calculate_constraint_errors(int offset,
					   const Eigen::Ref<const Eigen::Matrix3Xd> &xs)
{
  for (int mu = 1; mu < Nbeads; mu++) {
    constraint_errors(mu-1+offset) = (xs.col(mu)-xs.col(mu-1)).norm()-bondlength;
  }
  
  return;
  
}

/* -------------------------------------------------------------------------- */
/* Initialse G in G*eta = P (only call once). */
/* -------------------------------------------------------------------------- */
void NoTether::set_G()
{

  std::vector<T> coefficients;
  init_G_coeffsmatrix(0,coefficients);

  Gmunu.setFromTriplets(coefficients.begin(),coefficients.end());

  Gmunu_solver->analyzePattern(Gmunu);

  return;

}




void NoTether::compute_effective_kappa()
{


  bDets(Nbeads-1) = 1.0;
  tDets(0) = 1.0;

  bDets(Nbeads-2) = 2;
  tDets(1) = 2;

  set_bdets_and_tdets(0);

  /*
  Eigen::MatrixXd Gtmp(Nbeads-1,Nbeads-1);

  Gtmp.setZero();

  int i = 0;

  Gtmp(i,i) = 2.0;
  Gtmp(i,i+1) = -costhetas(i);

  for (i = 1; i < Nbeads-2; i++) {
    Gtmp(i,i) = 2.0;
    Gtmp(i,i-1) = -costhetas(i-1);
    Gtmp(i,i+1) = -costhetas(i);
  }

  Gtmp(i,i) = 2.0;
  Gtmp(i,i-1) = -costhetas(i-1);

  std::cout << " Ginv = " << std::endl;

  Gtmp = Gtmp.inverse();
  */
  double gDet = tDets(Nbeads-1);


  /*
  for (int i = 1; i < Nbeads-1; i++) {
    std::cout << "correct Ginv(i,i-1) = " << Gtmp(i,i-1) << ", "
	      << "incorrect Ginve(i,i-1) = " << costhetas(i-1)*tDets(i-1)*bDets(i+1)/gDet
	      << std::endl;
  }
  */
  


  for (int i = 0; i < Nbeads-2; i++) {

    k_effs(i) = (kappa + temp*bondlength*costhetas(i)*tDets(i)*bDets(i+2)/gDet
		 )/(bondlength*bondlength);
    //    k_effs(i) = kappa/(bondlength*bondlength);

    
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


  Hhat_solver->analyzePattern(Hhat);

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

  jacob_solver->analyzePattern(dCdlambda);

  return;

}


void NoTether::compute_noise()
{


  set_rhs_of_G(0);
  Gmunu_solver->factorize(Gmunu);

  dummy_for_noise =  Gmunu_solver->solve(rhs_of_G);

  update_noise(0);
  
}


void NoTether::compute_tension(const Eigen::Ref<const Eigen::Matrix3Xd> &Fs)
{



  set_rhs_of_Hhat(0,Fs);

  Hhat_solver->factorize(Hhat);
  tension =  Hhat_solver->solve(rhs_of_Hhat);
  
}
/*
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
*/
int NoTether::correct_tension(Eigen::Ref<Eigen::Matrix3Xd> xs,
			      const Eigen::Ref<const Eigen::Matrix3Xd> &Fs,
			      double Delta_t,int itermax,double tolerance)
{

  // set C_mu and dC_mu/dlambda_nu
  final_integrate(xs,Fs,Delta_t,0,NONE);
  calculate_constraint_errors(0,xs);
  update_dCdlambda(Delta_t,0);
  
  
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


  while (negative_tension_change.norm()/tension.norm() > tolerance && count <= itermax) {



    final_integrate(xs,Fs,Delta_t,0,NONE);
    calculate_constraint_errors(0,xs);
    update_dCdlambda(Delta_t,0);

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
 
}
