#include "polymer.hpp"
#include "input.hpp"


#include <iostream>
#include <cmath>
#include <stdexcept>

#define SMALL 1e-14


namespace BeadRodPmer {
/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
Polymer::Polymer(const std::vector<std::string> & splitvec)
  : gen(0), dist(-0.5,0.5), flag_x0(false),flag_xN(false),
    flag_initdoubleteth(false),equilibration_steps(0)
{


  int nargs = splitvec.size();
  int iarg = 0;

  int number_of_nuc_beads = 0;
  int number_of_nuc_strengths = 0;
  int number_of_nuc_maxs = 0;
  int number_of_nuc_widths = 0;
  int seed;
  
  while (iarg < nargs) {
    if (splitvec[iarg] == "beads") {
      input::isInt(splitvec[iarg+1],Nbeads,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "bondlength") {
      input::isDouble(splitvec[iarg+1],bondlength,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "zeta_para") {
      input::isDouble(splitvec[iarg+1],zpara,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "zeta_perp") {
      input::isDouble(splitvec[iarg+1],zperp,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "temp") {
      input::isDouble(splitvec[iarg+1],temp,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "bending") {
      input::isDouble(splitvec[iarg+1],kappa,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "x0") {
      input::isDouble(splitvec[iarg+1],x0(0),splitvec[iarg]);
      input::isDouble(splitvec[iarg+2],x0(1),splitvec[iarg]);
      input::isDouble(splitvec[iarg+3],x0(2),splitvec[iarg]);
      flag_x0 = true;
      iarg += 4;
    } else if (splitvec[iarg] == "xN") {
      input::isDouble(splitvec[iarg+1],xN(0),splitvec[iarg]);
      input::isDouble(splitvec[iarg+2],xN(1),splitvec[iarg]);
      input::isDouble(splitvec[iarg+3],xN(2),splitvec[iarg]);
      flag_xN = true;
      iarg += 4;
    } else if (splitvec[iarg] == "seed") {
      input::isInt(splitvec[iarg+1],seed,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "equilibration_steps") {
      input::isInt(splitvec[iarg+1],equilibration_steps,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "init_ends") {
      input::isDouble(splitvec[iarg+1],initspringK,splitvec[iarg]);
      input::isDouble(splitvec[iarg+2],initdt,splitvec[iarg]);
      input::isDouble(splitvec[iarg+3],inittolerance,splitvec[iarg]);
      input::isInt(splitvec[iarg+4],equilibration_steps,splitvec[iarg]);
      flag_initdoubleteth = true;
      iarg += 5;
    } else if (splitvec[iarg] == "nuc_beads") {
      input::isInt(splitvec[iarg+1],number_of_nuc_beads,splitvec[iarg]);
      nuc_beads.resize(number_of_nuc_beads);
      for (int i = 1; i <= number_of_nuc_beads; i++) 
	input::isInt(splitvec[iarg+1+i],nuc_beads[i-1],splitvec[iarg]);
      iarg += number_of_nuc_beads+2;
    } else if (splitvec[iarg] == "nuc_strengths") {
      input::isInt(splitvec[iarg+1],number_of_nuc_strengths,splitvec[iarg]);
      nuc_strengths.resize(number_of_nuc_strengths);
      for (int i = 1; i <= number_of_nuc_strengths; i++) 
	input::isDouble(splitvec[iarg+1+i],nuc_strengths[i-1],splitvec[iarg]);
      iarg += number_of_nuc_strengths+2;
    } else if (splitvec[iarg] == "nuc_maxs") {
      input::isInt(splitvec[iarg+1],number_of_nuc_maxs,splitvec[iarg]);
      nuc_maxs.resize(number_of_nuc_maxs);
      for (int i = 1; i <= number_of_nuc_maxs; i++) 
	input::isDouble(splitvec[iarg+1+i],nuc_maxs[i-1],splitvec[iarg]);
      iarg += number_of_nuc_maxs+2;
    } else if (splitvec[iarg] == "nuc_widths") {
      input::isInt(splitvec[iarg+1],number_of_nuc_widths,splitvec[iarg]);
      nuc_widths.resize(number_of_nuc_widths);
      for (int i = 1; i <= number_of_nuc_widths; i++) 
	input::isDouble(splitvec[iarg+1+i],nuc_widths[i-1],splitvec[iarg]);
      iarg += number_of_nuc_widths+2;
    } else {
      throw std::runtime_error("Error: invalid argument for build_polymer.");
    }
  }


  if (!flag_x0) 
    throw std::runtime_error("Need to specify x0 for polymer.");
  if (!flag_xN) 
    throw std::runtime_error("Need to specify xN for polymer.");
  if (!flag_initdoubleteth) 
    throw std::runtime_error("Need to specify initialisation tolerance for polymer.");

  if   (number_of_nuc_beads != number_of_nuc_strengths)
    throw std::runtime_error("Error: need same number of nuc_beads as nuc_strengths.");
  if   (number_of_nuc_beads != number_of_nuc_maxs)
    throw std::runtime_error("Error: need same number of nuc_beads as nuc_maxs.");
  if   (number_of_nuc_beads != number_of_nuc_widths)
    throw std::runtime_error("Error: need same number of nuc_beads as nuc_widths.");
  if   (number_of_nuc_strengths != number_of_nuc_maxs)
    throw std::runtime_error("Error: need same number of nuc_strengths as nuc_maxs.");
  if   (number_of_nuc_strengths != number_of_nuc_widths)
    throw std::runtime_error("Error: need same number of nuc_strengths as nuc_widths.");
  if   (number_of_nuc_maxs != number_of_nuc_widths)
    throw std::runtime_error("Error: need same number of nuc_maxs as nuc_widths.");  

  
  
    
  resize(Nbeads);

  tmp_bonds.resize(Eigen::NoChange,Nbeads-1);
  tmp_xs.resize(Eigen::NoChange,Nbeads);


  costhetas.resize(Nbeads-2);
  k_effs.resize(Nbeads-2);
  end_inverses.resize(6);
  end_inverses.setZero();

  gen.seed(seed);


  
}



/* ---------------------------------------------------------------------------- */
/* Compute the unit vectors tangent to the beads and rods, and then friction. */
/* Also computes cos(theta_i) relevant for potential force. */
/* Call only at the start of the time step (NOT the midstep). */
/* ---------------------------------------------------------------------------- */
void Polymer::compute_tangents_and_friction()
{

  double tangentnorm;


  bonds.col(0) = (xs.col(1)-xs.col(0))/bondlength;
  
  tangentnorm = sqrt(bonds.col(0).dot(bonds.col(0)));
  tangents.col(0) = bonds.col(0)/tangentnorm;


  single_inv_friction(0);


  for (int i = 1; i < Nbeads-1; i++) {

    bonds.col(i) = (xs.col(i+1)-xs.col(i))/bondlength;
    costhetas(i-1) = bonds.col(i).dot(bonds.col(i-1));
      
    tangentnorm = sqrt((bonds.col(i)+bonds.col(i-1)
			).dot(bonds.col(i)+bonds.col(i-1)));


    tangents.col(i) = (bonds.col(i)+bonds.col(i-1))/tangentnorm;

    single_inv_friction(i);

  }
  

  tangentnorm = sqrt(bonds.col(Nbeads-2).dot(bonds.col(Nbeads-2)));
  tangents.col(Nbeads-1) = bonds.col(Nbeads-2)/tangentnorm;

  single_inv_friction(Nbeads-1);
  
  return;
  
}



/* ---------------------------------------------------------------------------- */
/* Generate unprojected noises with correct variance. */
/* ---------------------------------------------------------------------------- */
void Polymer::set_unprojected_noise(double Delta_t)
{

  double dum = sqrt(24*temp/Delta_t);
  double para = sqrt(zpara);
  double perp = sqrt(zperp);


  
  Eigen::Vector3d rands;
  
  for (int i = 0; i < Nbeads; i++) {

    for (int j = 0; j < 3; j++ ) {

      rands(j) = dist(gen);
    }
    
    unprojected_noises.col(i) = dum*(perp*rands
				      +(para-perp)*(tangents.col(i).dot(rands))
				      *tangents.col(i));

  }

  return;
}






/* -------------------------------------------------------------------------- */
/* Compute forces on the particles. */
/* -------------------------------------------------------------------------- */
void Polymer::compute_uc_forces()
{

  Fpots.col(0) = (end_inverses({0,1,2}).dot(bonds.col(0))*bonds.col(0)
		   -1*k_effs(0)*(bonds.col(1)-costhetas(0)*bonds.col(0)));

  Fpots.col(0) -= end_inverses({0,1,2});

  Fpots.col(1) = (-end_inverses({0,1,2}).dot(bonds.col(0))*bonds.col(0)
		   +k_effs(0)*(bonds.col(1)-costhetas(0)*bonds.col(0)
			       - (bonds.col(0)-costhetas(0)*bonds.col(1)))
		   - k_effs(1)*(bonds.col(2)-costhetas(1)*bonds.col(1)));


  Fpots.col(1) += end_inverses({0,1,2});

  			
  for (int k = 2; k < Nbeads-2; k++) {

    Fpots.col(k) = (k_effs(k-1)*(bonds.col(k)-costhetas(k-1)*bonds.col(k-1)
				  -(bonds.col(k-1)-costhetas(k-1)*bonds.col(k)))
		     -k_effs(k)*(bonds.col(k+1)-costhetas(k)*bonds.col(k))
		     +k_effs(k-2)*(bonds.col(k-2)-costhetas(k-2)*bonds.col(k-1)));
    
  }

  int k = Nbeads-2;
  
  Fpots.col(k) = (-end_inverses({3,4,5}).dot(bonds.col(k))*bonds.col(k)
		   +k_effs(k-1)*(bonds.col(k)-costhetas(k-1)*bonds.col(k-1)
				 -(bonds.col(k-1)-costhetas(k-1)*bonds.col(k)))
		   +k_effs(k-2)*(bonds.col(k-2)-costhetas(k-2)*bonds.col(k-1)));
  
  Fpots.col(k) += end_inverses({3,4,5});

  k = Nbeads-1;

  Fpots.col(k) = (end_inverses({3,4,5}).dot(bonds.col(k-1))*bonds.col(k-1)
		   +k_effs(k-2)*(bonds.col(k-2)-costhetas(k-2)*bonds.col(k-1)));


  Fpots.col(k) -= end_inverses({3,4,5});

  return;
  
}


void Polymer::single_inv_friction(int i)
{
  double tmp = 1./zpara-1./zperp;
  frictions[i] = tmp*tangents.col(i)*tangents.col(i).transpose();
  
  for (int j = 0; j < 3; j++ )

    frictions[i](j,j) += 1./zperp;

  return;

}

void Polymer::add_external_force(const Eigen::Vector3d &dFdX, int i)
{
  Fpots.col(i) += -dFdX;
  return;

}

double Polymer::Hhat_diag_val(int mu)
{


  double tmp1 = tangents.col(mu).dot(bonds.col(mu));
  tmp1 = tmp1*tmp1;


  double tmp2 =  tangents.col(mu+1).dot(bonds.col(mu));
  tmp2 = tmp2*tmp2;


  return 2./zperp + (1./zpara-1./zperp)*(tmp1+tmp2);

}



double Polymer::dCdlambda_diag_val(int mu)
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
double Polymer::Hhat_loweroff_val(int mu)
{

  double tmp1 = bonds.col(mu-1).dot(bonds.col(mu));
  double tmp2 = bonds.col(mu-1).dot(tangents.col(mu));
  double tmp3 = tangents.col(mu).dot(bonds.col(mu));

  return -1./zperp*tmp1-(1./zpara-1./zperp)*tmp2*tmp3;

}

double Polymer::dCdlambda_loweroff_val(int mu)
{

  double tmp1 = bonds.col(mu-1).dot(tmp_bonds.col(mu));
  double tmp2 = bonds.col(mu-1).dot(tangents.col(mu));
  double tmp3 = tangents.col(mu).dot(tmp_bonds.col(mu));

  return (-1./zperp*tmp1-(1./zpara-1./zperp)*tmp2*tmp3)/tmp_bonds.col(mu).norm();

}


double Polymer::dCdlambda_upperoff_val(int mu)
{

  double tmp1 = bonds.col(mu).dot(tmp_bonds.col(mu-1));
  double tmp2 = tmp_bonds.col(mu-1).dot(tangents.col(mu));
  double tmp3 = tangents.col(mu).dot(bonds.col(mu));

  return (-1./zperp*tmp1-(1./zpara-1./zperp)*tmp2*tmp3)/tmp_bonds.col(mu-1).norm();

}



double Polymer::Hhat_endblocks(int first, int second,int mu)
{
  double extra = 0;

  if (first == second) extra = 1;
  
  return (extra/zperp + (1./zpara-1./zperp)*tangents(first,mu)
	  *tangents(second,mu));

}


double Polymer::Hhat_leftside(int first)
{

  double tmp =  tangents.col(0).dot(bonds.col(0))*tangents(first,0);

  
  return -1*(1./zpara-1./zperp)*tmp - bonds(first,0)/zperp;

}


double Polymer::dCdlambda_leftside(int first)
{

  double tmp =  tangents.col(0).dot(tmp_bonds.col(0))*tangents(first,0);

  
  return (-1*(1./zpara-1./zperp)*tmp - tmp_bonds(first,0)/zperp)/tmp_bonds.col(0).norm();

}



double Polymer::Hhat_bottomside(int first)
{
  
  double tmp =  tangents.col(Nbeads-1).dot(bonds.col(Nbeads-2))*tangents(first,Nbeads-1);
  
  
  return (1./zpara-1./zperp)*tmp + bonds(first,Nbeads-2)/zperp;
  
}


double Polymer::dCdlambda_bottomside(int first)
{
  
  double tmp =  tangents.col(Nbeads-1).dot(tmp_bonds.col(Nbeads-2))
    *tangents(first,Nbeads-1);
  
  
  return ((1./zpara-1./zperp)*tmp + tmp_bonds(first,Nbeads-2)/zperp)
    /tmp_bonds.col(Nbeads-2).norm();
  
}
  
  
/* ---------------------------------------------------------------------------- */
/* RHS of G*eta = P. To be used in derived classes. */
/* ---------------------------------------------------------------------------- */
void Polymer::set_rhs_of_G(int offset)
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
void Polymer::init_G_coeffsmatrix(int offset,std::vector<T> &coeffs)
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
void Polymer::update_G(int offset)
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


void Polymer::set_bdets_and_tdets(int offset)
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
void Polymer::set_rhs_of_Hhat(int offset)
{

  
  
  for (int i = 0; i< Nbeads-1; i++) {


    rhs_of_Hhat(i+offset) = bonds.col(i).dot(frictions[i+1]*(Fpots.col(i+1)
								  +noises.col(i+1))
					     -frictions[i]*(Fpots.col(i)
								 +noises.col(i)));
    
  }

  return;
  
}

/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix H hat. */
/* -------------------------------------------------------------------------- */
void Polymer::init_Hhat_coeffsmatrix(int offset,std::vector<T> & coeffs)
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
void Polymer::update_Hhat(int offset)
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
void Polymer::init_dCdlambda_coeffsmatrix(int offset,std::vector<T> &coeffs)
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

void Polymer::update_dCdlambda(double Delta_t,int offset)
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

void Polymer::update_noise(int offset)
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





void Polymer::initial_integrate(double Delta_t,int offset, PTYPE cflag)
{
  
  double tmp = Delta_t/2.0;

  double tangentnorm;

  
  int i = 0;
  tmp_xs.col(0) = xs.col(0);

  t_forces.col(0) = tension(offset)*bonds.col(0);

  if (cflag == SINGLE || cflag == DOUBLE) t_forces.col(0) -= tension({0,1,2});

  xs.col(0) += tmp*frictions[0]*(Fpots.col(0)+noises.col(0)+t_forces.col(0));

  i = 1;
  tmp_xs.col(i) = xs.col(i);


  t_forces.col(i) = tension(offset+i)*bonds.col(i)-tension(offset+i-1)*bonds.col(i-1);
  
  
  xs.col(i) += tmp*frictions[i]*(Fpots.col(i)+noises.col(i)+t_forces.col(i));

  bonds.col(i-1) = (xs.col(i)-xs.col(i-1))/bondlength;
  tmp_bonds.col(i-1) = bonds.col(i-1);
  
  tangentnorm = sqrt(bonds.col(i-1).dot(bonds.col(i-1)));
  tangents.col(i-1) = bonds.col(i-1)/tangentnorm;


  // compute friction with new values of R on bead i-1
  single_inv_friction(i-1);
  
  for (i = 2; i < Nbeads-1; i++) {

    tmp_xs.col(i) = xs.col(i);

    t_forces.col(i) = tension(offset+i)*bonds.col(i)-tension(offset-1+i)*bonds.col(i-1);


    xs.col(i) += tmp*frictions[i]*(Fpots.col(i)+noises.col(i)+t_forces.col(i));

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


  xs.col(i) += tmp*frictions[i]*(Fpots.col(i)+noises.col(i)+t_forces.col(i));

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
void Polymer::final_integrate(double Delta_t, int offset,PTYPE cflag)
{
  
  double tmp = Delta_t;


  t_forces.col(0)= tension(offset)*bonds.col(0);

  if (cflag == SINGLE || cflag == DOUBLE) {
    t_forces.col(0) -= tension({0,1,2});
  }

  xs.col(0) = tmp_xs.col(0) + tmp*frictions[0]*(Fpots.col(0)+noises.col(0)+t_forces.col(0));


    
  for (int i = 1; i < Nbeads-1; i++) {

    t_forces.col(i) = tension(offset+i)*bonds.col(i)-tension(offset-1+i)*bonds.col(i-1);

    xs.col(i) = tmp_xs.col(i) + tmp*frictions[i]*(Fpots.col(i)+noises.col(i)+t_forces.col(i));
    tmp_bonds.col(i-1) = (xs.col(i)-xs.col(i-1))/bondlength;
    
  }

  int i = Nbeads -1;

  t_forces.col(i) = -tension(offset-1+i)*bonds.col(i-1);

  if (cflag == DOUBLE) {
    int ind = offset+Nbeads-1;
    t_forces.col(i) -= tension({ind,ind+1,ind+2});
  }

  xs.col(i) = tmp_xs.col(i) + tmp*frictions[i]*(Fpots.col(i)+noises.col(i)+t_forces.col(i));
  tmp_bonds.col(i-1) = (xs.col(i)-xs.col(i-1))/bondlength;

  
  return;
}


  
  
  
void Polymer::calculate_constraint_errors(int offset)
{
  for (int mu = 1; mu < Nbeads; mu++) {
    constraint_errors(mu-1+offset) = (xs.col(mu)-xs.col(mu-1)).norm()-bondlength;
  }
  
  return;
  
}

void Polymer::rescale_positions(bool fromfront)
{
  
  if (fromfront) {
    for (int i = 1; i < Nbeads; i++) {
      xs.col(i) = xs.col(i-1) + bondlength*bonds.col(i-1);
    }
  } else {
    for (int i = Nbeads-2; i >= 0; i--) {
      xs.col(i) = xs.col(i+1) - bondlength*bonds.col(i-1);
    }
  }
  return;
  
}

int Polymer::get_Nbeads() const
{
  return Nbeads;
}


double Polymer::get_temp() const
{
  return temp;
}

double Polymer::get_zpara() const
{
  return zpara;
}

double Polymer::get_zperp() const
{
  return zperp;
}

double Polymer::get_bondlength() const
{
  return bondlength;
}

double Polymer::get_timescale(double dt) const
{
  return dt*temp/(bondlength*zperp);
}



Eigen::Vector3d Polymer::get_x0() const
{
  return x0;
}

Eigen::Vector3d Polymer::get_xN() const
{
  return xN;
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
void Polymer::init_atoms_rand()
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
    
    bonds(0,mu) = diffbasis(0)*u(0) + diffbasis(1)*v(0) + diffbasis(2)*w(0);
    bonds(1,mu) = diffbasis(0)*u(1) + diffbasis(1)*v(1) + diffbasis(2)*w(1);
    bonds(2,mu) = diffbasis(0)*u(2) + diffbasis(1)*v(2) + diffbasis(2)*w(2);
    
    u = bonds.col(mu);

  
  }

  xs.col(0) = x0;
;


  for (int mu = 0; mu< Nbeads-1; mu++) {
    xs.col(mu+1) = xs.col(mu) + bondlength*bonds.col(mu);
  }

  // shift so that centre of mass is halfway between the two end points specified

  Eigen::Vector3d com = {0,0,0};

  for (int i = 0; i < Nbeads; i++) {
    com += xs.col(i);
  }

  com /= Nbeads;

  for (int i = 0; i < Nbeads; i++) {
    xs.col(i) = xs.col(i) + (x0 + xN)/2 - com;
  }  
    
  return;

}

/* -------------------------------------------------------------------------- */
/* Call to set line configuration so that polymer intersects initial points
   x0 and xN with same centre of mass as the initial points. */
/* -------------------------------------------------------------------------- */
void Polymer::init_atoms_line()
{

  double length_init = (xN-x0).norm();
  
  double t1 = 0.5*(1-(Nbeads-1)*bondlength/length_init);
  for (int i = 0; i < Nbeads; i++) 
    xs.col(i) = x0 + (xN-x0)*(t1+i*bondlength/length_init);
    
  return;

}

void Polymer::init_atoms_caret()
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
      xs.col(i) = x0 + 2*alpha*tp;
    } else {
      xs.col(i) = xN + 2*beta*(1-tp);
    }
  }
  return;

}


  
};
