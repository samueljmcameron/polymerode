#include "polymer.hpp"
#include "input.hpp"


#include <iostream>
#include <cmath>
#include <stdexcept>


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

  
  
    
  atoms = std::vector<Atom> (Nbeads);
  bonds = std::vector<Bond>(Nbeads-1);
  tmpbonds = std::vector<Bond>(Nbeads-1);
  Rtmp.resize(Nbeads);
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


  bonds[0].rod = (atoms[1].R-atoms[0].R)/bondlength;
  
  tangentnorm = sqrt(bonds[0].rod.dot(bonds[0].rod));
  atoms[0].tangent = bonds[0].rod/tangentnorm;


  single_inv_friction(0);


  for (int i = 1; i < Nbeads-1; i++) {

    bonds[i].rod = (atoms[i+1].R-atoms[i].R)/bondlength;
    costhetas(i-1) = bonds[i].rod.dot(bonds[i-1].rod);
      
    tangentnorm = sqrt((bonds[i].rod+bonds[i-1].rod
			).dot(bonds[i].rod+bonds[i-1].rod));


    atoms[i].tangent = (bonds[i].rod+bonds[i-1].rod)/tangentnorm;

    single_inv_friction(i);

  }
  

  tangentnorm = sqrt(bonds[Nbeads-2].rod.dot(bonds[Nbeads-2].rod));
  atoms[Nbeads-1].tangent = bonds[Nbeads-2].rod/tangentnorm;

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
    
    atoms[i].unprojected_noise = dum*(perp*rands
				      +(para-perp)*(atoms[i].tangent.dot(rands))
				      *atoms[i].tangent);

  }

  return;
}






/* -------------------------------------------------------------------------- */
/* Compute forces on the particles. */
/* -------------------------------------------------------------------------- */
void Polymer::compute_uc_forces()
{

  atoms[0].Fpot = (end_inverses(0)*bonds[0].rod*(bonds[0].rod(0))
		   +end_inverses(1)*bonds[0].rod*(bonds[0].rod(1))
		   +end_inverses(2)*bonds[0].rod*(bonds[0].rod(2))
		   -1*k_effs(0)*(bonds[1].rod-costhetas(0)*bonds[0].rod));

  atoms[0].Fpot(0) += -1*end_inverses(0);
  atoms[0].Fpot(1) += -1*end_inverses(1);
  atoms[0].Fpot(2) += -1*end_inverses(2);

  atoms[1].Fpot = (-1*end_inverses(0)*bonds[0].rod*(bonds[0].rod(0))
		-end_inverses(1)*bonds[0].rod*(bonds[0].rod(1))
		-end_inverses(2)*bonds[0].rod*(bonds[0].rod(2))
		+k_effs(0)*(bonds[1].rod-costhetas(0)*bonds[0].rod
			    - (bonds[0].rod-costhetas(0)*bonds[1].rod))
		- k_effs(1)*(bonds[2].rod-costhetas(1)*bonds[1].rod));


  atoms[1].Fpot(0) += end_inverses(0);
  atoms[1].Fpot(1) += end_inverses(1);
  atoms[1].Fpot(2) += end_inverses(2);

  

  			
  for (int k = 2; k < Nbeads-2; k++) {

    atoms[k].Fpot = (k_effs(k-1)*(bonds[k].rod-costhetas(k-1)*bonds[k-1].rod
				  -(bonds[k-1].rod-costhetas(k-1)*bonds[k].rod))
		     -k_effs(k)*(bonds[k+1].rod-costhetas(k)*bonds[k].rod)
		     +k_effs(k-2)*(bonds[k-2].rod-costhetas(k-2)*bonds[k-1].rod));
    
  }

  int k = Nbeads-2;
  
  atoms[k].Fpot = (-1*end_inverses(3)*bonds[k].rod*bonds[k].rod(0)
		   -end_inverses(4)*bonds[k].rod*bonds[k].rod(1)
		   -end_inverses(5)*bonds[k].rod*bonds[k].rod(2)
		   +k_effs(k-1)*(bonds[k].rod-costhetas(k-1)*bonds[k-1].rod
				 -(bonds[k-1].rod-costhetas(k-1)*bonds[k].rod))
		   +k_effs(k-2)*(bonds[k-2].rod-costhetas(k-2)*bonds[k-1].rod));

  atoms[k].Fpot(0) += end_inverses(3);
  atoms[k].Fpot(1) += end_inverses(4);
  atoms[k].Fpot(2) += end_inverses(5);

  k = Nbeads-1;

  atoms[k].Fpot = (end_inverses(3)*bonds[k-1].rod*bonds[k-1].rod(0)
		   +end_inverses(4)*bonds[k-1].rod*bonds[k-1].rod(1)
		   +end_inverses(5)*bonds[k-1].rod*bonds[k-1].rod(2)
		   +k_effs(k-2)*(bonds[k-2].rod-costhetas(k-2)*bonds[k-1].rod));


  atoms[k].Fpot(0) -= end_inverses(3);
  atoms[k].Fpot(1) -= end_inverses(4);
  atoms[k].Fpot(2) -= end_inverses(5);
  

  return;
  
}


void Polymer::single_inv_friction(int i)
{
  double tmp = 1./zpara-1./zperp;
  atoms[i].friction = tmp*atoms[i].tangent*atoms[i].tangent.transpose();
  
  for (int j = 0; j < 3; j++ )

    atoms[i].friction(j,j) += 1./zperp;

  return;

}

void Polymer::add_external_force(const std::vector<double> &dFdX, int i)
{
  atoms[i].Fpot(0) += -dFdX[0];
  atoms[i].Fpot(1) += -dFdX[1];
  atoms[i].Fpot(2) += -dFdX[2];
  return;

}

double Polymer::Hhat_diag_val(int mu)
{


  double tmp1 = atoms[mu].tangent.dot(bonds[mu].rod);
  tmp1 = tmp1*tmp1;


  double tmp2 =  atoms[mu+1].tangent.dot(bonds[mu].rod);
  tmp2 = tmp2*tmp2;


  return 2./zperp + (1./zpara-1./zperp)*(tmp1+tmp2);

}



double Polymer::dCdlambda_diag_val(int mu)
{


  double tmp1 = atoms[mu].tangent.dot(tmpbonds[mu].rod);
  tmp1 = tmp1*atoms[mu].tangent.dot(bonds[mu].rod);


  double tmp2 =  atoms[mu+1].tangent.dot(tmpbonds[mu].rod);
  tmp2 = tmp2*atoms[mu+1].tangent.dot(bonds[mu].rod);

  double tmp3 = tmpbonds[mu].rod.dot(bonds[mu].rod);


  return (2*tmp3/zperp + (1./zpara-1./zperp)*(tmp1+tmp2))/tmpbonds[mu].rod.norm();

}


/* NOte that there is a typo in the paper I'm basing this off of
   (Montesi et al, 2005) , which is why there is a difference between
   my off-diagonal H values and theirs. In their equation 34, anywhere
   there is a mu+1, it should be a mu-1. */
double Polymer::Hhat_loweroff_val(int mu)
{

  double tmp1 = bonds[mu-1].rod.dot(bonds[mu].rod);
  double tmp2 = bonds[mu-1].rod.dot(atoms[mu].tangent);
  double tmp3 = atoms[mu].tangent.dot(bonds[mu].rod);

  return -1./zperp*tmp1-(1./zpara-1./zperp)*tmp2*tmp3;

}

double Polymer::dCdlambda_loweroff_val(int mu)
{

  double tmp1 = bonds[mu-1].rod.dot(tmpbonds[mu].rod);
  double tmp2 = bonds[mu-1].rod.dot(atoms[mu].tangent);
  double tmp3 = atoms[mu].tangent.dot(tmpbonds[mu].rod);

  return (-1./zperp*tmp1-(1./zpara-1./zperp)*tmp2*tmp3)/tmpbonds[mu].rod.norm();

}


double Polymer::dCdlambda_upperoff_val(int mu)
{

  double tmp1 = bonds[mu].rod.dot(tmpbonds[mu-1].rod);
  double tmp2 = tmpbonds[mu-1].rod.dot(atoms[mu].tangent);
  double tmp3 = atoms[mu].tangent.dot(bonds[mu].rod);

  return (-1./zperp*tmp1-(1./zpara-1./zperp)*tmp2*tmp3)/tmpbonds[mu-1].rod.norm();

}



double Polymer::Hhat_endblocks(int first, int second,int mu)
{
  double extra = 0;

  if (first == second) extra = 1;
  
  return (extra/zperp + (1./zpara-1./zperp)*atoms[mu].tangent(first)
	  *atoms[mu].tangent(second));

}


double Polymer::Hhat_leftside(int first)
{

  double tmp =  atoms[0].tangent.dot(bonds[0].rod)*atoms[0].tangent(first);

  
  return -1*(1./zpara-1./zperp)*tmp - bonds[0].rod(first)/zperp;

}


double Polymer::dCdlambda_leftside(int first)
{

  double tmp =  atoms[0].tangent.dot(tmpbonds[0].rod)*atoms[0].tangent(first);

  
  return (-1*(1./zpara-1./zperp)*tmp - tmpbonds[0].rod(first)/zperp)/tmpbonds[0].rod.norm();

}



double Polymer::Hhat_bottomside(int first)
{
  
  double tmp =  atoms[Nbeads-1].tangent.dot(bonds[Nbeads-2].rod)*atoms[Nbeads-1].tangent(first);
  
  
  return (1./zpara-1./zperp)*tmp + bonds[Nbeads-2].rod(first)/zperp;
  
}


double Polymer::dCdlambda_bottomside(int first)
{
  
  double tmp =  atoms[Nbeads-1].tangent.dot(tmpbonds[Nbeads-2].rod)
    *atoms[Nbeads-1].tangent(first);
  
  
  return ((1./zpara-1./zperp)*tmp + tmpbonds[Nbeads-2].rod(first)/zperp)
    /tmpbonds[Nbeads-2].rod.norm();
  
}
  
  
/* ---------------------------------------------------------------------------- */
/* RHS of G*eta = P. To be used in derived classes. */
/* ---------------------------------------------------------------------------- */
void Polymer::set_rhs_of_G(int offset)
{
  
  
  
  for (int i = 0; i < Nbeads-1; i++) {
    rhs_of_G(i+offset) = bonds[i].rod.dot(atoms[i+1].unprojected_noise
					  -atoms[i].unprojected_noise);
    
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


    rhs_of_Hhat(i+offset) = bonds[i].rod.dot(atoms[i+1].friction*(atoms[i+1].Fpot
								  +atoms[i+1].noise)
					     -atoms[i].friction*(atoms[i].Fpot
								 +atoms[i].noise));
    
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
  atoms[i].noise = (atoms[i].unprojected_noise
		    + dummy_for_noise(i+offset)*bonds[i].rod);

  
  
  for (i = 1; i < Nbeads-1; i++) {
    atoms[i].noise = (atoms[i].unprojected_noise
		      + dummy_for_noise(i+offset)*bonds[i].rod
		      - dummy_for_noise(i+offset-1)*bonds[i-1].rod);
  }


  i = Nbeads-1;

  atoms[i].noise = (atoms[i].unprojected_noise
		    - dummy_for_noise(i+offset-1)*bonds[i-1].rod);
  return;

}





void Polymer::initial_integrate(double Delta_t,int offset, PTYPE cflag)
{
  
  double tmp = Delta_t/2.0;

  double tangentnorm;

  
  int i = 0;
  Rtmp[0] = atoms[0].R;

  atoms[0].t_force = tension(offset)*bonds[0].rod;

  if (cflag == SINGLE || cflag == DOUBLE) atoms[0].t_force -= tension({0,1,2});

  atoms[0].R += tmp*atoms[0].friction*(atoms[0].Fpot+atoms[0].noise+atoms[0].t_force);

  i = 1;
  Rtmp[i] = atoms[i].R;


  atoms[i].t_force = tension(offset+i)*bonds[i].rod-tension(offset+i-1)*bonds[i-1].rod;
  
  
  atoms[i].R += tmp*atoms[i].friction*(atoms[i].Fpot+atoms[i].noise+atoms[i].t_force);

  bonds[i-1].rod = (atoms[i].R-atoms[i-1].R)/bondlength;
  tmpbonds[i-1].rod = bonds[i-1].rod;
  
  tangentnorm = sqrt(bonds[i-1].rod.dot(bonds[i-1].rod));
  atoms[i-1].tangent = bonds[i-1].rod/tangentnorm;


  // compute friction with new values of R on bead i-1
  single_inv_friction(i-1);
  
  for (i = 2; i < Nbeads-1; i++) {

    Rtmp[i] = atoms[i].R;

    atoms[i].t_force = tension(offset+i)*bonds[i].rod-tension(offset-1+i)*bonds[i-1].rod;


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

  atoms[i].t_force = -tension(offset-1+i)*bonds[i-1].rod;

  if (cflag == DOUBLE) {
    int ind = offset+Nbeads-1;
    atoms[i].t_force -= tension({ind,ind+1,ind+2});
  }


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
void Polymer::final_integrate(double Delta_t, int offset,PTYPE cflag)
{
  
  double tmp = Delta_t;


  atoms[0].t_force= tension(offset)*bonds[0].rod;

  if (cflag == SINGLE || cflag == DOUBLE) {
    atoms[0].t_force -= tension({0,1,2});
  }

  atoms[0].R = Rtmp[0] + tmp*atoms[0].friction*(atoms[0].Fpot+atoms[0].noise+atoms[0].t_force);


    
  for (int i = 1; i < Nbeads-1; i++) {

    atoms[i].t_force = tension(offset+i)*bonds[i].rod-tension(offset-1+i)*bonds[i-1].rod;

    atoms[i].R = Rtmp[i] + tmp*atoms[i].friction*(atoms[i].Fpot+atoms[i].noise+atoms[i].t_force);
    tmpbonds[i-1].rod = (atoms[i].R-atoms[i-1].R)/bondlength;
    
  }

  int i = Nbeads -1;

  atoms[i].t_force = -tension(offset-1+i)*bonds[i-1].rod;

  if (cflag == DOUBLE) {
    int ind = offset+Nbeads-1;
    atoms[i].t_force -= tension({ind,ind+1,ind+2});
  }

  atoms[i].R = Rtmp[i] + tmp*atoms[i].friction*(atoms[i].Fpot+atoms[i].noise+atoms[i].t_force);
  tmpbonds[i-1].rod = (atoms[i].R-atoms[i-1].R)/bondlength;

  
  return;
}


  
  
  
void Polymer::calculate_constraint_errors(int offset)
{
  for (int mu = 1; mu < Nbeads; mu++) {
    constraint_errors(mu-1+offset) = (atoms[mu].R-atoms[mu-1].R).norm()-bondlength;
  }
  
  return;
  
}

void Polymer::rescale_positions(bool fromfront)
{
  
  if (fromfront) {
    for (int i = 1; i < Nbeads; i++) {
      atoms[i].R = atoms[i-1].R + bondlength*bonds[i-1].rod;
    }
  } else {
    for (int i = Nbeads-2; i >= 0; i--) {
      atoms[i].R = atoms[i+1].R - bondlength*bonds[i-1].rod;
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
};
