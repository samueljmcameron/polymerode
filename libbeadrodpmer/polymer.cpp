#include "polymer.hpp"
#include "input.hpp"


#include <iostream>
#include <cmath>
#include <stdexcept>


#define SMALL 1e-14

namespace BeadRodPmer {
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
/* Constructor */
/* -------------------------------------------------------------------------- */
Polymer::Polymer(const std::vector<std::string> & splitvec)
  : gen(0), dist(-0.5,0.5)
{


  int nargs = splitvec.size();
  int iarg = 0;

  int number_of_nuc_beads;
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
    } else if (splitvec[iarg] == "initialise") {
      input::isDouble(splitvec[iarg+1],x0(0),splitvec[iarg]);
      input::isDouble(splitvec[iarg+2],x0(1),splitvec[iarg]);
      input::isDouble(splitvec[iarg+3],x0(2),splitvec[iarg]);
      iarg += 4;
    } else if (splitvec[iarg] == "seed") {
      input::isInt(splitvec[iarg+1],seed,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "nuc_beads") {
      input::isInt(splitvec[iarg+1],number_of_nuc_beads,splitvec[iarg]);
      nuc_beads.resize(number_of_nuc_beads);
      for (int i = 1; i <= number_of_nuc_beads; i++) 
	input::isInt(splitvec[iarg+1+i],nuc_beads[i-1],splitvec[iarg]);
      iarg += number_of_nuc_beads+2;
    } else {
      throw std::runtime_error("Error: invalid argument for build_polymer.");
    }
  }
    
  atoms = std::vector<Atom> (Nbeads);
  bonds = std::vector<Bond>(Nbeads-1);
  tmpbonds = std::vector<Bond>(Nbeads-1);
  Rtmp.resize(Nbeads);
  costhetas.resize(Nbeads-2);
  k_effs.resize(Nbeads-2);

  end_inverses.resize(6);
  gen.seed(seed);

  init_atoms();

  
}

/* -------------------------------------------------------------------------- */
/* Destructor */
/* -------------------------------------------------------------------------- */
Polymer::~Polymer()
{
}


/* -------------------------------------------------------------------------- */
/* Call to set configuration to a straight rod along x axis. */
/* -------------------------------------------------------------------------- */
void Polymer::init_atoms()
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

  xN(0) = atoms[Nbeads-1].R(0);
  xN(1) = atoms[Nbeads-1].R(1);
  xN(2) = atoms[Nbeads-1].R(2);
    
  return;

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

};

