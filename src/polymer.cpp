#include "polymer.hpp"
#include "input.hpp"

#include <iostream>
#include <cmath>
#include <stdexcept>


#define SMALL 1e-14
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
Polymer::Polymer(std::vector<std::string> splitvec)
  : gen(0), dist(-0.5,0.5)
{


  int nargs = splitvec.size();
  int iarg = 0;

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
    } else {
      throw std::runtime_error("Error: invalid argument for build_polymer.");
    }
  }
    
  atoms = std::vector<Atom> (Nbeads);

  bonds = std::vector<Bond>(Nbeads-1);
  
  Rtmp.resize(Nbeads);
  
  rhs_of_G.resize(Nbeads+5);
  dummy_for_noise.resize(Nbeads+5);
  Gmunu.resize(Nbeads+5,Nbeads+5);
  Hhat.resize(Nbeads+5,Nbeads+5);
  rhs_of_Hhat.resize(Nbeads+5);
  
  tension.resize(Nbeads+5);


  tDets.resize(Nbeads+3);
  bDets.resize(Nbeads+3);

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


/* ---------------------------------------------------------------------------- */
/* Compute the unit vectors tangent to the bonds (u) and beads (tu). */
/* ---------------------------------------------------------------------------- */
void Polymer::compute_tangents_and_rods_and_friction()
{

  double unorm;
  double tangentnorm;


  unorm = sqrt((atoms[1].R-atoms[0].R).dot(atoms[1].R-atoms[0].R));
  //unorm = bondlength;

  bonds[0].rod = (atoms[1].R - atoms[0].R)/unorm;

  tangentnorm = sqrt(bonds[0].rod.dot(bonds[0].rod));
  atoms[0].tangent = bonds[0].rod/tangentnorm;


  single_inv_friction(0);


  for (int i = 1; i < Nbeads-1; i++) {

    unorm = sqrt((atoms[i+1].R-atoms[i].R).dot(atoms[i+1].R-atoms[i].R));


    bonds[i].rod = (atoms[i+1].R - atoms[i].R)/unorm;

    costhetas(i-1) = bonds[i].rod.dot(bonds[i-1].rod);
      
    tangentnorm = sqrt((bonds[i].rod+bonds[i-1].rod
			).dot(bonds[i].rod+bonds[i-1].rod));


    atoms[i].tangent = (bonds[i].rod+bonds[i-1].rod)/tangentnorm;

    single_inv_friction(i);

  }
  

  tangentnorm = sqrt(bonds[Nbeads-2].rod.dot(bonds[Nbeads-2].rod));
  atoms[Nbeads-1].tangent = bonds[Nbeads-2].rod;

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


/* ---------------------------------------------------------------------------- */
/* RHS of G*eta = P. */
/* ---------------------------------------------------------------------------- */
void Polymer::set_rhs_of_G()
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
void Polymer::set_G()
{

  std::vector<T> coefficients = init_G_coeffsmatrix();

  Gmunu.setFromTriplets(coefficients.begin(),coefficients.end());

  Gmunu_solver.analyzePattern(Gmunu);

  return;

}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix G. */
/* -------------------------------------------------------------------------- */
std::vector<T> Polymer::init_G_coeffsmatrix()
{
  // just set lower diagonal for some reason..?
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
void Polymer::update_G()
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

void Polymer::compute_effective_kappa()
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
void Polymer::set_rhs_of_Hhat()
{


  int offset = 2;

  Eigen::Vector3d tmp = atoms[0].friction*atoms[0].F;

  rhs_of_Hhat(-2+offset) = tmp(0);
  rhs_of_Hhat(-1+offset) = tmp(1);
  rhs_of_Hhat(0+offset) = tmp(2);
  
  for (int i = 0; i< Nbeads-1; i++) {


    rhs_of_Hhat(i+offset+1) = bonds[i].rod.dot(atoms[i+1].friction*atoms[i+1].F
					     -atoms[i].friction*atoms[i].F);
    
  }

  tmp = atoms[Nbeads-1].friction*atoms[Nbeads-1].F;
  rhs_of_Hhat(offset+Nbeads) = tmp(0);
  rhs_of_Hhat(offset+Nbeads+1) = tmp(1);
  rhs_of_Hhat(offset+Nbeads+2) = tmp(2);
  
  return;
  
}



/* -------------------------------------------------------------------------- */
/* Compute forces on the particles. */
/* -------------------------------------------------------------------------- */
void Polymer::compute_uc_forces()
{

  compute_effective_kappa();

  atoms[0].F = (end_inverses(0)*bonds[0].rod*(bonds[0].rod(0))
		+end_inverses(1)*bonds[0].rod*(bonds[0].rod(1))
		+end_inverses(2)*bonds[0].rod*(bonds[0].rod(2))
		-1*k_effs(0)*(bonds[1].rod-costhetas(0)*bonds[0].rod)
		+ atoms[0].noise);

  atoms[0].F(0) += -1*end_inverses(0);
  atoms[0].F(1) += -1*end_inverses(1);
  atoms[0].F(2) += -1*end_inverses(2);

  atoms[1].F = (-1*end_inverses(0)*bonds[0].rod*(bonds[0].rod(0))
		-end_inverses(1)*bonds[0].rod*(bonds[0].rod(1))
		-end_inverses(2)*bonds[0].rod*(bonds[0].rod(2))
		+k_effs(0)*(bonds[1].rod-costhetas(0)*bonds[0].rod
			    - (bonds[0].rod-costhetas(0)*bonds[1].rod))
		- k_effs(1)*(bonds[2].rod-costhetas(1)*bonds[1].rod)
		+ atoms[1].noise);


  atoms[0].F(0) += end_inverses(0);
  atoms[0].F(1) += end_inverses(1);
  atoms[0].F(2) += end_inverses(2);

  

  			
  for (int k = 2; k < Nbeads-2; k++) {

    atoms[k].F = (k_effs(k-1)*(bonds[k].rod-costhetas(k-1)*bonds[k-1].rod
			       -(bonds[k-1].rod-costhetas(k-1)*bonds[k].rod))
		  -k_effs(k)*(bonds[k+1].rod-costhetas(k)*bonds[k].rod)
		  +k_effs(k-2)*(bonds[k-2].rod-costhetas(k-2)*bonds[k-1].rod)
		  +atoms[k].noise);
    
  }

  int k = Nbeads-2;
  
  atoms[k].F = (-1*end_inverses(3)*bonds[k].rod*bonds[k].rod(0)
		-end_inverses(4)*bonds[k].rod*bonds[k].rod(1)
		-end_inverses(5)*bonds[k].rod*bonds[k].rod(2)
		+k_effs(k-1)*(bonds[k].rod-costhetas(k-1)*bonds[k-1].rod
			     -(bonds[k-1].rod-costhetas(k-1)*bonds[k].rod))
		+k_effs(k-2)*(bonds[k-2].rod-costhetas(k-2)*bonds[k-1].rod)
		+atoms[k].noise);

  atoms[k].F(0) += end_inverses(3);
  atoms[k].F(1) += end_inverses(4);
  atoms[k].F(2) += end_inverses(5);

  k = Nbeads-1;

  atoms[k].F = (end_inverses(3)*bonds[k-1].rod*bonds[k-1].rod(0)
		+end_inverses(4)*bonds[k-1].rod*bonds[k-1].rod(1)
		+end_inverses(5)*bonds[k-1].rod*bonds[k-1].rod(2)
		+k_effs(k-2)*(bonds[k-2].rod-costhetas(k-2)*bonds[k-1].rod)
		+atoms[k].noise);
  
  
  return;
  
}



/* -------------------------------------------------------------------------- */
/* Initialise matrix H hat. */
/* -------------------------------------------------------------------------- */
void Polymer::set_Hhat()
{

  std::vector<T> coefficients = init_Hhat_coeffsmatrix();

  Hhat.setFromTriplets(coefficients.begin(),coefficients.end());

  Hhat_solver.analyzePattern(Hhat);

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

double Polymer::Hhat_bottomside(int first)
{

  double tmp =  atoms[Nbeads-1].tangent.dot(bonds[Nbeads-2].rod)*atoms[Nbeads-1].tangent(first);

  
  return (1./zpara-1./zperp)*tmp + bonds[Nbeads-2].rod(first)/zperp;

}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix M. */
/* -------------------------------------------------------------------------- */
std::vector<T> Polymer::init_Hhat_coeffsmatrix()
{
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
void Polymer::update_Hhat()
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




void Polymer::compute_noise()
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
		      + dummy_for_noise(i)*bonds[i].rod
		      - dummy_for_noise(i-1)*bonds[i-1].rod);
  }


  i = Nbeads-1;

  atoms[i].noise = (atoms[i].unprojected_noise
		    - dummy_for_noise(i-1)*bonds[i-1].rod);

  atoms[i].noise(0) -= dummy_for_noise(offset+1+Nbeads-1);
  atoms[i].noise(1) -= dummy_for_noise(offset+1+Nbeads);
  atoms[i].noise(2) -= dummy_for_noise(offset+1+Nbeads+1);
  
}


void Polymer::compute_tension()
{
  set_rhs_of_Hhat();
  Hhat_solver.factorize(Hhat);

  tension =  Hhat_solver.solve(rhs_of_Hhat);
  
}


void Polymer::initial_integrate(double Delta_t)
{
  
  double tmp = Delta_t/2.0;
  int offset = 2;
  
  Rtmp[0] = atoms[0].R;

  atoms[0].F += tension(offset+1)*bonds[0].rod;

  atoms[0].F(0) -= tension(offset-2);
  atoms[0].F(1) -= tension(offset-1);
  atoms[0].F(2) -= tension(offset);

  atoms[0].R += tmp*atoms[0].friction*atoms[0].F;
  
  for (int i = 1; i < Nbeads-1; i++) {

    Rtmp[i] = atoms[i].R;


    atoms[i].F += tension(offset+1+i)*bonds[i].rod-tension(offset+i)*bonds[i-1].rod;


    atoms[i].R += tmp*atoms[i].friction*atoms[i].F;
  }

  int i = Nbeads -1;

  Rtmp[i] = atoms[i].R;
  
  atoms[i].F += -tension(offset+i)*bonds[i-1].rod;

  atoms[i].F(0) -= tension(offset+Nbeads);
  atoms[i].F(1) -= tension(offset+Nbeads+1);
  atoms[i].F(2) -= tension(offset+Nbeads+2);


  atoms[i].R += tmp*atoms[i].friction*atoms[i].F;  

  
  return;
}


void Polymer::final_integrate(double Delta_t)
{
  
  double tmp = Delta_t;

  int offset=2;

  
  atoms[0].F += tension(offset+1)*bonds[0].rod;

  atoms[0].F(0) -= tension(offset-2);
  atoms[0].F(1) -= tension(offset-1);
  atoms[0].F(2) -= tension(offset);
  

  atoms[0].R = Rtmp[0] + tmp*atoms[0].friction*atoms[0].F;


    
  for (int i = 1; i < Nbeads-1; i++) {

    atoms[i].F += tension(offset+1+i)*bonds[i].rod-tension(offset+i)*bonds[i-1].rod;


    atoms[i].R = Rtmp[i] + tmp*atoms[i].friction*atoms[i].F;

    
    
  }

  int i = Nbeads -1;

  atoms[i].F += -tension(offset+i)*bonds[i-1].rod;

  atoms[i].F(0) -= tension(offset+Nbeads);
  atoms[i].F(1) -= tension(offset+Nbeads+1);
  atoms[i].F(2) -= tension(offset+Nbeads+2);

  atoms[i].R = Rtmp[i] + tmp*atoms[i].friction*atoms[i].F;


  for (int mu = 0; mu < Nbeads -1; mu ++) {
    std::cout << "mu = " << mu << ", constraint is not zero but "
	      << bonds[mu].rod.dot(atoms[mu+1].R-atoms[mu].R) << std::endl;
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

int Polymer::get_Nbeads()
{
  return Nbeads;
}


double Polymer::get_temp()
{
  return temp;
}

double Polymer::get_zpara()
{
  return zpara;
}

double Polymer::get_zperp()
{
  return zperp;
}

double Polymer::get_bondlength()
{
  return bondlength;
}

double Polymer::get_timescale(double dt)
{
  return dt*temp/(bondlength*zperp);
}
