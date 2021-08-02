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
  
  rhs_of_G.resize(Nbeads-1);
  dummy_for_noise.resize(Nbeads-1);
  Gmunu.resize(Nbeads-1,Nbeads-1);
  Hhat.resize(Nbeads-1,Nbeads-1);
  rhs_of_Hhat.resize(Nbeads-1);
  
  tension.resize(Nbeads-1);


  tDets.resize(Nbeads);
  bDets.resize(Nbeads);

  costhetas.resize(Nbeads-2);
  k_effs.resize(Nbeads-2);

  std::cout << seed << std::endl;
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
  std::cout << dist(gen) << std::endl;
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

  for (int i = 0; i < Nbeads-1; i++) {
    rhs_of_G(i) = bonds[i].rod.dot(atoms[i+1].unprojected_noise
				   -atoms[i].unprojected_noise);

  }
  return;
  
}


/* -------------------------------------------------------------------------- */
/* Initialse G in G*eta = P (only call once). */
/* -------------------------------------------------------------------------- */
void Polymer::set_G()
{

  std::vector<T> coefficients = init_G_coeffsmatrix();

  Gmunu.setFromTriplets(coefficients.begin(),coefficients.end());

  solver.analyzePattern(Gmunu);

  return;

}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix G. */
/* -------------------------------------------------------------------------- */
std::vector<T> Polymer::init_G_coeffsmatrix()
{
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
void Polymer::update_G()
{
  // update first col by hand

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

void Polymer::compute_effective_kappa()
{

  tDets(0) = 1.0;
  bDets(Nbeads-1) = 1.0;

  tDets(1) = 2;
  bDets(Nbeads-2) = 2;

  int i,j2;
  for (i = 1; i < Nbeads - 1; i++ ) {

    j2 = Nbeads-1-i;

    tDets(i+1) = 2*tDets(i)-costhetas(i-1)*costhetas(i-1)*tDets(i-1);

    bDets(j2-1) = 2*bDets(j2)-costhetas(j2-1)*costhetas(j2-1)*bDets(j2+1);


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
void Polymer::set_rhs_of_Hhat()
{

  for (int i = 0; i< Nbeads-1; i++) {


    rhs_of_Hhat(i) = bonds[i].rod.dot(atoms[i+1].friction*atoms[i+1].F
				      -atoms[i].friction*atoms[i].F);
    
  }
  return;
  
}



/* -------------------------------------------------------------------------- */
/* Compute forces on the particles. */
/* -------------------------------------------------------------------------- */
void Polymer::compute_uc_forces()
{

  compute_effective_kappa();


  atoms[0].F = -1*k_effs(0)*(bonds[1].rod-costhetas(0)*bonds[0].rod) + atoms[0].noise;

  atoms[1].F = (k_effs(0)*(bonds[1].rod-costhetas(0)*bonds[0].rod
			   - (bonds[0].rod-costhetas(0)*bonds[1].rod))
		- k_effs(1)*(bonds[2].rod-costhetas(1)*bonds[1].rod)
		+ atoms[1].noise);

  			
  for (int k = 2; k < Nbeads-2; k++) {

    atoms[k].F = (k_effs(k-1)*(bonds[k].rod-costhetas(k-1)*bonds[k-1].rod
			       -(bonds[k-1].rod-costhetas(k-1)*bonds[k].rod))
		  -k_effs(k)*(bonds[k+1].rod-costhetas(k)*bonds[k].rod)
		  +k_effs(k-2)*(bonds[k-2].rod-costhetas(k-2)*bonds[k-1].rod)
		  +atoms[k].noise);
    
  }

  int k = Nbeads-2;
  
  atoms[k].F = (k_effs(k-1)*(bonds[k].rod-costhetas(k-1)*bonds[k-1].rod
			     -(bonds[k-1].rod-costhetas(k-1)*bonds[k].rod))
		+k_effs(k-2)*(bonds[k-2].rod-costhetas(k-2)*bonds[k-1].rod)
		+atoms[k].noise);

  k = Nbeads-1;

  atoms[k].F = (k_effs(k-2)*(bonds[k-2].rod-costhetas(k-2)*bonds[k-1].rod)
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

  solver.analyzePattern(Hhat);

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

/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix M. */
/* -------------------------------------------------------------------------- */
std::vector<T> Polymer::init_Hhat_coeffsmatrix()
{
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
void Polymer::update_Hhat()
{
  // update first col by hand
  Hhat.coeffRef(0,0) = Hhat_diag_val(0);

  Hhat.coeffRef(1,0) = Hhat_loweroff_val(1);

  Hhat.coeffRef(Nbeads-2,Nbeads-2) = Hhat_diag_val(Nbeads-2);


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
  return;
}




void Polymer::compute_noise()
{
  set_rhs_of_G();
  solver.factorize(Gmunu);

  dummy_for_noise =  solver.solve(rhs_of_G);

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


void Polymer::compute_tension()
{
  set_rhs_of_Hhat();
  solver.factorize(Hhat);

  tension =  solver.solve(rhs_of_Hhat);
  
}


void Polymer::initial_integrate(double Delta_t)
{
  
  double tmp = Delta_t/2.0;

  Rtmp[0] = atoms[0].R;

  atoms[0].F += tension(0)*bonds[0].rod;

  atoms[0].R += tmp*atoms[0].friction*atoms[0].F;

  
  for (int i = 1; i < Nbeads-1; i++) {

    Rtmp[i] = atoms[i].R;


    atoms[i].F += tension(i)*bonds[i].rod-tension(i-1)*bonds[i-1].rod;


    atoms[i].R += tmp*atoms[i].friction*atoms[i].F;

  }

  int i = Nbeads -1;

  Rtmp[i] = atoms[i].R;
  
  atoms[i].F += -tension(i-1)*bonds[i-1].rod;


  atoms[i].R += tmp*atoms[i].friction*atoms[i].F;  

  return;
}


void Polymer::final_integrate(double Delta_t)
{
  
  double tmp = Delta_t;

  atoms[0].F += tension(0)*bonds[0].rod;

  atoms[0].R = Rtmp[0] + tmp*atoms[0].friction*atoms[0].F;


    
  for (int i = 1; i < Nbeads-1; i++) {

    atoms[i].F += tension(i)*bonds[i].rod-tension(i-1)*bonds[i-1].rod;


    atoms[i].R = Rtmp[i] + tmp*atoms[i].friction*atoms[i].F;

    
  }

  int i = Nbeads -1;

  atoms[i].F += -tension(i-1)*bonds[i-1].rod;


  atoms[i].R = Rtmp[i] + tmp*atoms[i].friction*atoms[i].F;  

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
