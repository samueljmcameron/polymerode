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
