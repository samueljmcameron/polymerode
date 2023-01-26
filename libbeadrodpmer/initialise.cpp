
#include <random>
#include <exception>

#include "no_tether.hpp"

#include "initialise.hpp"

#define SMALL 1e-14


namespace BeadRodPmer {
namespace Initialise {
/* -------------------------------------------------------------------------- */
/* Initialise a configuration (stored in xs vector) which relaxes from a
   caret (^) while holding the ends beads near endpoints (specified in istr),
   using the same parameters as whatever polymer is passed to it.
/* -------------------------------------------------------------------------- */

void init_atoms_relaxed_caret(const std::vector<std::string> &v_istr, const Polymer &pmer,
			      Eigen::Ref<Eigen::Matrix3Xd> xs)
{



  Eigen::Vector3d x0,xN;

  double springK,dt,tol;
  int seed;
  int nsteps,itermax;


  // read in input from istr
  if (v_istr.size() != 20)
    throw std::invalid_argument("incorrect number of arguments for polymer initialisation.");
  
  if (v_istr[0] != "x0")
    throw std::invalid_argument("missing x0 argument for polymer initialisation.");
  x0(0) = std::stod(v_istr[1]);
  x0(1) = std::stod(v_istr[2]);
  x0(2) = std::stod(v_istr[3]);

  if (v_istr[4] != "xN")
    throw std::invalid_argument("missing xN argument for polymer initialisation.");
  xN(0) = std::stod(v_istr[5]);
  xN(1) = std::stod(v_istr[6]);
  xN(2) = std::stod(v_istr[7]);

  if (v_istr[8] != "springK")
    throw std::invalid_argument("missing springK argument for polymer initialisation.");
  springK = std::stod(v_istr[9]);

  if (v_istr[10] != "dt")
    throw std::invalid_argument("missing dt argument for polymer initialisation.");
  dt = std::stod(v_istr[11]);


  if (v_istr[12] != "nsteps")
    throw std::invalid_argument("missing nsteps argument for polymer initialisation.");
  nsteps = std::stoi(v_istr[13]);

  if (v_istr[14] != "tol")
    throw std::invalid_argument("missing tol argument for polymer initialisation.");
  tol = std::stod(v_istr[15]);


  if (v_istr[16] != "itermax")
    throw std::invalid_argument("missing tol argument for polymer initialisation.");
  itermax = std::stoi(v_istr[17]);

  if (v_istr[18] != "seed")
    throw std::invalid_argument("missing tol argument for polymer initialisation.");
  seed = std::stoi(v_istr[19]);
  
  
  // transfer arguments over to a no-tethered polymer (except seed)
  std::vector<std::string> v_line
    = {"beads",std::to_string(pmer.get_Nbeads()),
       "bondlength",std::to_string(pmer.get_bondlength()),
       "zeta_para",std::to_string(pmer.get_zpara()),
       "zeta_perp",std::to_string(pmer.get_zperp()),
       "bending",std::to_string(pmer.get_bending()),
       "temp",std::to_string(pmer.get_temp()),
       "seed",std::to_string(seed)};

  
  NoTether NoTeth(v_line);

  // using the same seed may be problematic in some cases.. probably not as its only for
  //  initialisation?
  init_atoms_caret(xs,x0,xN,NoTeth.get_bondlength(),seed);
  
  Eigen::Matrix3Xd Fs(3,NoTeth.get_Nbeads());

  int last_bead = Fs.cols()-1;

  NoTeth.setup(xs);

  double t = 0;

  int iterations;

  Fs.setZero();
  Fs.col(0) = -springK*(xs.col(0)-x0);
  Fs.col(last_bead) = -springK*(xs.col(last_bead)-xN);
  
  NoTeth.first_step(xs,Fs,dt);

  Fs.setZero();
  Fs.col(0) = -springK*(xs.col(0)-x0);
  Fs.col(last_bead) = -springK*(xs.col(last_bead)-xN);

  iterations = NoTeth.second_step(xs,Fs,dt,itermax);

  if (iterations > itermax)
      throw std::runtime_error("Unable to initialise polymer (failure at time "
			       + std::to_string(t)
			       + std::string("). Consider reducing timestep (currently ")
			       + std::to_string(dt)
			       + std::string(") or maybe increasing itermax (currently ")
			       + std::to_string(itermax) + std::string(")."));

  t += dt;


  while ((xs.col(last_bead)-xN).norm() > tol
	 || (xs.col(0)-x0).norm() > tol
	 || t < nsteps*dt) {


    Fs.setZero();
    Fs.col(0) = -springK*(xs.col(0)-x0);
    Fs.col(last_bead) = -springK*(xs.col(last_bead)-xN);
    
    NoTeth.first_step(xs,Fs,dt);
    
    Fs.setZero();
    Fs.col(0) = -springK*(xs.col(0)-x0);
    Fs.col(last_bead) = -springK*(xs.col(last_bead)-xN);
    
    iterations = NoTeth.second_step(xs,Fs,dt,itermax);
    
    if (iterations > itermax)
      throw std::runtime_error("Unable to initialise polymer (failure at time "
			       + std::to_string(t)
			       + std::string("). Consider reducing timestep (currently ")
			       + std::to_string(dt)
			       + std::string(") or maybe increasing itermax (currently ")
			       + std::to_string(itermax) + std::string(")."));
    
    
    t += dt;


  }
  
  return ;
  
  
  
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
/* Call to set the xs vector to a random configuration with centre of mass
   given by (x0 + xN)/2, and the distance between points in xs is given by
   bondlength. */
/* -------------------------------------------------------------------------- */
void init_atoms_rand(Eigen::Ref<Eigen::Matrix3Xd> xs,
		     Eigen::Vector3d x0,Eigen::Vector3d xN,
		     double bondlength, int seed)
{

  Eigen::Vector3d u(1,0,0);
  Eigen::Vector3d v(0,1,0);
  Eigen::Vector3d w(0,0,1);

  Eigen::Vector3d diffbasis(0,0,0);


  std::mt19937 gen;

  std::uniform_real_distribution<double> dist;

  gen.seed(seed);

  
  
  double tmp,swich;
  double dum;

  int Nbeads = xs.cols();

  Eigen::Matrix3Xd bonds(3,Nbeads);

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
/* Call to set line configuration so that polymer intersects points
   x0 and xN with same centre of mass as the initial points.  */
/* -------------------------------------------------------------------------- */
void init_atoms_line(Eigen::Ref<Eigen::Matrix3Xd> xs,
		     Eigen::Vector3d x0,Eigen::Vector3d xN,
		     double bondlength)
{

  
  double length_init = (xN-x0).norm();
  int Nbeads = xs.cols();
  
  double t1 = 0.5*(1-(Nbeads-1)*bondlength/length_init);
  for (int i = 0; i < Nbeads; i++) 
    xs.col(i) = x0 + (xN-x0)*(t1+i*bondlength/length_init);
    
  return;

}
/* -------------------------------------------------------------------------- */
/* Call to set the xs vector to a caret (^) configuration with centre of mass
   given by (x0 + xN)/2, and the distance between points in xs is given by
   bondlength. */
/* -------------------------------------------------------------------------- */
void init_atoms_caret(Eigen::Ref<Eigen::Matrix3Xd> xs,
		      Eigen::Vector3d x0,Eigen::Vector3d xN,
		      double bondlength,int seed)
{

  std::mt19937 gen;

  std::uniform_real_distribution<double> dist;

  gen.seed(seed);

  int Nbeads = xs.cols();

  Eigen::Matrix3Xd bonds(3,Nbeads);


  

  
  Eigen::Vector3d nhat;

  nhat(0) = 2*dist(gen);
  nhat(1) = 2*dist(gen);
  nhat(2) = 2*dist(gen);
  
  Eigen::Vector3d dd = xN-x0;

  if (dd.norm() > bondlength*(xs.cols()-1)) {
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
  if (xs.cols() % 2 == 0) {
    fakelength = (xs.cols()-2)*bondlength;
  } else {
    fakelength = (xs.cols()-1)*bondlength;
  }
  

  double gamma = 0.5*sqrt(fakelength*fakelength-dd.squaredNorm());

  Eigen::Vector3d alpha = 0.5*dd+gamma*nhat;
  Eigen::Vector3d beta = -0.5*dd+gamma*nhat;


  double tp = 0;
  double dtp = 1.0*bondlength/fakelength;

  
  for (int i = 0; i < xs.cols(); i++) {
    
    tp = i*dtp;

    // integer division rounds toward zero
    if (i < xs.cols()/2) {
      xs.col(i) = x0 + 2*alpha*tp;
    } else {
      xs.col(i) = xN + 2*beta*(1-tp);
    }
  }
  return;

}


  
}
}
