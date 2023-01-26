#include "polymer.hpp"


#include <iostream>
#include <cmath>
#include <stdexcept>



namespace BeadRodPmer {
/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
Polymer::Polymer(const std::vector<std::string> & v_line)
  : gen(0), dist(-0.5,0.5)
{

  if (v_line.size() != 14)
    throw std::invalid_argument("incorrect number of arguments for polymer/no_tether.");


  if (v_line[0] != "beads")
    throw std::invalid_argument("missing bondlength argument for polymer/no_tether.");
  Nbeads = std::stoi(v_line[1]);

  if (v_line[2] != "bondlength")
    throw std::invalid_argument("missing bondlength argument for polymer/no_tether.");
  bondlength = std::stod(v_line[3]);

  if (v_line[4] != "zeta_para")
    throw std::invalid_argument("missing zeta_para argument for polymer/no_tether.");
  zpara = std::stod(v_line[5]);


  if (v_line[6] != "zeta_perp")
    throw std::invalid_argument("missing zeta_perp argument for polymer/no_tether.");
  zperp = std::stod(v_line[7]);


  if (v_line[8] != "bending")
    throw std::invalid_argument("missing bending argument for polymer/no_tether.");
  kappa = std::stod(v_line[9]);


  if (v_line[10] != "temp")
    throw std::invalid_argument("missing temp argument for polymer/no_tether.");
  temp = std::stod(v_line[11]);

  if (v_line[12] != "seed")
    throw std::invalid_argument("missing seed argument for polymer/no_tether.");
  seed = std::stoi(v_line[13]);

  resize(Nbeads);
  costhetas.resize(Nbeads-2);
  gen.seed(seed);


  
}



/* ---------------------------------------------------------------------------- */
/* Compute the unit vectors tangent to the beads and rods, and then friction. */
/* Also computes cos(theta_i) relevant for potential force. */
/* Call only at the start of the time step (NOT the midstep). */
/* ---------------------------------------------------------------------------- */
void Polymer::compute_tangents_and_friction(const Eigen::Ref<const Eigen::Matrix3Xd> &xs)
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

double Polymer::get_bending() const
{
  return kappa;
}

int Polymer::get_seed() const
{
  return seed;
}

  
double Polymer::get_timescale(double dt) const
{
  return dt*temp/(bondlength*zperp);
}






  
};
