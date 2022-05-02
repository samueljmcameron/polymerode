
#include <iostream>

#include "utilfuncs.hpp"



void utilFuncs::calc_A_1(Eigen::Vector3d & A_1, const std::vector<Bond> & bonds,
			 const std::vector<Atom> & atoms,
			 const Eigen::VectorXd & tension)
{
  int offset = 2;
  A_1 =  atoms[0].friction*(utilFuncs::e_x*tension(offset-2)
			    +utilFuncs::e_y*tension(offset-1)
			    +utilFuncs::e_z*tension(offset)
			    -bonds[0].rod*tension(offset+1));
  return;

}


void utilFuncs::calc_A_i(Eigen::Vector3d & A_i, const std::vector<Bond> & bonds,
			 const std::vector<Atom> & atoms,
			 const Eigen::VectorXd & tension, int i)
{
  int offset = 2;


  A_i =  atoms[i].friction*(bonds[i-1].rod*tension(offset+i)
			    - bonds[i].rod*tension(offset+i+1));
  return;
}


void utilFuncs::calc_A_N(Eigen::Vector3d & A_N, const std::vector<Bond> & bonds,
			 const std::vector<Atom> & atoms,
			 const Eigen::VectorXd & tension)
{
  int offset = 2;
  int Nbeads = atoms.size();
  A_N =  atoms[Nbeads-1].friction*(utilFuncs::e_x*tension(offset+Nbeads)
				   +utilFuncs::e_y*tension(offset+Nbeads+1)
				   +utilFuncs::e_z*tension(offset+Nbeads+2)
				   +bonds[Nbeads-2].rod*tension(offset+Nbeads-1));
  return;
}

void utilFuncs::calc_S_mu(Eigen::Vector3d & S_mu,
			  const Eigen::Vector3d & B_mup1,
			  const Eigen::Vector3d & B_mu,
			  const Eigen::Vector3d & bond,
			  double bondlength)
{
  S_mu = (B_mup1 -B_mu - bond*(bond.dot(B_mup1-B_mu)))/bondlength;

  return;
}



double utilFuncs::compute_rhs_helper(const Eigen::Vector3d & A_mu,
				     const Eigen::Vector3d & A_mu_p1,
				     const Eigen::Vector3d & B_mu,
				     const Eigen::Vector3d & B_mu_p1,
				     const std::vector<Bond> & bonds,
				     double prefactor,double Delta_t,
				     double bondlength, int mu )

{
  double tmp;
  tmp = (A_mu_p1.dot(B_mu_p1-B_mu)
	 -A_mu_p1.dot(bonds[mu-1].rod)*(B_mu_p1-B_mu).dot(bonds[mu-1].rod)
	 )/bondlength*Delta_t*prefactor;
  

  tmp  -= (A_mu.dot(B_mu_p1-B_mu)
	   -A_mu.dot(bonds[mu-1].rod)*(B_mu_p1-B_mu).dot(bonds[mu-1].rod)
	   )/bondlength*Delta_t*prefactor;


  return tmp; 
}
