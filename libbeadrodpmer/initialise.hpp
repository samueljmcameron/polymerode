#ifndef BEADRODPMER_INITIALISE_HPP
#define BEADRODPMER_INITIALISE_HPP

#include <vector>
#include <string>
#include <Eigen/Core>


#include "polymer.hpp"

namespace BeadRodPmer {
namespace Initialise {
  
void init_atoms_relaxed_caret(const std::vector<std::string> &,const Polymer &,
			      Eigen::Ref<Eigen::Matrix3Xd>);

void init_atoms_rand(Eigen::Ref<Eigen::Matrix3Xd> ,Eigen::Vector3d ,Eigen::Vector3d ,
		     double , int );

void init_atoms_line(Eigen::Ref<Eigen::Matrix3Xd>,Eigen::Vector3d ,Eigen::Vector3d ,
		     double );
  
void init_atoms_caret(Eigen::Ref<Eigen::Matrix3Xd> , Eigen::Vector3d ,Eigen::Vector3d ,
		      double,int);

void init_atoms_equilibrium(Eigen::Ref<Eigen::Matrix3Xd> ,
			    Eigen::Vector3d ,Eigen::Vector3d ,
			    double, double ,int );

  
}
}

#endif
