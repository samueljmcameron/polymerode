#ifndef IOVTK_HPP
#define IOVTK_HPP

#include <string>
#include <vector>
#include "atom.hpp"
#include "bond.hpp"

#include <Eigen/Core>


namespace utilFuncs {


  const Eigen::Vector3d e_x = (Eigen::Vector3d() << 1, 0 , 0).finished();
  const Eigen::Vector3d e_y = (Eigen::Vector3d() << 0, 1 , 0).finished();
  const Eigen::Vector3d e_z = (Eigen::Vector3d() << 0, 0 , 1).finished();
  

  void calc_A_1(Eigen::Vector3d & , const std::vector<Bond> & ,
		const std::vector<Atom> & ,
		const Eigen::VectorXd & );

  void calc_A_i(Eigen::Vector3d & , const std::vector<Bond> & ,
		const std::vector<Atom> & ,
		const Eigen::VectorXd & , int );

  
  void calc_A_N(Eigen::Vector3d & , const std::vector<Bond> &,
		const std::vector<Atom> & ,
		const Eigen::VectorXd & );

  void calc_S_mu(Eigen::Vector3d & , const Eigen::Vector3d & ,
		 const Eigen::Vector3d & ,const Eigen::Vector3d & ,
		 double);

  double compute_rhs_helper(const Eigen::Vector3d &, const Eigen::Vector3d &,
			    const Eigen::Vector3d &, const Eigen::Vector3d &,
			    const std::vector<Bond> & ,
			    double ,double,double, int);

  
  
}



#endif
