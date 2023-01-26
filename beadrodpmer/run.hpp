#ifndef RUN_HPP
#define RUN_HPP

#include <string>
#include <Eigen/Core>

#include "polymer.hpp"

void run(BeadRodPmer::Polymer & ,Eigen::Ref<Eigen::Matrix3Xd> ,
	 int ,double , std::string , int );

#endif
