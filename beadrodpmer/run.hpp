#ifndef RUN_HPP
#define RUN_HPP

#include <string>
#include <Eigen/Core>

#include "polymer.hpp"

void run(BeadRodPmer::Polymer & ,Eigen::Ref<Eigen::Matrix3Xd> ,
	 int ,double , std::string , int ,bool histogram=false,int nbins = 0,
	 int bin_every = -1);

#endif
