#ifndef BEADRODPMER_INITIALISE_HPP
#define BEADRODPMER_INITIALISE_HPP

#include <vector>
#include <string>


namespace BeadRodPmer {
namespace Initialise {

void init_atoms(const std::vector<std::string> & splitvec,
		Eigen::Ref<Eigen::Matrix3Xd> xs_to_set,
		double springK,double dt, double tolerance,
		int equilibration_steps);
}
}

#endif
