#include "atom.hpp"

#include <vector>
#include <string>


namespace BeadRodPmer {
namespace Initialise {

void init_atoms(const std::vector<std::string> & splitvec,
		std::vector<Atom> &atoms_to_set,
		Eigen::Vector3d &x0, Eigen::Vector3d &xN,
		double springK,double dt, double tolerance,int bufsteps=10000);
}
}
