#include "atom.hpp"

#include <vector>
#include <string>


namespace BeadRodPmer {
namespace Initialise {

void init_atoms(const std::vector<std::string> & splitvec,
		std::vector<Atom> &atoms_to_set, double springK,
		double dt, double tolerance);
}
}
