#ifndef BEADRODPMER_INITIALISE_HPP
#define BEADRODPMER_INITIALISE_HPP

#include "atom.hpp"

#include <vector>
#include <string>


namespace BeadRodPmer {
namespace Initialise {

void init_atoms(const std::vector<std::string> & splitvec,
		Atom &atoms_to_set,
		double springK,double dt, double tolerance,
		int equilibration_steps);
}
}

#endif
