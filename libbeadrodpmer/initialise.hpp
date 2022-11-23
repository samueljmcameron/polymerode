#ifndef BEADRODPMER_INITIALISE_HPP
#define BEADRODPMER_INITIALISE_HPP

#include <vector>
#include <string>

namespace BeadRodPmer {
namespace Initialise {

template<typename Pmer>
void init_atoms(const Pmer &,Eigen::Ref<Eigen::Matrix3Xd> );
}
}

#endif
