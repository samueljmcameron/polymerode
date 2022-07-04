#ifndef RUN_HPP
#define RUN_HPP


#include "globalparams.hpp"
#include "double_tether.hpp"
#include "single_tether.hpp"
#include "no_tether.hpp"

void run(BeadRodPmer::GlobalParams& gp, BeadRodPmer::DoubleTether& pmer);
void run(BeadRodPmer::GlobalParams& gp, BeadRodPmer::SingleTether& pmer);
void run(BeadRodPmer::GlobalParams& gp, BeadRodPmer::NoTether& pmer);

#endif
