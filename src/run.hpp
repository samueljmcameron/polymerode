#ifndef RUN_HPP
#define RUN_HPP


#include "globalparams.hpp"
#include "double_tether.hpp"
#include "single_tether.hpp"
#include "no_tether.hpp"

void run(GlobalParams& gp, DoubleTether& pmer);
void run(GlobalParams& gp, SingleTether& pmer);
void run(GlobalParams& gp, NoTether& pmer);

#endif
