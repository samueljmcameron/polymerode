#ifndef BEADRODPMER_GLOBALPARAMS_HPP
#define BEADRODPMER_GLOBALPARAMS_HPP

#include <fstream>
#include <set>
#include <map>

#include "input.hpp"


namespace BeadRodPmer {
class GlobalParams {
public:
  GlobalParams(std::ifstream&,
	       std::map<std::string,std::string> const &,
	       std::string&,
	       std::vector<std::string> EO_globals = {"no_tether","single_tether","double_tether"});

  int molecules;
  int steps;
  double timestep;
  int dump_every;
  std::string dump_file;


  
};
};
#endif
