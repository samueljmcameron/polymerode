#ifndef BEADRODPMER_GLOBALPARAMS_HPP
#define BEADRODPMER_GLOBALPARAMS_HPP

#include <fstream>
#include <set>
#include <map>

#include "input.hpp"


namespace BeadRodPmer {
class GlobalParams {
private:

  std::set<std::string> pset {"molecules", "steps",
      "timestep", "dump_every", "dump_file"};
  
  int default_molecules;
  int default_steps;
  double default_timestep;
  int default_dump_every;
  std::string default_dump_file;
  
public:
  GlobalParams(std::ifstream&,
	       std::map<std::string,std::string> const &,
	       std::string&);

  int molecules;
  int steps;
  double timestep;
  int dump_every;
  std::string dump_file;


  
};
};
#endif
