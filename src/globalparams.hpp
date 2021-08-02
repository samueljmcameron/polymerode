#ifndef GLOBALPARAMS_HPP
#define GLOBALPARAMS_HPP

#include <fstream>
#include <set>

#include "input.hpp"

class GlobalParams {
private:

  std::set<std::string> pset {"molecules", "steps",
      "timestep", "dump_every", "dump_file"};
  
  const int default_molecules = 10;
  const int default_steps = 1000000;
  const double default_timestep = 0.001;
  const int default_dump_every = 1000000;
  const std::string default_dump_file{"dump.trj"};
  
public:
  GlobalParams(std::ifstream&,
	       std::map<std::string,std::string>,
	       std::string&);

  int molecules;
  int steps;
  double timestep;
  int dump_every;
  std::string dump_file;
  
};

#endif
