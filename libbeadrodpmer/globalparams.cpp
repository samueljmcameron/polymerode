#include <iostream>
#include "globalparams.hpp"


using namespace BeadRodPmer;
GlobalParams::GlobalParams(std::ifstream& input,
			   std::map<std::string,std::string> const& varMap,
			   std::string& lastline,
			   std::vector<std::string> EO_globals) :
  molecules(1),steps(1000000),timestep(0.001),
  dump_every(1000000),dump_file("dump.trj")
{

  std::string line;  
  if (input) {

    
    std::vector<std::string> splitvec;
    std::vector<std::string>::size_type expectedsize = 2;
    
    
    while(std::getline(input,line)) {
      
      std::size_t found = line.find_first_of("#");
      line = line.substr(0,found);
      
      if (line == "")  continue;
      
      
      splitvec = input::split_line(line);

      // check if keyword is global parameter, assume if not then
      // no longer defining global params

      int break_global = 0;

      // if build_solution, then global variable definitions are done
      for (auto &eo_glob : EO_globals)
	if (splitvec[0] == eo_glob) {
	  break_global = 1;
	  break;
	}

      if (break_global) break;

      
      if (splitvec.size() != expectedsize) {
	throw std::runtime_error("Error: invalid input file.");

      }


      input::convertVariable(splitvec[1],varMap);      
      if (splitvec[0] == "molecules") {
	input::isInt(splitvec[1],molecules,splitvec[0]);

      } else if (splitvec[0] == "steps") {
	input::isInt(splitvec[1],steps,splitvec[0]);
      } else if (splitvec[0] == "timestep") {
	input::isDouble(splitvec[1],timestep,splitvec[0]);
      } else if (splitvec[0] == "dump_every") {
	input::isInt(splitvec[1],dump_every,splitvec[0]);
      } else if (splitvec[0] == "dump_file") {
	dump_file = splitvec[1];
      }
    }
    
  }

  lastline = line;
  return; 
}
