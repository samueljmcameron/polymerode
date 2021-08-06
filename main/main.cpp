
#include "PolyConfig.h"
#include "polymer.hpp"
#include "globalparams.hpp"
#include "input.hpp"
#include "run.hpp"
#include "runtests.hpp"

#include <Eigen/SparseLU>

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <map>



int main(int argc, char* argv[])
{
  std::cout << argv[0] << " version " << poly_VERSION_MAJOR << "."
	    << poly_VERSION_MINOR << std::endl;


  std::ifstream infile;

  std::map<std::string,std::string> variables;

  std::string simulation_type = "run";

  int iarg = 1;  
  while(iarg < argc) {
    if (strcmp(argv[iarg],"-in") == 0) {
      if (iarg+1 == argc) {
	std::cerr << "Error: input flag '-in' specified, but no file given."
		  << std::endl;
	return EXIT_FAILURE;
      }
      infile.open(argv[iarg+1]);
      iarg += 2;
      
    } else if (strcmp(argv[iarg],"-var") == 0) {
      
      if (iarg + 2 >= argc) {
	std::cerr << "Error: invalid command line variable specification."
		  << std::endl;
	return EXIT_FAILURE;
      }
      variables[argv[iarg+1]] = argv[iarg+2];
      iarg += 3;
    } else if (strcmp(argv[iarg],"-testnoise") == 0) {
      simulation_type = "test";
      iarg += 1;
    } else {
      std::cerr << "Error: invalid command line variable specification."
		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  
  if (not infile.is_open()) {
    std::cerr << "Error: need to specify input file." << std::endl;
    return EXIT_FAILURE;
  }

  std::string line;
  std::vector<std::string> splitvec;
    
  GlobalParams gp(infile,variables,line);

  line = line.substr(0,line.find_first_of("#"));
  splitvec = input::split_line(line);
  

  while (std::getline(infile,line) &&
	 splitvec[0] != "build_polymer") {

    line = line.substr(0,line.find_first_of("#"));
    splitvec = input::split_line(line);
    
  }

  splitvec.erase(splitvec.begin());

  for (auto &c : splitvec)
    input::convertVariable(c,variables);
  
  
  Polymer pmer(splitvec);
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();




  
  if (simulation_type == "run") {
    std::cout << "Running simulation of polymer." << std::endl;
    run(gp,pmer);
  } else if (simulation_type == "test") {
    std::cout << "Testing simulation of polymer." << std::endl;
    test_noise(gp,pmer);
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Run time = "
	    << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6
	    << "seconds." << std::endl;  
  return 0;
  
}

