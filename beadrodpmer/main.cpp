
#include "BeadRodPmerConfig.h"

#include "input.hpp"
#include "polymer.hpp"
#include "no_tether.hpp"
#include "single_tether.hpp"
#include "double_tether.hpp"
#include "initialise.hpp"
#include "iovtk.hpp"

#include "run.hpp"

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <map>



int main(int argc, char* argv[])
{
  std::cout << argv[0] << " version " << beadrodpmer_VERSION_MAJOR << "."
	    << beadrodpmer_VERSION_MINOR << std::endl;


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

    }/* else if (strcmp(argv[iarg],"-endtype") == 0) {
      endtype = argv[iarg+1];
      iarg += 2;
      } */
    else {
      std::cerr << "Error: invalid command line variable specification."
		<< std::endl;
      return EXIT_FAILURE;
    }
  }


  std::unique_ptr<BeadRodPmer::Polymer> pmer;
  
  
  if (not infile.is_open()) {
    std::cerr << "Error: need to specify input file." << std::endl;
    return EXIT_FAILURE;
  };

  
  std::string line;
  std::vector<std::string> v_line;
  
  std::string polymertype;
  while (std::getline(infile,line)) {
    
    if (line == "" || line[0] == '#') continue;
    
    input::convertVariable(line,variables);
    
    v_line = input::split_line(line);
    
    polymertype = v_line[0];
    v_line.erase(v_line.begin());
    
    if (polymertype == "double_tether") {
      pmer = std::make_unique<BeadRodPmer::DoubleTether>(v_line);
    } else if (polymertype == "single_tether") {
      pmer = std::make_unique<BeadRodPmer::SingleTether>(v_line);	
    } else if (polymertype == "no_tether") {
      pmer = std::make_unique<BeadRodPmer::NoTether>(v_line);
    } else {
      std::cerr << "Error: Incorrect polymer type." << std::endl;
      return EXIT_FAILURE;
    }
    
    break;
    
  }

  Eigen::Matrix3Xd xs(3,pmer->get_Nbeads());



  while (std::getline(infile,line)) {
    
    if (line == "" || line[0] == '#') continue;
    
    input::convertVariable(line,variables);
    
    v_line = input::split_line(line);

    if (v_line[0] == "initialise" ) {
      v_line.erase(v_line.begin());

      // initialise the xs vector (doing some unnecessary but illustrative slicing)
      BeadRodPmer::Initialise::init_atoms_relaxed_caret(v_line,*pmer,
							xs.middleCols(0,pmer->get_Nbeads()));

      
    } else if (v_line[0] == "read") {
      BeadRodPmer::ioVTK::readVTKPolyData(xs,v_line.at(1));

    } else  {
      std::cerr << " Invalid argument in input file" << std::endl;
      return EXIT_FAILURE;
    }

    break;
  }
  
  double dt;
  std::string dump_fname;
  int dump_every;

  int count = 0; // count that all arguments are included
  while (std::getline(infile,line)) {
    
    if (line == "" || line[0] == '#') continue;
    
    input::convertVariable(line,variables);
    
    v_line = input::split_line(line);

    if (v_line.size() != 2) {
      std::cerr << " Invalid argument in input file" << std::endl;
      return EXIT_FAILURE;
    }
  
    if (v_line[0] == "timestep") {
      dt = std::stod(v_line[1]);
      count += 1;
    } else if (v_line[0] == "dump_file") {
      dump_fname = v_line[1];
      count += 1;
    } else if (v_line[0] == "dump_every") {
      dump_every = std::stoi(v_line[1]);
      count += 1;
    } else {
      std::cerr << " Invalid argument in input file" << std::endl;
      return EXIT_FAILURE;
    }
    if (count == 3) break;
  }


  while (std::getline(infile,line)) {

    if (line == "" || line[0] == '#') continue;
    
    input::convertVariable(line,variables);
    
    v_line = input::split_line(line);


    if (v_line[0] == "run") {
      std::cout << "reading run" << std::endl;
      if (v_line.size() == 2) {
	int numsteps = std::stoi(v_line[1]);
	std::cout << "Running simulation of " << polymertype << " polymer." << std::endl;
	
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	
	run(*pmer,xs,numsteps,dt,dump_fname,dump_every);

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Run time = "
		  << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6
		  << "seconds." << std::endl;  
      } else if (v_line.size() == 4) {
	int numsteps = std::stoi(v_line[1]);
	int nbins = std::stoi(v_line[2]);
	int bin_every = std::stoi(v_line[3]);
	std::cout << "Running simulation of " << polymertype
		  << " polymer while also binning probability distribution of cosines."
		  << std::endl;
	
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	
	run(*pmer,xs,numsteps,dt,dump_fname,dump_every,true,nbins,bin_every);

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Run time = "
		  << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6
		  << "seconds." << std::endl;  

      } else {
	std::cerr << " Invalid argument in input file" << std::endl;
	return EXIT_FAILURE;
      }

	
    } else {
      std::cerr << " Invalid argument in input file" << std::endl;
    }
    break;

  }
    

  return 0;
  
}

