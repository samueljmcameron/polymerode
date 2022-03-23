#ifndef IOVTK_HPP
#define IOVTK_HPP

#include <string>
#include <vector>

#include "polymer.hpp"


namespace ioVTK {

  void writeVTKPolyData(std::string,const Polymer &);
  void readVTKPolyData(Polymer &, std::string);
  
  void writeVTKcollectionHeader(const std::string);
  void writeVTKcollectionMiddle(const std::string,
				const std::string, const double);
  void writeVTKcollectionFooter(const std::string);
  void restartVTKcollection(const std::string );
  

  
}



#endif
