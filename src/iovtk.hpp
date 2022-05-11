#ifndef IOVTK_HPP
#define IOVTK_HPP

#include <string>
#include <vector>


namespace ioVTK {

  template <typename pMer>
  void writeVTKPolyData(std::string,const pMer &);
  template <typename pMer>
  void readVTKPolyData(pMer &, std::string);
  
  void writeVTKcollectionHeader(const std::string);
  void writeVTKcollectionMiddle(const std::string,
				const std::string, const double);
  void writeVTKcollectionFooter(const std::string);
  void restartVTKcollection(const std::string );
  

  
}



#endif
