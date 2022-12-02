#ifndef BEADRODPMER_IOVTK_HPP
#define BEADRODPMER_IOVTK_HPP

#include <string>
#include <vector>

#include <Eigen/Core>

namespace BeadRodPmer {
namespace ioVTK {


  void writeVTKPolyData(std::string,const Eigen::Ref<Eigen::Matrix3Xd>,
			bool connect_atoms=true);

  void readVTKPolyData(Eigen::Ref<Eigen::Matrix3Xd>, std::string);
  
  void writeVTKcollectionHeader(const std::string);
  void writeVTKcollectionMiddle(const std::string,
				const std::string, const double);
  void writeVTKcollectionFooter(const std::string);
  void restartVTKcollection(const std::string );
  
}

};

#endif
