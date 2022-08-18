
#ifndef BEADRODPMER_ATOM_HPP
#define BEADRODPMER_ATOM_HPP

#include <Eigen/Core>


namespace BeadRodPmer {
class Atom {
public:
  Atom() {};
  ~Atom() {};

  
  Eigen::Vector3d R;
  Eigen::Vector3d Fpot, t_force;
  Eigen::Vector3d noise;
  Eigen::Vector3d unprojected_noise;
  Eigen::Vector3d tangent;
  Eigen::Matrix3d friction;



};
};
#endif