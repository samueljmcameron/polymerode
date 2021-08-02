
#ifndef ATOM_HPP
#define ATOM_HPP

#include <Eigen/Core>

class Atom {
public:
  Atom() {};
  ~Atom() {};

  
  Eigen::Vector3d R;
  Eigen::Vector3d F;
  Eigen::Vector3d noise;
  Eigen::Vector3d unprojected_noise;
  Eigen::Vector3d tangent;
  Eigen::Matrix3d friction;



};

#endif
