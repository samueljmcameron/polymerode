
#ifndef BEADRODPMER_ATOM_HPP
#define BEADRODPMER_ATOM_HPP

#include <Eigen/Core>


namespace BeadRodPmer {
class Atom {
public:


  void resize(int Nbeads) {
    xs.resize(Eigen::NoChange,Nbeads);
    Fpots.resize(Eigen::NoChange,Nbeads);
    t_forces.resize(Eigen::NoChange,Nbeads);
    noises.resize(Eigen::NoChange,Nbeads);
    unprojected_noises.resize(Eigen::NoChange,Nbeads);
    tangents.resize(Eigen::NoChange,Nbeads);
    bonds.resize(Eigen::NoChange,Nbeads-1);
    frictions.resize(Nbeads);
  }
  
  Eigen::Matrix<double, 3, Eigen::Dynamic> xs;


  Eigen::Matrix<double, 3, Eigen::Dynamic> Fpots,t_forces;
  Eigen::Matrix<double, 3, Eigen::Dynamic> noises,unprojected_noises;
  Eigen::Matrix<double, 3, Eigen::Dynamic> tangents;

  Eigen::Matrix<double, 3, Eigen::Dynamic> bonds;
  
  std::vector<Eigen::Matrix3d> frictions;



};
};
#endif
