#include <Eigen/Core>
#include <iostream>

#include "double_tether.hpp"

#include "initialise.hpp"


namespace BeadRodPmer {
namespace Initialise {
/* -------------------------------------------------------------------------- */
/* This function creates a polymer with endpoints at specific locations in a
   way that is not unreasonably asymmetric (as it would be if one end were
   tethered). 
/* -------------------------------------------------------------------------- */

template <typename Pmer>
void init_atoms(const Pmer &pmer, Eigen::Ref<Eigen::Matrix3Xd> xs)
{

  Eigen::Matrix3Xd Fs(3,pmer.get_Nbeads());

  
  NoTether noTethSlice = pmer.make_NoTether();
  noTethSlice.init_atoms_caret(xs);
  
  Eigen::Vector3d x0 = noTethSlice.get_x0();
  Eigen::Vector3d xN = noTethSlice.get_xN();

  
  // delete all nucleation sites
  noTethSlice.nuc_beads.clear();

  // then just add a single bead at the end of the polymer to guide it to
  // the desired location

  int b_index = noTethSlice.get_Nbeads()-1;


  noTethSlice.nuc_beads.push_back(0);
  noTethSlice.nuc_beads.push_back(b_index);



  // and add derivative of potential
  
  std::vector<Eigen::Vector3d> dFdX_is;

  double t = 0;



  dFdX_is.push_back(noTethSlice.initspringK*(xs.col(0)-x0));
  
  dFdX_is.push_back(noTethSlice.initspringK*(xs.col(b_index)-xN));



  noTethSlice.setup(xs);
  
  noTethSlice.NoTether::single_step(xs,Fs,t,noTethSlice.initdt,dFdX_is);
  t += noTethSlice.initdt;


  while ((xs.col(b_index)-xN).norm() > noTethSlice.inittolerance
	 || (xs.col(0)-x0).norm() > noTethSlice.inittolerance
	 || t < noTethSlice.equilibration_steps*noTethSlice.initdt) {


    dFdX_is[0] = noTethSlice.initspringK*(xs.col(0)-x0);

    dFdX_is[1] = noTethSlice.initspringK*(xs.col(b_index)-xN);


    noTethSlice.NoTether::single_step(xs,Fs,t,noTethSlice.initdt,dFdX_is);
    t += noTethSlice.initdt;


  }
  
  return ;
  
  
  
}
template void init_atoms(const NoTether &, Eigen::Ref<Eigen::Matrix3Xd>);
template void init_atoms(const SingleTether &, Eigen::Ref<Eigen::Matrix3Xd>);
template void init_atoms(const DoubleTether &, Eigen::Ref<Eigen::Matrix3Xd>);
}
}
