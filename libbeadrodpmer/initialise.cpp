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
void init_atoms(const Pmer &pmer, Eigen::Ref<Eigen::Matrix3Xd> xs_to_set)
{

  NoTether noTethSlice = pmer.make_NoTether();
  noTethSlice.init_atoms_caret();
  
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



  dFdX_is.push_back(noTethSlice.initspringK*(noTethSlice.xs.col(0)-x0));
  
  dFdX_is.push_back(noTethSlice.initspringK*(noTethSlice.xs.col(b_index)-xN));



  noTethSlice.setup();
  
  noTethSlice.NoTether::single_step(t,noTethSlice.initdt,dFdX_is);
  t += noTethSlice.initdt;




  while ((noTethSlice.xs.col(b_index)-xN).norm() > noTethSlice.inittolerance
	 || (noTethSlice.xs.col(0)-x0).norm() > noTethSlice.inittolerance
	 || t < noTethSlice.equilibration_steps*noTethSlice.initdt) {


    dFdX_is[0] = noTethSlice.initspringK*(noTethSlice.xs.col(0)-x0);

    dFdX_is[1] = noTethSlice.initspringK*(noTethSlice.xs.col(b_index)-xN);


    noTethSlice.NoTether::single_step(t,noTethSlice.initdt,dFdX_is);
    t += noTethSlice.initdt;

  }

  for (int index = 0; index < noTethSlice.get_Nbeads(); index ++ ) 
    xs_to_set.col(index) = noTethSlice.xs.col(index);
  
  return ;
  
  
  
}
template void init_atoms(const NoTether &, Eigen::Ref<Eigen::Matrix3Xd>);
template void init_atoms(const SingleTether &, Eigen::Ref<Eigen::Matrix3Xd>);
template void init_atoms(const DoubleTether &, Eigen::Ref<Eigen::Matrix3Xd>);
}
}
