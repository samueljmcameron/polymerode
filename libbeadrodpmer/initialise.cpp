#include "no_tether.hpp"
#include <Eigen/Core>
#include <iostream>

#include "initialise.hpp"

namespace BeadRodPmer {
namespace Initialise {
/* -------------------------------------------------------------------------- */
/* This function creates a polymer with endpoints at specific locations in a
   way that is not unreasonably asymmetric (as it would be if one end were
   tethered). */
/* -------------------------------------------------------------------------- */

void init_atoms(const std::vector<std::string> & splitvec,
		Eigen::Ref<Eigen::Matrix3Xd> xs_to_set,
		double springK,double dt, double tolerance,
		int equilibration_steps)
{

  
  BeadRodPmer::NoTether pmer(splitvec);

  pmer.init_atoms_caret();
  
  Eigen::Vector3d x0 = pmer.get_x0();
  Eigen::Vector3d xN = pmer.get_xN();

  
  // delete all nucleation sites
  pmer.nuc_beads.clear();

  // then just add a single bead at the end of the polymer to guide it to
  // the desired location

  int b_index = pmer.get_Nbeads()-1;


  pmer.nuc_beads.push_back(0);
  pmer.nuc_beads.push_back(b_index);



  // and add derivative of potential
  
  std::vector<Eigen::Vector3d> dFdX_is;

  double t = 0;



  dFdX_is.push_back(springK*(pmer.xs.col(0)-x0));
  
  dFdX_is.push_back(springK*(pmer.xs.col(b_index)-xN));



  pmer.setup();
  
  pmer.single_step(t,dt,dFdX_is);
  t += dt;




  while ((pmer.xs.col(b_index)-xN).norm() > tolerance
	 || (pmer.xs.col(0)-x0).norm() > tolerance || t < equilibration_steps*dt) {


    dFdX_is[0] = springK*(pmer.xs.col(0)-x0);

    dFdX_is[1] = springK*(pmer.xs.col(b_index)-xN);


    pmer.single_step(t,dt,dFdX_is);
    t += dt;

  }

  for (int index = 0; index < pmer.get_Nbeads(); index ++ ) 
    xs_to_set.col(index) = pmer.xs.col(index);
  
  return ;
  
  
  
}

}
}
