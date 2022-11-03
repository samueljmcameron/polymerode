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
		std::vector<Atom> &atoms_to_set,
		double springK,double dt, double tolerance,
		int equilibration_steps)
{

  
  BeadRodPmer::NoTether pmer(splitvec);

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
  
  std::vector<std::vector<double>> dFdX_is;

  double t = 0;



  dFdX_is.push_back(
		    {springK*(pmer.atoms[0].R(0)-x0(0)),
		     springK*(pmer.atoms[0].R(1)-x0(1)),
		     springK*(pmer.atoms[0].R(2)-x0(2))}
		    );
  
  dFdX_is.push_back(
		    {springK*(pmer.atoms[b_index].R(0)-xN(0)),
		     springK*(pmer.atoms[b_index].R(1)-xN(1)),
		     springK*(pmer.atoms[b_index].R(2)-xN(2))}
		    );




  pmer.compute_tangents_and_friction();
  
  pmer.set_Hhat();
  pmer.set_dCdlambda();
  pmer.set_G();
  pmer.single_step(t,dt,dFdX_is);
  t += dt;




  while ((pmer.atoms[b_index].R-xN).norm() > tolerance
	 || (pmer.atoms[0].R-x0).norm() > tolerance || t < equilibration_steps*dt) {


    dFdX_is[0][0] = springK*(pmer.atoms[0].R(0)-x0(0));
    dFdX_is[0][1] = springK*(pmer.atoms[0].R(1)-x0(1));
    dFdX_is[0][2] = springK*(pmer.atoms[0].R(2)-x0(2));

    dFdX_is[1][0] = springK*(pmer.atoms[b_index].R(0)-xN(0));
    dFdX_is[1][1] = springK*(pmer.atoms[b_index].R(1)-xN(1));
    dFdX_is[1][2] = springK*(pmer.atoms[b_index].R(2)-xN(2));

    pmer.single_step(t,dt,dFdX_is);
    t += dt;

  }

  for (int index = 0; index < pmer.get_Nbeads(); index ++ ) 
    atoms_to_set.at(index).R = pmer.atoms[index].R;
  
  return ;
  
  
  
}

}
}
