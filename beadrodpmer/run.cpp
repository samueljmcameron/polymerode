#include "run.hpp"
#include "iovtk.hpp"

void run(BeadRodPmer::Polymer & pmer,Eigen::Ref<Eigen::Matrix3Xd> xs,
	 int nsteps,double dt, std::string dump_name, int dump_every)
{

  Eigen::Matrix3Xd Fs(3,pmer.get_Nbeads());

  int Nbeads = pmer.get_Nbeads();
  

  int itermax = 50;
  double t = 0;

  // save at time t = 0 always.

  std::string collection_name = dump_name + std::string(".pvd");
  BeadRodPmer::ioVTK::writeVTKcollectionHeader(collection_name);

  
  std::string fname = dump_name + std::string("_") + std::to_string(0) + std::string(".vtp");
  BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs);
  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  pmer.setup(xs);

  int iterations;

  // time-step
  for (int i = 1; i <= nsteps; i++) {

    Fs.setZero();
    pmer.first_step(xs,Fs,dt);

    Fs.setZero();
    iterations = pmer.second_step(xs,Fs,dt,itermax);
    t += dt;

    
    if (iterations >= itermax) {
      
      t -= dt;
      
      fname = dump_name + std::string("_") + std::to_string(i-1) + std::string(".vtp");
  
      BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs);

      BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);
      BeadRodPmer::ioVTK::writeVTKcollectionFooter(collection_name);
      
      throw std::runtime_error("Polymer simulation failed to converge at step "
			       + std::to_string(i) + std::string("."));

    }

    

    if (i % dump_every == 0) {
      fname = dump_name + std::string("_") + std::to_string(i) + std::string(".vtp");
  
      BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs);

      BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);
      
    }
  }


  BeadRodPmer::ioVTK::writeVTKcollectionFooter(collection_name);


  return;
}

