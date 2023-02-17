#include "run.hpp"
#include "iovtk.hpp"
#include <fstream>
#include <iostream>

void bin_costhetas(Eigen::Ref<Eigen::MatrixXd> hist,BeadRodPmer::Polymer &pmer)
{

  double spacing = 2.0/hist.cols();

  int ibin;

  for (int mu = 0; mu < hist.rows(); mu++) {

    ibin = static_cast<int> ((pmer.bonds.col(mu+1).dot(pmer.bonds.col(mu))+1)/spacing);
    if (ibin == hist.cols()) ibin -= 1;

    hist(mu,ibin) += 1.0;

  }

    
}


void run(BeadRodPmer::Polymer & pmer,Eigen::Ref<Eigen::Matrix3Xd> xs,
	 int nsteps,double dt, std::string dump_name, int dump_every,
	 bool histogram,int nbins, int bin_every)
{

  Eigen::Matrix3Xd Fs(3,pmer.get_Nbeads());

  int Nbeads = pmer.get_Nbeads();

  std::vector<double> end2end_vec;

  
  int itermax = 50;
  double t = 0;

  // save at time t = 0 always.

  std::string collection_name = dump_name + std::string(".pvd");
  BeadRodPmer::ioVTK::writeVTKcollectionHeader(collection_name);

  
  std::string fname = dump_name + std::string("_") + std::to_string(0) + std::string(".vtp");
  BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs);
  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);

  

  pmer.setup(xs);


  Eigen::MatrixXd hist(pmer.get_Nbeads()-2,nbins);
  hist.setZero();

  
  int iterations;
  int ibin;
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

    if (histogram && i % bin_every == 0) {
      bin_costhetas(hist,pmer);
    }

    if (i % bin_every == 0) {
      end2end_vec.push_back((xs.col(Nbeads-1)-xs.col(0)).squaredNorm());
    }
    
  }

  if (histogram) {
    std::string histName = dump_name + std::string(".hist");
    std::ofstream histFile(histName);
    if (histFile.is_open()) {
      histFile << hist;
    } else {
      std::cerr << "could not open histogram file." << std::endl;
    }
  }

  std::string end2endName = dump_name + std::string(".end2end");
  std::ofstream end2endFile(end2endName);
  if (end2endFile.is_open()) {
    double end2end_variance = 0.0;
    double end2end_average = 0.0;
    for (auto val : end2end_vec)
      end2end_average += val;

    end2end_average /= end2end_vec.size();

    for (auto val : end2end_vec)
      end2end_variance += (val-end2end_average)*(val-end2end_average);

    end2end_variance /= (end2end_vec.size()-1);

    end2endFile << "# average variance number-of-samples" << std::endl
		<< end2end_average << " " << end2end_variance << " "
		<< end2end_vec.size() << std::endl;
  } else {
    std::cerr << "could not open end-to-end file." << std::endl;
  }
  
  
  BeadRodPmer::ioVTK::writeVTKcollectionFooter(collection_name);


  return;
}

