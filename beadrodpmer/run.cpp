#include "run.hpp"
#include "iovtk.hpp"
#include "initialise.hpp"
#include <chrono>
#include <iostream>
#include <functional>

class movingEnds
{
private:
  
  const double rad,omega;
  
  const Eigen::Vector3d x0,xN;
  
  Eigen::Vector3d dX0dt_val, X0_t_val;

  
public:
  movingEnds(double rad, double omega,
	     const Eigen::Vector3d &x0,
	     const Eigen::Vector3d &xN)
    : rad(rad),omega(omega),x0(x0),xN(xN) {};

  
  Eigen::Vector3d X0_t(double t) {
    X0_t_val(0) = x0(0) + rad*(cos(omega*t)-1);
    X0_t_val(1) = x0(1) + rad*sin(omega*t);
    X0_t_val(2) = x0(2);
    return X0_t_val;
  }

  Eigen::Vector3d dX0dt(double t)
  {
    dX0dt_val(0) = -rad*omega*sin(omega*t);
    dX0dt_val(1) = rad*omega*cos(omega*t);
    dX0dt_val(2) = 0;
    return dX0dt_val;
  }



  Eigen::Vector3d XN_t(double t) {
    return xN;
  }

  Eigen::Vector3d dXNdt(double t) {
    return {0,0,0};
  }
   
};




void run(BeadRodPmer::GlobalParams& gp, BeadRodPmer::DoubleTether& pmer)
{

  // making the arrays too big to demonstrate that there is minimal slow down
  // when slicing an array (necessary to show if I want multiple polymers
  // interacting)
  Eigen::Matrix3Xd xs(3,pmer.get_Nbeads()+20);
  Eigen::Matrix3Xd Fs(3,pmer.get_Nbeads()+20);

  int Nbeads = pmer.get_Nbeads();
  int offset = 5;

  
  BeadRodPmer::Initialise::init_atoms(pmer,xs.middleCols(offset,Nbeads));

  int numsteps = gp.steps;

  double dt = gp.timestep;

  int dump_every = gp.dump_every;

  double t = 0;

  movingEnds move(0.1,100,xs.middleCols(offset,Nbeads).col(0),xs.middleCols(offset,Nbeads).col(pmer.get_Nbeads()-1));

  auto X0_t = std::bind(&movingEnds::X0_t,&move,std::placeholders::_1);
  auto XN_t = std::bind(&movingEnds::XN_t,&move,std::placeholders::_1);
  auto dX0dt = std::bind(&movingEnds::dX0dt,&move,std::placeholders::_1);
  auto dXNdt = std::bind(&movingEnds::dXNdt,&move,std::placeholders::_1);

  std::vector<Eigen::Vector3d> dFdX_is;

  for (int index = 0; index < pmer.nuc_beads.size(); index ++ ) 
    dFdX_is.push_back({0,0,0});


  // save at time t = 0 always.

  std::string collection_name = gp.dump_file + std::string(".pvd");
  BeadRodPmer::ioVTK::writeVTKcollectionHeader(collection_name);

  
  std::string fname = gp.dump_file + std::string("_") + std::to_string(0) + std::string(".vtp");
  BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs.middleCols(offset,Nbeads));
  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  pmer.setup(xs.middleCols(offset,Nbeads));

  pmer.single_step(xs.middleCols(offset,Nbeads),Fs.middleCols(offset,Nbeads),t,dt,dFdX_is,X0_t,XN_t,dX0dt,dXNdt);
  t += dt;

  if (dump_every == 1) {

    fname = gp.dump_file + std::string("_") + std::to_string(1) + std::string(".vtp");
  
    BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs.middleCols(offset,Nbeads));

    BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  }

  // continue the remaining numstep - 2 steps
  for (int i = 2; i <= numsteps; i++) {


    pmer.single_step(xs.middleCols(offset,Nbeads),Fs.middleCols(offset,Nbeads),t,dt,dFdX_is,X0_t,XN_t,dX0dt,dXNdt);
    t += dt;
    if (i % dump_every == 0) {
      fname = gp.dump_file + std::string("_") + std::to_string(i) + std::string(".vtp");
  
      BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs.middleCols(offset,Nbeads));

      BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);
      
    }
      
  }


  BeadRodPmer::ioVTK::writeVTKcollectionFooter(collection_name);


  return;
}



void run(BeadRodPmer::GlobalParams& gp, BeadRodPmer::SingleTether& pmer)
{

  // making the arrays too big to demonstrate that there is minimal slow down
  // when slicing an array (necessary to show if I want multiple polymers
  // interacting)
  Eigen::Matrix3Xd xs(3,pmer.get_Nbeads()+20);
  Eigen::Matrix3Xd Fs(3,pmer.get_Nbeads()+20);

  int Nbeads = pmer.get_Nbeads();
  int offset = 5;

  
  BeadRodPmer::Initialise::init_atoms(pmer,xs.middleCols(offset,Nbeads));
  int numsteps = gp.steps;

  double dt = gp.timestep;

  int dump_every = gp.dump_every;

  double t = 0;

  // since only fixing one end, the XN argument of the move class is ignored
  // so putting in the zero vector {0,0,0} as a placeholder
  movingEnds move(0.1,100,xs.middleCols(offset,Nbeads).col(0),{0,0,0});


  auto X0_t = std::bind(&movingEnds::X0_t,&move,std::placeholders::_1);
  auto dX0dt = std::bind(&movingEnds::dX0dt,&move,std::placeholders::_1);

  

  std::vector<Eigen::Vector3d> dFdX_is;

  for (int index = 0; index < pmer.nuc_beads.size(); index ++ ) 
    dFdX_is.push_back({0,0,0});
  
  // save at time t = 0 always.

  std::string collection_name = gp.dump_file + std::string(".pvd");
  BeadRodPmer::ioVTK::writeVTKcollectionHeader(collection_name);

  
  std::string fname = gp.dump_file + std::string("_") + std::to_string(0) + std::string(".vtp");
  BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs.middleCols(offset,Nbeads));
  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  pmer.setup(xs.middleCols(offset,Nbeads));
  pmer.single_step(xs.middleCols(offset,Nbeads),Fs.middleCols(offset,Nbeads),t,dt,dFdX_is,X0_t,dX0dt);
  t += dt;

  if (dump_every == 1) {

    fname = gp.dump_file + std::string("_") + std::to_string(1) + std::string(".vtp");
  
    BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs.middleCols(offset,Nbeads));

    BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  }

  // continue the remaining numstep - 2 steps
  for (int i = 2; i <= numsteps; i++) {


    pmer.single_step(xs.middleCols(offset,Nbeads),Fs.middleCols(offset,Nbeads),t,dt,dFdX_is,X0_t,dX0dt);
    t += dt;

    if (i % dump_every == 0) {
      fname = gp.dump_file + std::string("_") + std::to_string(i) + std::string(".vtp");
  
      BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs.middleCols(offset,Nbeads));

      BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);
      
    }
      
  }


  BeadRodPmer::ioVTK::writeVTKcollectionFooter(collection_name);



  return;
}

void run(BeadRodPmer::GlobalParams& gp, BeadRodPmer::NoTether& pmer)
{

  // making the arrays too big to demonstrate that there is minimal slow down
  // when slicing an array (necessary to show if I want multiple polymers
  // interacting)
  Eigen::Matrix3Xd xs(3,pmer.get_Nbeads()+20);
  Eigen::Matrix3Xd Fs(3,pmer.get_Nbeads()+20);

  int Nbeads = pmer.get_Nbeads();
  int offset = 5;
  
  BeadRodPmer::Initialise::init_atoms(pmer,xs.middleCols(offset,Nbeads));
  int numsteps = gp.steps;

  double dt = gp.timestep;

  int dump_every = gp.dump_every;

  double t = 0;


  std::vector<Eigen::Vector3d> dFdX_is;

  for (int index = 0; index < pmer.nuc_beads.size(); index ++ ) 
    dFdX_is.push_back({0,0,0});
  
  // save at time t = 0 always.

  std::string collection_name = gp.dump_file + std::string(".pvd");
  BeadRodPmer::ioVTK::writeVTKcollectionHeader(collection_name);

  
  std::string fname = gp.dump_file + std::string("_") + std::to_string(0) + std::string(".vtp");
  BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs.middleCols(offset,Nbeads));
  BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  pmer.setup(xs.middleCols(offset,Nbeads));
  pmer.single_step(xs.middleCols(offset,Nbeads),Fs.middleCols(offset,Nbeads),t,dt,dFdX_is);
  t += dt;

  if (dump_every == 1) {

    fname = gp.dump_file + std::string("_") + std::to_string(1) + std::string(".vtp");
  
    BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs.middleCols(offset,Nbeads));

    BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);


  }

  // continue the remaining numstep - 2 steps
  for (int i = 2; i <= numsteps; i++) {

    
    pmer.single_step(xs.middleCols(offset,Nbeads),Fs.middleCols(offset,Nbeads),t,dt,dFdX_is);
    t += dt;

    if (i % dump_every == 0) {
      fname = gp.dump_file + std::string("_") + std::to_string(i) + std::string(".vtp");
  
      BeadRodPmer::ioVTK::writeVTKPolyData(fname,xs.middleCols(offset,Nbeads));

      BeadRodPmer::ioVTK::writeVTKcollectionMiddle(collection_name,fname,t);
      
    }
      
  }


  BeadRodPmer::ioVTK::writeVTKcollectionFooter(collection_name);


  return;
}
