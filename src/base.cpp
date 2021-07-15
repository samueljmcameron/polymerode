#include "base.hpp"

#include <iostream>
#include <cmath>
#include <stdexcept>


/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
Base::Base(int Nin, double rodin, double Delta_tin, double zparain,
	   double zperpin, double tempin,double kappain,
	   Eigen::Vector3d x0in,int seedin)
  : Nbeads(Nin), rodlength(rodin), Delta_t(Delta_tin), zpara(zparain),
    zperp(zperpin), temp(tempin), kappa(kappain), seed(seedin)
{


  x0 = x0in;
  R_x.resize(Nbeads);
  R_y.resize(Nbeads);
  R_z.resize(Nbeads);

  Rtmp_x.resize(Nbeads);
  Rtmp_y.resize(Nbeads);
  Rtmp_z.resize(Nbeads);

  F_x.resize(Nbeads);
  F_y.resize(Nbeads);
  F_z.resize(Nbeads);
    
  u_x.resize(Nbeads-1);
  u_y.resize(Nbeads-1);
  u_z.resize(Nbeads-1);

  tu_x.resize(Nbeads);
  tu_y.resize(Nbeads);
  tu_z.resize(Nbeads);

  eta_x.resize(Nbeads);
  eta_y.resize(Nbeads);
  eta_z.resize(Nbeads);

  unproj_eta_x.resize(Nbeads);
  unproj_eta_y.resize(Nbeads);
  unproj_eta_z.resize(Nbeads);

  
  rhs_of_G.resize(Nbeads-1);
  dummy_for_noise.resize(Nbeads-1);
  Gmunu.resize(Nbeads-1,Nbeads-1);
  Hhat.resize(Nbeads-1,Nbeads-1);
  rhs_of_Hhat.resize(Nbeads-1);
  
  tension.resize(Nbeads-1);

  zetainvs.resize(Nbeads);  
  
  init_R();




  gen = std::mt19937(seed);
  dist = std::uniform_real_distribution<double>(-0.5,0.5) ;


  
}

/* -------------------------------------------------------------------------- */
/* Destructor */
/* -------------------------------------------------------------------------- */
Base::~Base()
{
}


/* -------------------------------------------------------------------------- */
/* Setting the initial configuration as a straight rod along x axis.
/* -------------------------------------------------------------------------- */
void Base::init_R()
{

  R_x(0) = 0;
  for (int i = 1; i< Nbeads; i++) {
    R_x(i) = R_x(i-1) + rodlength;
  }
  R_y.setZero();
  R_z.setZero();



  R_x.array() -= R_x(Nbeads-1) - x0(0);
  R_y.array() -= R_y(Nbeads-1) - x0(1);
  R_z.array() -= R_z(Nbeads-1) - x0(2);
  
  return;

}

/* ---------------------------------------------------------------------------- */
/* Compute the unit vectors tangent to the bonds (u) and beads (tu). */
/* ---------------------------------------------------------------------------- */
void Base::compute_tangents()
{

  double norm;

  norm = sqrt((R_x(1)-R_x(0))*(R_x(1)-R_x(0))
  	      +(R_y(1)-R_y(0))*(R_y(1)-R_y(0))
  	      +(R_z(1)-R_z(0))*(R_z(1)-R_z(0)));


  u_x(0) = (R_x(1)-R_x(0))/norm;
  u_y(0) = (R_y(1)-R_y(0))/norm;
  u_z(0) = (R_z(1)-R_z(0))/norm;

  tu_x(0) = u_x(0);
  tu_y(0) = u_y(0);
  tu_z(0) = u_z(0);


  for (int i = 1; i < Nbeads-1; i++) {

    norm = sqrt((R_x(i+1)-R_x(i))*(R_x(i+1)-R_x(i))
		+(R_y(i+1)-R_y(i))*(R_y(i+1)-R_y(i))
		+(R_z(i+1)-R_z(i))*(R_z(i+1)-R_z(i)));

    u_x(i) = (R_x(i+1)-R_x(i))/norm;
    u_y(i) = (R_y(i+1)-R_y(i))/norm;
    u_z(i) = (R_z(i+1)-R_z(i))/norm;

    norm = sqrt((u_x(i)+u_x(i-1))*(u_x(i)+u_x(i-1))
		+ (u_y(i)+u_y(i-1))*(u_y(i)+u_y(i-1))
		+ (u_z(i)+u_z(i-1))*(u_z(i)+u_z(i-1)));
    
    tu_x(i) = (u_x(i)+u_x(i-1))/norm;
    tu_y(i) = (u_y(i)+u_y(i-1))/norm;
    tu_z(i) = (u_z(i)+u_z(i-1))/norm;

    
  }

  tu_x(Nbeads-1) = u_x(Nbeads-2);
  tu_y(Nbeads-1) = u_y(Nbeads-2);
  tu_z(Nbeads-1) = u_z(Nbeads-2);
  
}

/* ---------------------------------------------------------------------------- */
/* Generate unprojected noises with correct variance. */
/* ---------------------------------------------------------------------------- */
void Base::set_unprojected_noise()
{

  double dum = sqrt(24*temp/Delta_t);
  double para = sqrt(zpara);
  double perp = sqrt(zperp);

  double rx,ry,rz;
  double dotting;
  
  for (int i = 0; i < Nbeads; i++) {
    rx = dist(gen);
    ry = dist(gen);
    rz = dist(gen);

    dotting = tu_x(i)*rx + tu_y(i)*ry + tu_z(i)*rz;
    
    unproj_eta_x(i) = dum*(perp*rx + (para-perp)*tu_x(i)*dotting);
    unproj_eta_y(i) = dum*(perp*ry + (para-perp)*tu_y(i)*dotting);
    unproj_eta_z(i) = dum*(perp*rz + (para-perp)*tu_z(i)*dotting);

  }

  return;
}


/* ---------------------------------------------------------------------------- */
/* RHS of G*eta = P. */
/* ---------------------------------------------------------------------------- */
void Base::set_rhs_of_G()
{

  for (int i = 0; i < Nbeads-1; i++) {
    rhs_of_G(i) = (unproj_eta_x(i+1)-unproj_eta_x(i))*u_x(i)
      +(unproj_eta_y(i+1)-unproj_eta_y(i))*u_y(i)
      +(unproj_eta_z(i+1)-unproj_eta_z(i))*u_z(i);
  }
  return;
  
}


/* -------------------------------------------------------------------------- */
/* Initialse G in G*eta = P (only call once). */
/* -------------------------------------------------------------------------- */
void Base::set_G()
{

  std::vector<T> coefficients = init_G_coeffsmatrix();

  Gmunu.setFromTriplets(coefficients.begin(),coefficients.end());

  solver.analyzePattern(Gmunu);

  return;

}


/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix G. */
/* -------------------------------------------------------------------------- */
std::vector<T> Base::init_G_coeffsmatrix()
{
  std::vector<T> coeffs;
  coeffs.push_back(T(0,0,2));
  for (int i = 1; i < Nbeads-1; i++) {
    coeffs.push_back(T(i,i,2));
    coeffs.push_back(T(i,i-1,-u_x(i)*u_x(i-1)-u_y(i)*u_y(i-1)
		       -u_z(i)*u_z(i-1)));

  }

  return coeffs;

  
}

/* -------------------------------------------------------------------------- */
/* Update matrix G. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void Base::update_G()
{
  // update first col by hand

  Gmunu.coeffRef(1,0) = (-u_x(1)*u_x(0)-u_y(1)*u_y(0)
			 -u_z(1)*u_z(0));


  // update middle columns
  for (int k = 1; k < Nbeads-2; k++) {


    int count = 0;
    for (SpMat::InnerIterator it(Gmunu,k); it; ++it) {
      if (count == 1) {
	it.valueRef() = (-u_x(k+1)*u_x(k)-u_y(k+1)*u_y(k)
			 -u_z(k+1)*u_z(k));

	

      }
      count += 1;
    }
  }
  return;
}


/* ---------------------------------------------------------------------------- */
/* RHS of Hhat*lambda = Q. */
/* ---------------------------------------------------------------------------- */
void Base::set_rhs_of_Hhat()
{

  double tmp1,tmp2,tmp3;

  for (int i = 0; i< Nbeads-1; i++) {

    zetainv_i = zetainvs[i+1];
    tmp1 = zetainv_i(0,0)*F_x(i+1) + zetainv_i(0,1)*F_y(i+1) + zetainv_i(0,2)*F_z(i+1);
    tmp2 = zetainv_i(1,0)*F_x(i+1) + zetainv_i(1,1)*F_y(i+1) + zetainv_i(1,2)*F_z(i+1);
    tmp3 = zetainv_i(2,0)*F_x(i+1) + zetainv_i(2,1)*F_y(i+1) + zetainv_i(2,2)*F_z(i+1);

    zetainv_i = zetainvs[i];
    tmp1 -= zetainv_i(0,0)*F_x(i) + zetainv_i(0,1)*F_y(i) + zetainv_i(0,2)*F_z(i);
    tmp2 -= zetainv_i(1,0)*F_x(i) + zetainv_i(1,1)*F_y(i) + zetainv_i(1,2)*F_z(i);
    tmp3 -= zetainv_i(2,0)*F_x(i) + zetainv_i(2,1)*F_y(i) + zetainv_i(2,2)*F_z(i);
    
    rhs_of_Hhat(i) = u_x(i)*tmp1+u_y(i)*tmp2 + u_z(i)*tmp3;

  }
  return;
  
}

void Base::set_zetainv()
{
  for (int i = 0; i < Nbeads; i++) {
    set_zetainv_i(i);

    zetainvs[i] = zetainv_i;
  }
  return;

}
void Base::set_zetainv_i(int i)
{

  double tmp = 1./zpara-1./zperp;
  
  zetainv_i(0,0) = 1./zperp + tmp*tu_x(i)*tu_x(i);
  zetainv_i(0,1) = tmp*tu_x(i)*tu_y(i);
  zetainv_i(0,2) = tmp*tu_x(i)*tu_z(i);

  zetainv_i(1,0) = tmp*tu_y(i)*tu_x(i);
  zetainv_i(1,1) = 1./zperp + tmp*tu_y(i)*tu_y(i);
  zetainv_i(1,2) = tmp*tu_y(i)*tu_z(i);

  zetainv_i(2,0) = tmp*tu_z(i)*tu_x(i);
  zetainv_i(2,1) = tmp*tu_z(i)*tu_y(i);
  zetainv_i(2,2) = 1./zperp + tmp*tu_z(i)*tu_z(i);

  return;
  

  
}


/* -------------------------------------------------------------------------- */
/* Compute forces on the particles. */
/* -------------------------------------------------------------------------- */
void Base::compute_uc_forces()
{
  F_x(0) = Force_0(R_x)+eta_x(0);
  
  F_y(0) = Force_0(R_y)+eta_y(0);
  F_z(0) = Force_0(R_z)+eta_z(0);
  
  F_x(1) = Force_1(R_x)+eta_x(1);
  F_y(1) = Force_1(R_y)+eta_y(1);
  F_z(1) = Force_1(R_z)+eta_z(1);
  
  for (int i = 2; i < Nbeads-2; i++) {
    F_x(i) = Force_mid(R_x,i)+eta_x(i);
    F_y(i) = Force_mid(R_y,i)+eta_y(i);
    F_z(i) = Force_mid(R_z,i)+eta_z(i);
    
  }
  
  
  F_x(Nbeads-2) = Force_Nm2(R_x)+eta_x(Nbeads-2);
  F_y(Nbeads-2) = Force_Nm2(R_y)+eta_y(Nbeads-2);
  F_z(Nbeads-2) = Force_Nm2(R_z)+eta_z(Nbeads-2);
  
  F_x(Nbeads-1) = Force_Nm1(R_x)+eta_x(Nbeads-1);
  F_y(Nbeads-1) = Force_Nm1(R_y)+eta_y(Nbeads-1);
  F_z(Nbeads-1) = Force_Nm1(R_z)+eta_z(Nbeads-1);
  return;
  
}

/* -------------------------------------------------------------------------- */
/* Helper functions for force calculation. */

double Base::Force_mid(Eigen::VectorXd &u,int i)
{

  double tmp1 = u_x(i)*u_x(i-1)+u_y(i)*u_y(i-1)+u_z(i)*u_z(i-1);
  double tmp2 = u_x(i+1)*u_x(i)+u_y(i+1)*u_y(i)+u_z(i+1)*u_z(i);
  double tmp3 = u_x(i-2)*u_x(i-1)+u_y(i-2)*u_y(i-1)+u_z(i-2)*u_z(i-1);

  return kappa/(rodlength*rodlength)*(u(i)-tmp1*u(i-1)-u(i+1)+tmp2*u(i)
				      +u(i-2)-tmp3*u(i-1)-u(i-1)+tmp1*u(i));
  
  /*return kappa/rodlength*(R(i-1) - R(i-2) - R_x(i) + R_x(i-1)
    + R_x(i+1) - R_x(i) - R_x(i+2) + R_x(i+1));*/
}


double Base::Force_0(Eigen::VectorXd &u)
{

  double tmp2 = u_x(0+1)*u_x(0)+u_y(0+1)*u_y(0)+u_z(0+1)*u_z(0);

  return kappa/(rodlength*rodlength)*(-u(0+1)+tmp2*u(0));
  
  /*  return kappa/rodlength*(- R(2) + R(1)); */
}

double Base::Force_1(Eigen::VectorXd &u)
{
  double tmp1 = u_x(1)*u_x(1-1)+u_y(1)*u_y(1-1)+u_z(1)*u_z(1-1);
  double tmp2 = u_x(1+1)*u_x(1)+u_y(1+1)*u_y(1)+u_z(1+1)*u_z(1);


  return kappa/(rodlength*rodlength)*(u(1)-tmp1*u(1-1)-u(1+1)+tmp2*u(1)
				      -u(1-1)+tmp1*u(1));
      
  /*return kappa/rodlength*(-R(1) + R(0) + R(2) - R(1) - R(3) + R(2)); */
}

double Base::Force_Nm2(Eigen::VectorXd &u)
{

  double tmp1 = u_x(Nbeads-2)*u_x(Nbeads-2-1)+u_y(Nbeads-2)*u_y(Nbeads-2-1)+u_z(Nbeads-2)*u_z(Nbeads-2-1);
  double tmp3 = u_x(Nbeads-2-2)*u_x(Nbeads-2-1)+u_y(Nbeads-2-2)*u_y(Nbeads-2-1)+u_z(Nbeads-2-2)*u_z(Nbeads-2-1);

  return kappa/(rodlength*rodlength)*(u(Nbeads-2)-tmp1*u(Nbeads-2-1)
				      +u(Nbeads-2-2)-tmp3*u(Nbeads-2-1)-u(Nbeads-2-1)+tmp1*u(Nbeads-2));
  
  /*return kappa/rodlength*(R(Nbeads-3) - R(Nbeads-4) - R(Nbeads-2) + R(Nbeads-3)
    +R(Nbeads-1) - R(Nbeads-2));*/
}

double Base::Force_Nm1(Eigen::VectorXd &u)
{

  double tmp3 = u_x(Nbeads-1-2)*u_x(Nbeads-1-1)+u_y(Nbeads-1-2)*u_y(Nbeads-1-1)+u_z(Nbeads-1-2)*u_z(Nbeads-1-1);

  return kappa/(rodlength*rodlength)*(+u(Nbeads-1-2)-tmp3*u(Nbeads-1-1));
  
  /*return kappa/rodlength*(R(Nbeads-2) - R(Nbeads-3));*/
}

/* -------------------------------------------------------------------------- */



/* -------------------------------------------------------------------------- */
/* Initialise matrix M in M * tension = C. Only call at start of simulation. */
/* -------------------------------------------------------------------------- */
void Base::set_Hhat()
{

  std::vector<T> coefficients = init_Hhat_coeffsmatrix();

  Hhat.setFromTriplets(coefficients.begin(),coefficients.end());

  solver.analyzePattern(Hhat);

  return;

}


double Base::Hhat_diag_val(int mu)
{

  double tmp1 = (tu_x(mu)*u_x(mu)+tu_y(mu)*u_y(mu)+tu_z(mu)*u_z(mu));
  tmp1 = tmp1*tmp1;

  double tmp2 = (tu_x(mu+1)*u_x(mu)+tu_y(mu+1)*u_y(mu)+tu_z(mu+1)*u_z(mu));

  tmp2 = tmp2*tmp2;

  return 2./zperp + (1./zpara-1./zperp)*(tmp1+tmp2);

}

double Base::Hhat_loweroff_val(int mu)
{

  double tmp1 = (u_x(mu-1)*u_x(mu)+u_y(mu-1)*u_y(mu)+u_z(mu-1)*u_z(mu));
  double tmp2 = (u_x(mu-1)*tu_x(mu)+u_y(mu-1)*tu_y(mu)+u_z(mu-1)*tu_z(mu));  
  double tmp3 = (tu_x(mu)*u_x(mu)+tu_y(mu)*u_y(mu)+tu_z(mu)*u_z(mu));

  return -1./zperp*tmp1-(1./zpara-1./zperp)*tmp2*tmp3;

}






/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix M. */
/* -------------------------------------------------------------------------- */
std::vector<T> Base::init_Hhat_coeffsmatrix()
{
  std::vector<T> coeffs;
  coeffs.push_back(T(0,0,Hhat_diag_val(0)));
  for (int i = 1; i < Nbeads-1; i++) {
    coeffs.push_back(T(i,i,Hhat_diag_val(i)));
    coeffs.push_back(T(i,i-1,Hhat_loweroff_val(i)));

  }

  return coeffs;

  
}


/* -------------------------------------------------------------------------- */
/* Update matrix M. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void Base::update_Hhat()
{
  // update first col by hand
  Hhat.coeffRef(0,0) = Hhat_diag_val(0);

  Hhat.coeffRef(1,0) = Hhat_loweroff_val(1);

  Hhat.coeffRef(Nbeads-2,Nbeads-2) = Hhat_diag_val(Nbeads-2);


  // update middle columns
  for (int k = 1; k < Nbeads-2; k++) {

    int count = 0;
    for (SpMat::InnerIterator it(Gmunu,k); it; ++it) {
      if (count == 0) {
	it.valueRef() = Hhat_diag_val(k);
      } else {
	it.valueRef() = Hhat_loweroff_val(k+1);
      }

      count += 1;
    }
  }
  return;
}




void Base::compute_noise()
{
  set_rhs_of_G();
  solver.factorize(Gmunu);

  dummy_for_noise =  solver.solve(rhs_of_G);

  int i = 0;
  eta_x(i) = unproj_eta_x(i) + dummy_for_noise(i)*u_x(i);
  eta_y(i) = unproj_eta_y(i) + dummy_for_noise(i)*u_y(i);
  eta_z(i) = unproj_eta_z(i) + dummy_for_noise(i)*u_z(i);

  
  for (i = 1; i < Nbeads-1; i++) {
    eta_x(i) = (unproj_eta_x(i) + dummy_for_noise(i)*u_x(i)
		 -dummy_for_noise(i-1)*u_x(i-1));
    eta_y(i) = (unproj_eta_y(i) + dummy_for_noise(i)*u_y(i)
		 -dummy_for_noise(i-1)*u_y(i-1));
    eta_z(i) = (unproj_eta_z(i) + dummy_for_noise(i)*u_z(i)
		 -dummy_for_noise(i-1)*u_z(i-1));
  }


  i = Nbeads-1;


  eta_x(i) = unproj_eta_x(i) - dummy_for_noise(i-1)*u_x(i-1);
  eta_y(i) = unproj_eta_y(i) - dummy_for_noise(i-1)*u_y(i-1);
  eta_z(i) = unproj_eta_x(i) - dummy_for_noise(i-1)*u_z(i-1);

  
  
}


void Base::compute_tension()
{
  set_rhs_of_Hhat();
  solver.factorize(Hhat);

  tension =  solver.solve(rhs_of_Hhat);
  
}


void Base::initial_integrate()
{
  
  double tmp = Delta_t/2.0;
  
  zetainv_i = zetainvs[0];
  F_x(0) += tension(0)*u_x(0);
  F_y(0) += tension(0)*u_y(0);
  F_z(0) += tension(0)*u_z(0);

  Rtmp_x(0) = R_x(0);
  Rtmp_y(0) = R_y(0);
  Rtmp_z(0) = R_z(0);
  
  R_x(0) += tmp*(zetainv_i(0,0)*F_x(0)+zetainv_i(0,1)*F_y(0)+zetainv_i(0,2)*F_z(0));
  R_y(0) += tmp*(zetainv_i(1,0)*F_x(0)+zetainv_i(1,1)*F_y(0)+zetainv_i(1,2)*F_z(0));
  R_z(0) += tmp*(zetainv_i(2,0)*F_x(0)+zetainv_i(2,1)*F_y(0)+zetainv_i(2,2)*F_z(0));
  
  for (int i = 1; i < Nbeads-1; i++) {

    Rtmp_x(i) = R_x(i);
    Rtmp_y(i) = R_y(i);
    Rtmp_z(i) = R_z(i);
    
    zetainv_i = zetainvs[i];

    F_x(i) += tension(i)*u_x(i)-tension(i-1)*u_x(i-1);
    F_y(i) += tension(i)*u_y(i)-tension(i-1)*u_y(i-1);
    F_z(i) += tension(i)*u_z(i)-tension(i-1)*u_z(i-1);

    R_x(i) += tmp*(zetainv_i(0,0)*F_x(i)+zetainv_i(0,1)*F_y(i)+zetainv_i(0,2)*F_z(i));
    R_y(i) += tmp*(zetainv_i(1,0)*F_x(i)+zetainv_i(1,1)*F_y(i)+zetainv_i(1,2)*F_z(i));
    R_z(i) += tmp*(zetainv_i(2,0)*F_x(i)+zetainv_i(2,1)*F_y(i)+zetainv_i(2,2)*F_z(i));
  }

  int i = Nbeads -1;
  zetainv_i = zetainvs[i];

  F_x(i) += -tension(i-1)*u_x(i-1);
  F_y(i) += -tension(i-1)*u_y(i-1);
  F_z(i) += -tension(i-1)*u_z(i-1);


  Rtmp_x(i) = R_x(i);
  Rtmp_y(i) = R_y(i);
  Rtmp_z(i) = R_z(i);
  
  R_x(i) += tmp*(zetainv_i(0,0)*F_x(i)+zetainv_i(0,1)*F_y(i)+zetainv_i(0,2)*F_z(i));
  R_y(i) += tmp*(zetainv_i(1,0)*F_x(i)+zetainv_i(1,1)*F_y(i)+zetainv_i(1,2)*F_z(i));
  R_z(i) += tmp*(zetainv_i(2,0)*F_x(i)+zetainv_i(2,1)*F_y(i)+zetainv_i(2,2)*F_z(i));
  
  return;
}


void Base::final_integrate()
{
  
  double tmp = Delta_t;
  
  zetainv_i = zetainvs[0];
  F_x(0) += tension(0)*u_x(0);
  F_y(0) += tension(0)*u_y(0);
  F_z(0) += tension(0)*u_z(0);
  
  R_x(0) = Rtmp_x(0) + tmp*(zetainv_i(0,0)*F_x(0)+zetainv_i(0,1)*F_y(0)+zetainv_i(0,2)*F_z(0));
  R_y(0) = Rtmp_y(0) + tmp*(zetainv_i(1,0)*F_x(0)+zetainv_i(1,1)*F_y(0)+zetainv_i(1,2)*F_z(0));
  R_z(0) = Rtmp_z(0) + tmp*(zetainv_i(2,0)*F_x(0)+zetainv_i(2,1)*F_y(0)+zetainv_i(2,2)*F_z(0));
  
  for (int i = 1; i < Nbeads-1; i++) {

    zetainv_i = zetainvs[i];

    F_x(i) += tension(i)*u_x(i)-tension(i-1)*u_x(i-1);
    F_y(i) += tension(i)*u_y(i)-tension(i-1)*u_y(i-1);
    F_z(i) += tension(i)*u_z(i)-tension(i-1)*u_z(i-1);

    R_x(i) = Rtmp_x(i) + tmp*(zetainv_i(0,0)*F_x(i)+zetainv_i(0,1)*F_y(i)+zetainv_i(0,2)*F_z(i));
    R_y(i) = Rtmp_y(i) + tmp*(zetainv_i(1,0)*F_x(i)+zetainv_i(1,1)*F_y(i)+zetainv_i(1,2)*F_z(i));
    R_z(i) = Rtmp_z(i) + tmp*(zetainv_i(2,0)*F_x(i)+zetainv_i(2,1)*F_y(i)+zetainv_i(2,2)*F_z(i));
  }

  int i = Nbeads -1;
  zetainv_i = zetainvs[i];

  F_x(i) += -tension(i-1)*u_x(i-1);
  F_y(i) += -tension(i-1)*u_y(i-1);
  F_z(i) += -tension(i-1)*u_z(i-1);
  
  R_x(i) = Rtmp_x(i) + tmp*(zetainv_i(0,0)*F_x(i)+zetainv_i(0,1)*F_y(i)+zetainv_i(0,2)*F_z(i));
  R_y(i) = Rtmp_y(i) + tmp*(zetainv_i(1,0)*F_x(i)+zetainv_i(1,1)*F_y(i)+zetainv_i(1,2)*F_z(i));
  R_z(i) = Rtmp_z(i) + tmp*(zetainv_i(2,0)*F_x(i)+zetainv_i(2,1)*F_y(i)+zetainv_i(2,2)*F_z(i));
  
  return;
}

void Base::rescale_positions(bool fromfront)
{

  if (fromfront) {
    for (int i = 1; i < Nbeads; i++) {
      R_x(i) = R_x(i-1) + rodlength*u_x(i-1);
      R_y(i) = R_y(i-1) + rodlength*u_y(i-1);
      R_z(i) = R_z(i-1) + rodlength*u_z(i-1);
    }
  } else {
    for (int i = Nbeads-2; i >= 0; i--) {
      R_x(i) = R_x(i+1) - rodlength*u_x(i-1);
      R_y(i) = R_y(i+1) - rodlength*u_y(i-1);
      R_z(i) = R_z(i+1) - rodlength*u_z(i-1);
    }
  }
  return;

}
