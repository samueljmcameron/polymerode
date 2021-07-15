
#ifndef TESTCLASS_HPP
#define TESTCLASS_HPP

#include "base.hpp"

class TestClass : public Base {
public:
  // constructor
  TestClass(int,double,double,double,double,double,double,
	    Eigen::Vector3d,int);
  ~TestClass();  

  Eigen::Matrix3d projection(int, int);
  Eigen::VectorXd compute_boldn_x(int);
  Eigen::VectorXd compute_boldn_y(int);
  Eigen::VectorXd compute_boldn_z(int);
  
};
#endif
