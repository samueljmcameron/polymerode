
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>

template<class T>
class Matrix
{
  size_t rows;
  size_t cols;
  std::vector<T> m_Data;
  /*
  const double get_comp(size_t i, size_t j, size_t mu)
  {
    return (*this)(i,j)(mu);
  };
  */
  
public:
  Matrix(size_t rin=0, size_t cin=0) :
    rows(rin), cols(cin) {
    m_Data.resize(rows*cols);
  };
  
  T& operator()(size_t y, size_t x)
  {
    return m_Data[y * cols + x];
    
  };

  //  friend std::ostream& operator << (std::ostream&,
  //				    Matrix&);

  void resize(size_t r, size_t c)
  {
    rows = r;
    cols = c;
    m_Data.resize(rows*cols);
    return;
  };
  

};


/*
std::ostream& operator << (std::ostream& out, Matrix& data)
{

  out << std::endl << std::endl;;
  for (size_t i =0; i < data.rows; i++) {
    size_t j = 0;
    for (j = 0; j < data.cols-1; j++) {
      out << data.get_comp(i,j,0) << " " << data.get_comp(i,j,1)
	  << " " << data.get_comp(i,j,2) << ",\t";
    }
    
    out << data.get_comp(i,j,0) << " " << data.get_comp(i,j,1)
	<< " " << data.get_comp(i,j,2) << std::endl << std::endl;
  }
  
  return out;
}

*/



#endif
