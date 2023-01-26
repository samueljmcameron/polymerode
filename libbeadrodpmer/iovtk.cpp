/* -*- c++ -*- */

#include <stdexcept>
#include <fstream>

#include <algorithm>

#include "iovtk.hpp"


/* Simple function to write a vtk file with binary data. */

namespace BeadRodPmer {
void ioVTK::writeVTKPolyData(std::string fname,
			     const Eigen::Ref<Eigen::Matrix3Xd> xs,
			     bool connect_atoms)
/*============================================================================*/
/*
  Write scalar image data to a vtk (and paraview) compatible file
  (extension .vti).

  Parameters
  ----------

  fname : string
      Name of file to save with extension (either ".vtp" or ".pvtp").

  xs : beads to save
*/
/*============================================================================*/

{
  
  auto myfile = std::fstream(fname, std::ios::out);

  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + fname);
  }

  
  myfile << "<?xml version=\"1.0\"?>" << std::endl
	 << "<VTKFile type=\"PolyData\"  version=\"1.0\""
	 << " byte_order=\"LittleEndian\">" << std::endl
	 << "<PolyData>" << std::endl
	 << "<Piece NumberOfPoints=\"" << xs.cols()
	 << "\" NumberOfLines=\"1\">" << std::endl << "<Points>" << std::endl
	 << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
	 << std::endl;

  for (int i = 0; i < xs.cols(); i++) {
    myfile << xs.col(i)(0) << " " << xs.col(i)(1) << " "
	   << xs.col(i)(2) << " ";
    
  }
  myfile << std::endl << "</DataArray>" << std::endl
	 << "</Points>" << std::endl;

  if (connect_atoms) {
	myfile << "<Lines>" << std::endl
	       << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"
	       << std::endl;

	for (int i = 0; i < xs.cols(); i++) 
	  myfile << i << " ";
	myfile << std::endl << "</DataArray>" << std::endl
	       << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"
	       << std::endl << xs.cols() << std::endl
	       << "</DataArray>" << std::endl
	       << "</Lines>" << std::endl;
  }
  myfile << "</Piece>" << std::endl
	 << "</PolyData>" << std::endl
	 << "</VTKFile>" << std::endl;    
  myfile.close();
  
}


void ioVTK::writeVTKcollectionHeader(const std::string fname)
{
  auto myfile = std::ofstream(fname);
  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + fname);
  }
  
  myfile << "<?xml version=\"1.0\"?>" << std::endl
	 << "<VTKFile type=\"Collection\"  version=\"1.0\""
	 << " byte_order=\"LittleEndian\">" << std::endl
	 << "<Collection>" << std::endl;


  myfile.close();
  

}

void ioVTK::writeVTKcollectionMiddle(const std::string collectionfname,
				     const std::string filename,
				     const double time)
{


  size_t firstslash = filename.find_last_of("\\/");
  std::string file_no_path = filename;
  if (firstslash != std::string::npos) {
    file_no_path = filename.substr(firstslash+1);
  }

  auto myfile = std::fstream(collectionfname,std::ios_base::app);
  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + collectionfname);
  }
  
  myfile << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\""
	 << " file=\"" << file_no_path << "\"/>" << std::endl;

  myfile.close();
  
}


void ioVTK::writeVTKcollectionFooter(const std::string fname)
{
  auto myfile = std::fstream(fname, std::ios_base::app);
  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + fname);
  }
  
  myfile << "</Collection>" << std::endl
	 << "</VTKFile>";
  
  myfile.close();
}

void ioVTK::restartVTKcollection(const std::string fname)
/* restart the collection by copying the current file (up to the collection
   end) into a new file called "restart" + oldfname. */
   
{

  auto oldfile = std::ifstream(fname);
  if (not oldfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + fname);
  }



  std::vector<std::string> lines;
  std::string stopline = "";
  std::getline(oldfile,stopline);
  lines.push_back(stopline);

  while (std::getline(oldfile,stopline) && stopline != "</Collection>") {
    lines.push_back(stopline);
  }
  oldfile.close();

  if (stopline != "</Collection>") {
    throw std::runtime_error("Invalid restart file." );
  }



  auto newfile = std::ofstream(fname);
  if (not newfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + fname);
  }

  for (const auto &line : lines)
    newfile << line << std::endl;;
  


  newfile.close();

  return;
}


/* Simple function to write a vtk file with binary data. */

void ioVTK::readVTKPolyData(Eigen::Ref<Eigen::Matrix3Xd> xs,std::string fname)
/*============================================================================*/
/*
  Read poly data to a vtk (and paraview) compatible file
  (extension .vtp).

  Parameters
  ----------

  xs : beads to overwrite

  fname : string
      Name of file to read from with extension (either ".vti" or ".pvti").
*/
/*============================================================================*/

{


  std::string stopline = "";

  
  auto myfile = std::fstream(fname, std::ios::in);
  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + fname);
  }

  

  while (stopline != "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">") {
    std::getline(myfile,stopline);
  }

  std::getline(myfile,stopline);

  
  std::stringstream ss(stopline);
  
  for (int i = 0; i < xs.cols(); i++) {
    ss >> xs.col(i)(0);
    ss >> xs.col(i)(1);
    ss >> xs.col(i)(2);
  }

  myfile.close();
  
}

}
