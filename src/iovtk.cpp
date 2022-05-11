#include <stdexcept>
#include <fstream>

#include <algorithm>

#include "iovtk.hpp"

#include "double_tether.hpp"
#include "single_tether.hpp"
#include "no_tether.hpp"

/* Simple function to write a vtk file with binary data. */
template <typename pMer>
void ioVTK::writeVTKPolyData(std::string fname,const pMer &pmer)
/*============================================================================*/
/*
  Write scalar image data to a vtk (and paraview) compatible file
  (extension .vti).

  Parameters
  ----------

  fname : string
      Name of file to save with extension (either ".vtp" or ".pvtp").

  pmer : polymer object
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
	 << "<Piece NumberOfPoints=\"" << pmer.get_Nbeads()
	 << "\" NumberOfLines=\"1\">" << std::endl << "<Points>" << std::endl
	 << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
	 << std::endl;

  for (int i = 0; i < pmer.get_Nbeads(); i++) {
    myfile << pmer.atoms[i].R(0) << " " << pmer.atoms[i].R(1) << " "
	   << pmer.atoms[i].R(2) << " ";
    
  }
  myfile << std::endl << "</DataArray>" << std::endl
	 << "</Points>" << std::endl
	 << "<Lines>" << std::endl
	 << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"
	 << std::endl;

  for (int i = 0; i < pmer.get_Nbeads(); i++) 
    myfile << i << " ";
  myfile << std::endl << "</DataArray>" << std::endl
	 << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"
	 << std::endl << pmer.get_Nbeads() << std::endl
	 << "</DataArray>" << std::endl
    	 << "</Lines>" << std::endl
	 << "</Piece>" << std::endl
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

  size_t num_slashes = std::count(collectionfname.begin(), collectionfname.end(), '/');

  std::string prepath = "";
  for (int i = 0; i < num_slashes; i++)
    prepath += std::string("../");
  
  auto myfile = std::fstream(collectionfname,std::ios_base::app);
  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + collectionfname);
  }
  
  myfile << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\""
	 << " file=\"" << prepath+filename << "\"/>" << std::endl;

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

  while (stopline != "</Collection>") {
    lines.push_back(stopline);
    std::getline(oldfile,stopline);
  }
  oldfile.close();

  
  lines.erase(lines.begin());  

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
template <typename pMer>
void ioVTK::readVTKPolyData(pMer &pmer,std::string fname)
/*============================================================================*/
/*
  Read poly data to a vtk (and paraview) compatible file
  (extension .vtp).

  Parameters
  ----------

  pmer : polymer object

  fname : string
      Name of file to read from with extension (either ".vti" or ".pvti").
*/
/*============================================================================*/

{


  std::string stopline = "";

  
  auto myfile = std::fstream(fname, std::ios::in);

  

  while (stopline != "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">") {
    std::getline(myfile,stopline);
  }

  std::getline(myfile,stopline);

  
  std::stringstream ss(stopline);
  
  for (int i = 0; i < pmer.get_Nbeads(); i++) {
    ss >> pmer.atoms[i].R(0);
    ss >> pmer.atoms[i].R(1);
    ss >> pmer.atoms[i].R(2);
  }

  myfile.close();
  
}

template void ioVTK::readVTKPolyData(DoubleTether &pmer,std::string fname);
template void ioVTK::readVTKPolyData(SingleTether &pmer,std::string fname);
template void ioVTK::readVTKPolyData(NoTether &pmer,std::string fname);

template void ioVTK::writeVTKPolyData(std::string fname,const DoubleTether &pmer);
template void ioVTK::writeVTKPolyData(std::string fname,const SingleTether &pmer);
template void ioVTK::writeVTKPolyData(std::string fname,const NoTether &pmer);
