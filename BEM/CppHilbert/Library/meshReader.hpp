///////////////////////////////////////////////////////////////////////////////
/// \file meshReader.hpp
/// \brief This file provides functions to read the meshes used in HILBERT.
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <istream> 
#include <vector>
#include <set>
#include <string>

#include <Eigen/Core>
#include <Eigen/Dense>


/**
 *  This function reads a .dat-file and uses it contents to fill the
 *  corresponding matrix data.
 *
 *  @tparam T  Type of Eigen matrix corresponding to data.
 *  @param[in] filename  String with the name of the mesh to be read.
 *  @param[out] data Matrix containing the data read from file.
 */
template<typename T>
void readData(const std::string& filename, T & data)
{
  
  std::ifstream indata(filename);
  // check state
  if ( !indata ) {
    std::cout << "Could not open file '" << filename << " \n"
	      << "File does not exist!" << std::endl;
    exit( -1 );
  }
  
  std::vector<typename T::Scalar> values;
  std::string line;
  int rows = 0;
  // read every line from the stream
  while( std::getline(indata, line) )
    {
      std::stringstream dataStream(line);      
      std::string dataCell;
      // read every cell from the line that is seperated by space
      // and put it into the vector or strings
      while( std::getline(dataStream, dataCell, '	') ){
        values.push_back(std::stod(dataCell));
      }
      rows++;
    }
  int cols = values.size()/rows;
  
  data.resize(rows,cols);
  data = Eigen::Map<const Eigen::Matrix<typename T::Scalar,
					T::RowsAtCompileTime,
					T::ColsAtCompileTime,
					Eigen::RowMajor>>(values.data(),
							  rows, cols);
  
}

/**
 *  This function reads the .dat-files containing the mesh data and fills the
 *  matrices corresponding to elements and coordinates.
 *
 *  @param[in] filename string with the name of the mesh to be read.
 *  @param[out] coordinates Matrix containing the coordinates of each vertex
 *              of the boundary mesh.
 *  @param[out] elements Matrix of integers containing the indices of the
 *              vertices corresponding to each element of the boundary mesh.
 */
void readMesh(const std::string& filename, Eigen::MatrixXd & coordinates,
	      Eigen::MatrixXi & elements)
{
  readData<Eigen::MatrixXd>(filename + "_coordinates.dat", coordinates);
  readData<Eigen::MatrixXi>(filename + "_elements.dat", elements);
  // elements file has indexing starting from 1. Fix it!
  elements = elements - Eigen::MatrixXi::Ones(elements.rows(), elements.cols());

  // Print mesh information
  std::cout << std::string(80, '=') << std::endl
	    << std::string(27, ' ') << " READING MESH FROM FILE \n"
	    << std::string(80, '=') << std::endl;
  std::cout << "Input file : " << filename << std::endl; 
  std::cout<< "Created " << coordinates.rows() << " vertices "
	   << "with coordinates :\n" << coordinates << std::endl
	   << std::endl;
  std::cout<< "Created " << elements.rows() << " elements "
	   << ": \n" << elements << std::endl
	   << std::endl;
  std::cout << std::string(80, '=') << std::endl;
}

/**
 *  This function reads the .dat-files containing the mesh data and fills the
 *  matrices corresponding to elements and coordinates.
 *
 *  @param[in] filename string with the name of the mesh to be read.
 *  @param[out] coordinates Matrix containing the coordinates of each vertex
 *              of the boundary mesh.
 *  @param[out] elements Matrix of integers containing the indices of the
 *              vertices corresponding to each element of the boundary mesh.
 *  @param[out] triangles Matrix of integers containing the indices of the
 *              vertices corresponding to each triangle of the domain mesh.
 */
void readMesh(const std::string& filename, Eigen::MatrixXd & coordinates,
	      Eigen::MatrixXi & elements, Eigen::MatrixXi & triangles )
{
  readData<Eigen::MatrixXd>(filename + "_coordinates.dat", coordinates);
  readData<Eigen::MatrixXi>(filename + "_elements.dat", elements);
  readData<Eigen::MatrixXi>(filename + "_triangles.dat", triangles);
  // elements file has indexing starting from 1. Fix it!
  elements = elements - Eigen::MatrixXi::Ones(elements.rows(), elements.cols());
  triangles = triangles - Eigen::MatrixXi::Ones(triangles.rows(), triangles.cols());

  // Print mesh information
  std::cout << std::string(80, '=') << std::endl
	    << std::string(27, ' ') << " READING MESH FROM FILE \n"
	    << std::string(80, '=') << std::endl;
  std::cout << "Input file : " << filename << std::endl; 
  std::cout<< "Created " << coordinates.rows() << " vertices. "
	   << "with coordinates :\n" << coordinates << std::endl
	   << std::endl;
  std::cout<< "Created " << elements.rows() << " boundary elements. "
	   << ": \n" << elements << std::endl
	   << std::endl;
    std::cout<< "Created " << triangles.rows() << " triangles. "
	   << ": \n" << triangles << std::endl
	   << std::endl;
  std::cout << std::string(80, '=') << std::endl;
}
