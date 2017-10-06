///////////////////////////////////////////////////////////////////////////////
/// \file BoundaryMesh.hpp
/// \brief This class provides functions to read the meshes used in HILBERT and
///        handles the boundary mesh data.
///////////////////////////////////////////////////////////////////////////////
#ifndef BOUNDARY_MESH_HPP
#define BOUNDARY_MESH_HPP

#include <iostream>
#include <fstream>
#include <istream> 
#include <vector>
#include <set>
#include <string>

#include <Eigen/Core>
#include <Eigen/Dense>


class BoundaryMesh
{
  private:
  /// type for mesh vertices (2d points)
  typedef Eigen::Matrix<double, Eigen::Dynamic, 2>  coord_matrix_t;
  /// type for mesh elements (storing indices of the vertices conforming the given element)
  typedef Eigen::Matrix<int,    Eigen::Dynamic, 2>  elem_matrix_t;

  /// Class data
  coord_matrix_t coordinates_;
  elem_matrix_t  elements_;
  bool isInitialized_;


  public:
  BoundaryMesh( )
  {
    isInitialized_ = 0;
  }

  BoundaryMesh(const std::string& filename)
  {
    loadMesh(filename);
  }


  int numVertices() const
  {
      return coordinates_.rows();
  }

  int numElements() const
  {
      return elements_.rows();
  }
  
  /**
   *  This function returns the matrix containing the mesh points.
   *
   *  @return Matrix containing the coordinates of each vertex of the 
   *  boundary mesh.
   */
  coord_matrix_t getMeshVertices() const
  {
    assert(isInitialized_);
    
    return coordinates_;
  }

  
  /**
   *  This function returns the matrix containing the indices of the vertices 
   *  corresponding to each element of the boundary mesh.
   */
  elem_matrix_t getMeshElements() const
  {
    assert(isInitialized_);
    
    return elements_;
  }

  
  /**
   *  This function returns the 2dvector containing the coordinates of the 
   *  vertex i.
   */
  Eigen::Vector2d getVertex(int i) const
  {
    assert(isInitialized_);
    assert(i<elements_.rows());
    
    return coordinates_.row(i);
  }


  /**
   *  This function returns the coordinates of the 2 vertices on the 
   *  i-th element.
   */
  std::pair<Eigen::Vector2d, Eigen::Vector2d> getElementVertices(int i) const
  {
    assert(isInitialized_);
    assert(i<elements_.rows());
    
    return std::make_pair<Eigen::Vector2d,
			  Eigen::Vector2d>(coordinates_.row(elements_(i,0)),
					   coordinates_.row(elements_(i,1)) );

  }
  

  /**
   *  This function returns the index of the j-th vertex on the element i.
   */
  int getElementVertex(int i, int j) const
  {
    assert(isInitialized_);
    assert(i<elements_.rows());
    assert(j<2);
    
    return elements_(i,j);
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
  void loadMesh(const std::string& filename)
  {
    readData<coord_matrix_t>(filename + "_coordinates.dat", coordinates_);
    readData<elem_matrix_t>(filename + "_elements.dat", elements_);
    // elements file has indexing starting from 1. Fix it!
    elements_ = elements_ - Eigen::MatrixXi::Ones(elements_.rows(), elements_.cols());

    // Print mesh information
    std::cout << std::string(80, '=') << std::endl
	      << std::string(27, ' ') << " READING MESH FROM FILE \n"
	      << std::string(80, '=') << std::endl;
    std::cout << "Input file : " << filename << std::endl; 
    std::cout<< "Created " << coordinates_.rows() << " vertices "
	     << "with coordinates :\n" << coordinates_ << std::endl
	     << std::endl;
    std::cout<< "Created " << elements_.rows() << " elements "
	     << ": \n" << elements_ << std::endl
	     << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    isInitialized_ = 1;
  }


  private:
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


}; //end class

#endif




