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
  /// The two coordinates for vertices are stored in the rows of a matrix
  typedef Eigen::Matrix<double, Eigen::Dynamic, 2>  coord_matrix_t;
  /// The indices of endpoints of flat panels are stored in the rows of a matrix
  typedef Eigen::Matrix<int,    Eigen::Dynamic, 2>  elem_matrix_t;

  /// Class data
  coord_matrix_t coordinates_;
  elem_matrix_t  elements_;
  bool isInitialized_;

  public:
  /**
   * Construct from matrices of coordinates and elements
   */
  BoundaryMesh(const coord_matrix_t& coords, const elem_matrix_t& elems)
  {
    coordinates_ = coords;
    elements_    = elems;
    isInitialized_ = 1;
  }


  /**
   * Construct from data files
   */
  BoundaryMesh(const std::string& filename)
  {
    loadMeshFromFile(filename);
  }


  /**
   *  This function returns the number of vertices on the mesh.
   */
  int numVertices() const;


  /**
   *  This function returns the number of elements on the mesh
   */
  int numElements() const;
  
  /**
   *  This function returns the matrix containing the mesh points.
   *
   *  @return Matrix containing the coordinates of each vertex of the 
   *  boundary mesh.
   */
  const coord_matrix_t &getMeshVertices() const;

  
  /**
   *  This function returns the matrix containing the indices of the vertices 
   *  corresponding to each element of the boundary mesh.
   */
  const elem_matrix_t &getMeshElements() const;

  
  /**
   *  This function returns the 2dvector containing the coordinates of the 
   *  vertex i.
   */
  Eigen::Vector2d getVertex(int i) const;


  /**
   *  This function returns the coordinates of the 2 vertices on the 
   *  i-th element.
   */
  std::pair<Eigen::Vector2d, Eigen::Vector2d> getElementVertices(int i) const;
  

  /**
   *  This function returns the index of the j-th vertex on the element i.
   */
  int getElementVertex(int i, int j) const;
  
  
  /**
   *  This function reads the .dat-files containing the mesh data and fills the
   *  matrices corresponding to elements and coordinates.
   *
   *  @param[in] filename string with the name of the mesh to be read.
   */
  void loadMeshFromFile(const std::string& filename);

  
  /**
   *  This function writes two .dat-files containing the mesh data: elements 
   *  and coordinates.
   *
   *  @param[in] filename string with the name of the files to be written:
   *             filename_coordinates.dat and filename_elements.dat
   */
  void writeMeshToFile(const std::string& filename);

  
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
  void readData(const std::string& filename, T & data);


}; //end class

#endif




