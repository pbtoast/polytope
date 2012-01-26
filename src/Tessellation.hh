#ifndef POLYTOPE_TESSELLATION_HH
#define POLYTOPE_TESSELLATION_HH

#include <vector>

namespace polytope
{

//! \class Mesh - A basic descriptor class for a topologically-consistent 
//! arbitrary poly(gonal/hedral) mesh.
template<typename Real>
class Tessellation
{
  public:

  //! An array of (Dimension*numNodes) values containing components of 
  //! node positions. The components are stored in node-major order and 
  //! the 0th component of the ith node appears in nodes[Dimension*i].
  std::vector<Real> nodes;

  //! This two-dimensional array defines the cell-face topology of the 
  //! mesh. A cell has an arbitrary number of faces in 2D and 3D.
  //! cells[i][j] gives the index of the jth face of the ith cell.
  //! A negative face index indicates the actual face index is the 1's 
  //! complement of the value (~cells[i][j]) and the face is oriented
  //! witih an inward pointing normal for cells[i].
  std::vector<std::vector<int> > cells;

  //! This two-dimensional array defines the topology of the faces of the 
  //! mesh. A face has an arbitrary number of nodes in 3D and 2 nodes in 2D. 
  //! faces[i][j] gives the index of the jth node of the ith face.
  //! Nodes for a given face are arranged counterclockwise around the face
  //! viewed from the "positive" (outside) direction. 
  std::vector<std::vector<unsigned> > faces;

  //! An array of cell indices for each face, i.e., the cells that share
  //! the face.
  //! For a given cell there will be either 1 or 2 cells -- the cases with 1
  //! cell indicate a face on a boundary of the tessellation.
  std::vector<std::vector<unsigned> > faceCells;

  // This constructor should be provided automatically, but for some reason
  // clang++ on Darwin isn't doing it?
  Tessellation():
    nodes(),
    cells(),
    faces(),
    faceCells() {}

  private:

  // Disallowed.
  Tessellation(const Tessellation&);
  Tessellation& operator=(const Tessellation&);
};

}

#endif
