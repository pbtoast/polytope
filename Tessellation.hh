#ifndef POLYTOPE_TESSELLATION_HH
#define POLYTOPE_TESSELLATION_HH

#include <vector>

namespace polytope
{

//! \class Mesh - A basic descriptor class for a topologically-consistent 
//! arbitrary poly(gonal/hedral) mesh.
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
  std::vector<std::vector<int> > cells;

  //! This two-dimensional array defines the topology of the faces of the 
  //! mesh. A face has an arbitrary number of nodes in 3D and 2 nodes in 2D. 
  //! faces[i][j] gives the index of the jth node of the ith face.
  std::vector<std::vector<int> > faces;

  private:

  // Disallowed.
  Tessellation(const Tessellation&);
  Tessellation& operator=(const Tessellation&);
};

}

#endif
