#ifndef PLC_HH
#define PLC_HH

#include <vector>
typedef double Real;

//! \class PLC - A Piecewise Linear Complex in 3D, or a Planar Straight Line 
//! Graph (PSLG) in 2D.
class PLC
{
  public:

  //! An array of (Dimension*numNodes) values containing components of 
  //! node positions. The components are stored in node-major order and 
  //! the 0th component of the ith node appears in nodes[Dimension*i].
  std::vector<Real> nodes;

  //! This two-dimensional array defines the topology of the faces of the 
  //! piecewise linear complex. A face has an arbitrary number of nodes 
  //! in 3D and 2 nodes in 2D. faces[i][j] gives the index of the jth 
  //! node of the ith face.
  std::vector<std::vector<int> > faces;

  //! This array of size (Dimension*numHoles) contains components of 
  //! points identifying holes to be subtracted from the volume enclosed by 
  //! the PLC. Regions of the PLC containing a hole will not contain any 
  //! cells in their corresponding mesh. The components are stored in 
  //! point-major order and the 0th component of the ith point appears in 
  //! holes[Dimension*i].
  std::vector<Real> holes;
};


#endif
