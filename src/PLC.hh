#ifndef POLYTOPE_PLC_HH
#define POLYTOPE_PLC_HH

#include <vector>

namespace polytope
{

//! \class PLC - A Piecewise Linear Complex in 3D, or a Planar Straight Line 
//! Graph (PSLG) in 2D.
template<typename Real>
class PLC
{
  public:

  //! This two-dimensional array defines the topology of the facets of the 
  //! piecewise linear complex in terms of connections to generating points. 
  //! A facet has an arbitrary number of points in 3D and 2 points in 2D. 
  //! facets[i][j] gives the index of the jth generating point of the ith 
  //! facet.
  std::vector<std::vector<int> > facets;

  //! This array of size (Dimension*numHoles) contains components of 
  //! points identifying holes to be subtracted from the volume enclosed by 
  //! the PLC. Regions of the PLC containing a hole will not contain any 
  //! cells in their corresponding mesh. The components are stored in 
  //! point-major order and the 0th component of the ith point appears in 
  //! holes[Dimension*i].
  std::vector<Real> holes;

  //! Returns true if this PLC is empty, false otherwise.
  bool empty() const
  {
    return (facets.empty() and holes.empty());
  }
};

}

#endif
