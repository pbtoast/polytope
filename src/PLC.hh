#ifndef POLYTOPE_PLC_HH
#define POLYTOPE_PLC_HH

#include <vector>
#include <iostream>

namespace polytope
{

//! \class PLC - A Piecewise Linear Complex in 3D, or a Planar Straight Line 
//! Graph (PSLG) in 2D.
template<int Dimension, typename Real>
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

  //! Returns true if this PLC is valid (at first glance), false if it 
  //! is obviously invalid. This is not a rigorous check!
  bool valid() const
  {
    if (Dimension == 2)
    {
      // In 2D all facets must have at least 2 points.
      for (int f = 0; f < facets.size(); ++f)
      {
        if (facets[f].size() < 2)
          return false;
      }
    }
    else if (Dimension == 3)
    {
      // In 3D all facets must have at least 3 points.
      for (int f = 0; f < facets.size(); ++f)
      {
        if (facets[f].size() < 3)
          return false;
      }
    }
    return true;
  }

  //! output operator.
  friend std::ostream& operator<<(std::ostream& s, const PLC& plc)
  {
    s << "PLC (" << Dimension << "D):" << std::endl;
    s << plc.facets.size() << " facets:" << std::endl;
    for (int f = 0; f < plc.facets.size(); ++f)
    {
      s << " " << f << ": (";
      for (int p = 0; p < plc.facets[f].size(); ++p)
      {
        if (p < plc.facets[f].size()-1)
          s << plc.facets[f][p] << ", ";
        else
          s << plc.facets[f][p];
      }
      s << ")" << std::endl;
    }
    s << std::endl;
    s << plc.holes.size() << " holes:" << std::endl;
    for (int h = 0; h < plc.holes.size()/Dimension; ++h)
    {
      s << " " << h << ": "; 
      if (Dimension == 2)
        s << "(" << plc.holes[2*h] << ", " << plc.holes[2*h+1] << ")" << std::endl;
      else
      {
        ASSERT(Dimension == 3);
        s << "(" << plc.holes[3*h] << ", " << plc.holes[3*h+1] << ", " << plc.holes[3*h+2] << ")" << std::endl;
      }
    }

    return s;
  }

};

}

#endif
