//---------------------------------Spheral++----------------------------------//
// Triangle
// 
// An implemenation of the Tessellator interface that uses the Triangle
// library by Jonathan Shewchuk.
//----------------------------------------------------------------------------//
#ifndef __Polytope_Triangle__
#define __Polytope_Triangle__

#include <vector>
#include <cmath>

#include "Tessellator.hh"

namespace polytope {

template<typename Real>
class Triangle: public Tessellator<Real> 
{
  public:

  // Constructor, destructor.
  Triangle();
  ~Triangle();

  // Tessellate the given generators.
  void tessellate(const std::vector<Real>& points,
                  Tessellation<Real>& mesh) const;

  // Tessellate obeying the given boundaries.
  void tessellate(const std::vector<Real>& points,
                  const PLC<Real>& geometry,
                  Tessellation<Real>& mesh) const;

};

}

#endif
