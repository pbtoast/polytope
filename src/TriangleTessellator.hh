//---------------------------------Spheral++----------------------------------//
// TriangleTessellator
// 
// An implemenation of the Tessellator interface that uses the Triangle
// library by Jonathan Shewchuk.
//----------------------------------------------------------------------------//
#ifndef __Polytope_TriangleTessellator__
#define __Polytope_TriangleTessellator__

#ifdef HAVE_TRIANGLE

#include <vector>
#include <cmath>

#include "Tessellator.hh"

namespace polytope 
{

template<typename Real>
class TriangleTessellator: public Tessellator<2, Real> 
{
  public:

  // Constructor, destructor.
  TriangleTessellator();
  ~TriangleTessellator();

  // Tessellate the given generators.
  void tessellate(std::vector<Real>& points,
                  Tessellation<2, Real>& mesh) const;

  // Tessellate obeying the given boundaries.
  void tessellate(std::vector<Real>& points,
                  const PLC<2, Real>& geometry,
                  Tessellation<2, Real>& mesh) const;

  // This Tessellator handles PLCs!
  bool handlesPLCs() const { return true; }

};

}

#endif
#endif
