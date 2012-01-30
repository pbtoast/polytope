//---------------------------------Spheral++----------------------------------//
// TetgenTessellator
// 
// An implemenation of the Tessellator interface that uses the Tetgen
// library by Hang Si.
//----------------------------------------------------------------------------//
#ifndef __Polytope_TetgenTessellator__
#define __Polytope_TetgenTessellator__

#ifdef HAVE_TETGEN

#include <vector>
#include <cmath>

#include "Tessellator.hh"

namespace polytope {

template<typename Real>
class TetgenTessellator: public Tessellator<3, Real> 
{
  public:

  // Constructor, destructor.
  TetgenTessellator();
  ~TetgenTessellator();

  // Tessellate the given generators.
  void tessellate(const std::vector<Real>& points,
                  Tessellation<3, Real>& mesh) const;

  // Tessellate obeying the given boundaries.
  void tessellate(const std::vector<Real>& points,
                  const PLC<3, Real>& geometry,
                  Tessellation<3, Real>& mesh) const;

  // This Tessellator handles PLCs!
  bool handlesPLCs() const { return true; }

};

}

#endif
#endif
