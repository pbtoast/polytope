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
class TetgenTessellator: public Tessellator<Real> 
{
  public:

  // Constructor, destructor.
  TetgenTessellator();
  ~TetgenTessellator();

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
#endif
