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

template<typename RealType>
class TetgenTessellator: public Tessellator<3, RealType> 
{
  public:

  // Constructor, destructor.
  TetgenTessellator();
  ~TetgenTessellator();

  // Tessellate the given generators.
  void tessellate(const std::vector<RealType>& points,
                  Tessellation<3, RealType>& mesh) const;

  // Tessellate obeying the given boundaries.
  void tessellate(const std::vector<RealType>& points,
                  const PLC<3, RealType>& geometry,
                  Tessellation<3, RealType>& mesh) const;

  // This Tessellator handles PLCs!
  bool handlesPLCs() const { return true; }

};

}

#endif
#endif
