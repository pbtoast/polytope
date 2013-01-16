//------------------------------------------------------------------------
// TetgenTessellator
// 
// An implemenation of the Tessellator interface that uses the Tetgen
// library by Hang Si.
// By default tetgen is built assuming double coordinates, so we only
// provide that implementation as our RealType.
//------------------------------------------------------------------------
#ifndef __Polytope_TetgenTessellator__
#define __Polytope_TetgenTessellator__

#if HAVE_TETGEN

#include <vector>
#include <cmath>

#include "Tessellator.hh"

namespace polytope {

class TetgenTessellator: public Tessellator<3, double> {
  public:
  typedef double RealType;

  // Constructor, destructor.
  TetgenTessellator();
  ~TetgenTessellator();

  // Tessellate the given generators.
  void tessellate(const std::vector<double>& points,
                  Tessellation<3, double>& mesh) const;

  // Tessellate with a bounding box representing the boundaries.
  void tessellate(const std::vector<double>& points,
                  double* low,
                  double* high,
                  Tessellation<3, double>& mesh) const;

  // Tessellate obeying the given boundaries.
  void tessellate(const std::vector<double>& points,
                  const std::vector<double>& PLCpoints,
                  const PLC<3, double>& geometry,
                  Tessellation<3, double>& mesh) const;

  // This Tessellator handles PLCs!
  bool handlesPLCs() const { return true; }
};

}

#endif
#endif
