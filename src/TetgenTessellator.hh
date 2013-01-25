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
  //-------------------- Public interface --------------------
  typedef double RealType;

  // Constructor, destructor.
  TetgenTessellator(const bool directComputation = false);
  ~TetgenTessellator();

  // Tessellate the given generators (unbounded).
  void tessellate(const std::vector<RealType>& points,
                  Tessellation<3, RealType>& mesh) const;

  // Tessellate with a bounding box representing the boundaries.
  void tessellate(const std::vector<RealType>& points,
                  RealType* low,
                  RealType* high,
                  Tessellation<3, RealType>& mesh) const;

  // // Tessellate obeying the given boundaries.
  // void tessellate(const std::vector<RealType>& points,
  //                 const std::vector<RealType>& PLCpoints,
  //                 const PLC<3, RealType>& geometry,
  //                 Tessellation<3, RealType>& mesh) const;

  // This Tessellator handles PLCs!
  bool handlesPLCs() const { return false; }

  // Attributes.
  bool directComputation() const { return mDirectComputation; }
  void directComputation(const bool x) { mDirectComputation = x; }

private:
  //-------------------- Private interface --------------------
  bool mDirectComputation;

  // Internal method to compute the tessellation directly.
  void computeVoronoiNatively(const std::vector<RealType>& points,
                              Tessellation<3, RealType>& mesh) const;

  // Internal method to compute the tessellation by tetrahedralizing and computing the dual.
  void computeVoronoiThroughTetrahedralization(const std::vector<RealType>& points,
                                               Tessellation<3, RealType>& mesh) const;

  // Forbidden methods.
  TetgenTessellator();
};

}

#endif
#endif
