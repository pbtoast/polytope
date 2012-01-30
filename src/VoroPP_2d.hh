//---------------------------------Spheral++----------------------------------//
// VoroPP2d
// 
// An implemenation of the Tessellator interface that uses the 2D Voro++
// library.
//----------------------------------------------------------------------------//
#ifndef __Polytope_VoroPP_2d__
#define __Polytope_VoroPP_2d__

#include <vector>
#include <cmath>

#include "Tessellator.hh"

namespace polytope {

template<typename Real>
class VoroPP_2d: public Tessellator<2, Real> {

  //--------------------------- Public Interface ---------------------------//
public:

  // Constructors, destructor.
  VoroPP_2d(const Real xmin, const Real ymin,
            const Real xmax, const Real ymax,
            const unsigned nx = 20,
            const unsigned ny = 20,
            const Real degeneracy = 1.0e-14);
  ~VoroPP_2d();

  // Tessellate the given generators.
  virtual void tessellate(const std::vector<Real>& points,
                          Tessellation<2, Real>& mesh) const;

  // Tessellate obeying the given boundaries.
  virtual void tessellate(const std::vector<Real>& points,
                          const PLC<2, Real>& geometry,
                          Tessellation<2, Real>& mesh) const;

  // This Tessellator does not handle PLCs... yet.
  bool handlesPLCs() const { return false; }

  // Access our attributes.
  unsigned nx() const { return mNx; }
  unsigned ny() const { return mNy; }
  Real xmin() const { return mxmin; }
  Real ymin() const { return mymin; }
  Real xmax() const { return mxmax; }
  Real ymax() const { return mymax; }
  Real degeneracy() const { return std::sqrt(mDegeneracy2); }
  Real scale() const { return mScale; }

private:
  unsigned mNx, mNy;
  Real mxmin, mymin, mxmax, mymax, mDegeneracy2, mScale;

  // Forbidden methods.
  VoroPP_2d();
  VoroPP_2d(const VoroPP_2d&);
  VoroPP_2d& operator=(const VoroPP_2d&);
};

}

#endif
