//---------------------------------Spheral++----------------------------------//
// VoroPP2d
// 
// An implemenation of the Tessellator interface that uses the 2D Voro++
// library.
//----------------------------------------------------------------------------//
#ifndef __Polytope_VoroPP_2d__
#define __Polytope_VoroPP_2d__

#include <vector>

#include "Tessellator.hh"
#include "container_2d.hh"

namespace polytope {

template<typename Real>
class VoroPP_2d: public Tessellator<Real> {

  //--------------------------- Public Interface ---------------------------//
public:

  // Constructors, destructor.
  VoroPP(const Real xmin, const Real ymin,
         const Real xmax, const Real ymax,
         const unsigned nx = 20,
         const unsigned ny = 20,
         const double degeneracy = 1.0e-14);
  ~VoroPP();

  // Tessellate the given generators.
  virtual void tessellate(const std::vector<Real>& points,
                          Tessellation& mesh) const;

  // Tessellate obeying the given boundaries.
  virtual void tessellate(const std::vector<Real>& points,
                          const PLC& geometry,
                          Tessellation& mesh) const;

  // Compute which sub-region the given position is in.
  void subRegion(const Real px, const Real py, 
                 unsigned& i, unsigned& j, unsigned& k) const;

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
  VoroPP();
  VoroPP(const VoroPP&);
  VoroPP& operator=(const VoroPP&);
};

}

#endif
