#ifndef POLYTOPE_TESSELLATOR_HH
#define POLYTOPE_TESSELLATOR_HH

#include <vector>
#include "Tessellation.hh"
#include "PLC.hh"

namespace polytope
{

//! \class Tessellator - An abstract base class for objects that generate 
//! Voronoi and Voronoi-like tessellations for sets of points and/or 
//! geometries.
template<typename Real>
class Tessellator
{
  public:

  //! Default constructor.
  Tessellator() {}

  //! Destructor.
  virtual ~Tessellator() {}

  //! Generate a Voronoi tessellation for the given set of generator points.
  //! The coordinates of these points are stored in point-major order and 
  //! the 0th component of the ith point appears in points[Dimension*i].
  //! \param points A (Dimension*numPoints) array containing point coordinates.
  //! \param mesh This will store the resulting tessellation.
  virtual void tessellate(const std::vector<Real>& points,
                          Tessellation<Real>& mesh) const = 0;

  //! Generate a Voronoi-like tessellation for the given set of generator 
  //! points and a description of the geometry in which they exist.
  //! The coordinates of these points are stored in point-major order and 
  //! the 0th component of the ith point appears in points[Dimension*i].
  //! \param points A (Dimension*numPoints) array containing point coordinates.
  //! \param geometry A description of the geometry in Piecewise Linear Complex form.
  //! \param mesh This will store the resulting tessellation.
  virtual void tessellate(const std::vector<Real>& points,
                          const PLC<Real>& geometry,
                          Tessellation<Real>& mesh) const = 0;

  private:

  // Disallowed.
  Tessellator(const Tessellator&);
  Tessellator& operator=(const Tessellator&);

};

}

#endif
