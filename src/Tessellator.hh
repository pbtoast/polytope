#ifndef POLYTOPE_TESSELLATOR_HH
#define POLYTOPE_TESSELLATOR_HH

#include <vector>
#include "Tessellation.hh"
#include "Tessellator.hh"
#include "PLC.hh"

namespace polytope
{

//! \class Tessellator - An abstract base class for objects that generate 
//! Voronoi and Voronoi-like tessellations for sets of points and/or 
//! geometries.
template<int Dimension, typename Real>
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
                          Tessellation<Dimension, Real>& mesh) const = 0;

  //! Generate a Voronoi-like tessellation for the given set of generator 
  //! points and a description of the geometry in which they exist.
  //! The coordinates of these points are stored in point-major order and 
  //! the 0th component of the ith point appears in points[Dimension*i].
  //! This default implementation issues an error explaining that the 
  //! Tessellator does not support PLCs.
  //! \param points A (Dimension*numPoints) array containing point coordinates.
  //! \param geometry A description of the geometry in Piecewise Linear Complex form.
  //! \param mesh This will store the resulting tessellation.
  virtual void tessellate(const std::vector<Real>& points,
                          const PLC<Dimension, Real>& geometry,
                          Tessellation<Dimension, Real>& mesh) const
  {
    error("This Tessellator does not support boundaries specified as PLCs.");
  }

  //! Override this method to return true if this Tessellator supports 
  //! the description of a domain boundary using a PLC (as in the second 
  //! tessellate method, above), and false if it does not. Some algorithms 
  //! for tessellation do not naturally accommodate an explicit boundary 
  //! description, and Tessellators using these algorithms should override 
  //! this method to return false. A stub method for PLC-enabled
  //! tessellation is provided for convenience.
  //! This query mechanism prevents us from descending into the taxonomic 
  //! hell associated with elaborate inheritance hierarchies.
  virtual bool handlesPLCs() const = 0;

  private:

  // Disallowed.
  Tessellator(const Tessellator&);
  Tessellator& operator=(const Tessellator&);

};

}

#endif
