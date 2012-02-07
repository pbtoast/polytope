//----------------------------------------------------------------------------//
// ParallelTessellator
// 
// Provides a parallel tessellation.
//
// This class assumes the user provides:
// 1.  A serial Tessellator.
// 2.  The generators in parallel, distributed in any manner the user likes
//     so long as the generators are not degenerate: i.e., don't repeast the
//     same generator on different domains.
//
// Based on the parallel tessellation algorithm originally implemented in 
// Spheral++.
//----------------------------------------------------------------------------//
#ifndef __Polytope_ParallelTessellator__
#define __Polytope_ParallelTessellator__

#include "Tessellator.hh"

namespace polytope {

template<int Dimension, typename RealType>
class ParallelTessellator: public Tessellator<Dimension, RealType> {

  //--------------------------- Public Interface ---------------------------//
public:

  //! Constructor.
  //! \param serialTessellator A serial implementation of Tessellator.
  ParallelTessellator(const Tessellator& serialTessellator);
  ~ParallelTessellator();

  // Note the ParallelTesselator doesn't know which of the boundary treatments
  // are supported by the underlying serial Tessellator, so it's up to the user
  // to know which of these generic tessellate options they can use.

  //! Generate a Voronoi tessellation for the given set of generator points.
  //! The coordinates of these points are stored in point-major order and 
  //! the 0th component of the ith point appears in points[Dimension*i].
  //! \param points A (Dimension*numPoints) array containing point coordinates.
  //! \param mesh This will store the resulting tessellation.
//   virtual void tessellate(std::vector<RealType>& points,
//                           Tessellation<Dimension, RealType>& mesh) const;

  //! Generate a Voronoi tessellation for the given set of generator points
  //! with a bounding box specified by \a low and \a high. Here, low[i]
  //! contains the ith coordinate for the "lower-left-near" corner of the 
  //! bounding box in 2D or 3D, and high[i] contains the corresponding 
  //! opposite corner. The coordinates of these points are stored in 
  //! point-major order and the 0th component of the ith point appears in 
  //! points[Dimension*i].
  //! \param points A (Dimension*numPoints) array containing point coordinates.
  //! \param low The coordinates of the "lower-left-near" bounding box corner.
  //! \param high The coordinates of the "upper-right-far" bounding box corner.
  //! \param mesh This will store the resulting tessellation.
  virtual void tessellate(std::vector<RealType>& points,
                          RealType* low,
                          RealType* high,
                          Tessellation<Dimension, RealType>& mesh) const;

  //! Generate a Voronoi-like tessellation for the given set of generator 
  //! points and a description of the geometry in which they exist.
  //! The coordinates of these points are stored in point-major order and 
  //! the 0th component of the ith point appears in points[Dimension*i].
  //! This default implementation issues an error explaining that the 
  //! Tessellator does not support PLCs.
  //! \param points A (Dimension*numPoints) array containing point coordinates.
  //! \param geometry A description of the geometry in Piecewise Linear Complex form.
  //! \param mesh This will store the resulting tessellation.
//   virtual void tessellate(std::vector<RealType>& points,
//                           const PLC<Dimension, RealType>& geometry,
//                           Tessellation<Dimension, RealType>& mesh) const;

  //! Override this method to return true if this Tessellator supports 
  //! the description of a domain boundary using a PLC (as in the second 
  //! tessellate method, above), and false if it does not. Some algorithms 
  //! for tessellation do not naturally accommodate an explicit boundary 
  //! description, and Tessellators using these algorithms should override 
  //! this method to return false. A stub method for PLC-enabled
  //! tessellation is provided for convenience.
  //! This query mechanism prevents us from descending into the taxonomic 
  //! hell associated with elaborate inheritance hierarchies.
  virtual bool handlesPLCs() const { return mSerialTessellator.handlesPLCs(); }

private:
  Tessellator<Dimension, RealType>& mSerialTessellator;

  // Forbidden methods.
  ParallelTessellator();
  ParallelTessellator(const ParallelTessellator&);
  ParallelTessellator& operator=(const ParallelTessellator&);
};

}

#endif
