#ifndef POLYTOPE_TESSELLATOR_HH
#define POLYTOPE_TESSELLATOR_HH

#include <vector>
#include <float.h>
#include "Tessellation.hh"
#include "Tessellator.hh"
#include "OrphanageBase.hh"
#include "PLC.hh"
#include "ReducedPLC.hh"
#include "polytope_internal.hh"
#include "polytope_geometric_utilities.hh"
#include "polytope_tessellator_utilities.hh"

namespace polytope
{



//! \class Tessellator - An abstract base class for objects that generate 
//! Voronoi and Voronoi-like tessellations for sets of points and/or 
//! geometries.
template<int Dimension, typename RealType>
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
  virtual void tessellate(const std::vector<RealType>& points,
                          Tessellation<Dimension, RealType>& mesh) const 
  {
    error("This Tessellator does not support unbounded point tessellations.");
  }

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
  virtual void tessellate(const std::vector<RealType>& points,
                          RealType* low,
                          RealType* high,
                          Tessellation<Dimension, RealType>& mesh) const
  {
    error("This Tessellator does not support point tessellations with a bounding box.");
  }

  //! Generate a Voronoi-like tessellation for the given set of generator 
  //! points and a description of the geometry in which they exist.
  //! The coordinates of these points are stored in point-major order and 
  //! the 0th component of the ith point appears in points[Dimension*i].
  //! This default implementation issues an error explaining that the 
  //! Tessellator does not support PLCs.
  //! \param points A (Dimension*numPoints) array containing point coordinates.
  //! \param PLCpoints A (Dimension*n) array containing point coordinates for the PLC.
  //! \param geometry A description of the geometry in Piecewise Linear Complex form.
  //! \param mesh This will store the resulting tessellation.
  virtual void tessellate(const std::vector<RealType>& points,
                          const std::vector<RealType>& PLCpoints,
                          const PLC<Dimension, RealType>& geometry,
                          Tessellation<Dimension, RealType>& mesh) const
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

  virtual std::string name() const = 0;

  protected:

  //! This helper method creates a piecewise linear complex (PLC) 
  //! representing the bounding box containing the given points and 
  //! adds the corners of the bounding box to \a points.
  PLC<Dimension, RealType> boundingBox(std::vector<RealType>& points) const
  {
    // Find the minimum and maximum coordinates within the point set and 
    // their point indices.
    RealType Min[Dimension], Max[Dimension];
    for (int d = 0; d < Dimension; ++d)
    {
      Min[d] = FLT_MAX; 
      Max[d] = -FLT_MAX;
    }
    for (size_t i = 0; i < points.size()/Dimension; ++i)
    {
      for (int d = 0; d < Dimension; ++d)
      {
        Min[d] = std::min(Min[d], points[Dimension*i+d]);
        Max[d] = std::max(Max[d], points[Dimension*i+d]);
      }
    }
    // Now create the PLC and add the new generators for the corners.
    return boundingBox(Min, Max, points);
  }

  //! This helper method creates a piecewise linear complex (PLC) 
  //! representing the bounding box with the given "low" and "high"
  //! corners, and adds these corners as generator points to \a points.
  ReducedPLC<Dimension, RealType>
  boundingBox(RealType* low, RealType* high) const
  {
    ReducedPLC<Dimension, RealType> box;
    if (Dimension == 2)
    {
      // Add the new generators to points.
      box.points.push_back(low[0]);
      box.points.push_back(low[1]);

      box.points.push_back(high[0]);
      box.points.push_back(low[1]);

      box.points.push_back(high[0]);
      box.points.push_back(high[1]);

      box.points.push_back(low[0]);
      box.points.push_back(high[1]);

      // Construct the box.
      box.facets.resize(4);

      // -y face.
      box.facets[0].resize(2);
      box.facets[0][0] = 0;
      box.facets[0][1] = 1;

      // +x face.
      box.facets[1].resize(2);
      box.facets[1][0] = 1;
      box.facets[1][1] = 2;

      // +y face.
      box.facets[2].resize(2);
      box.facets[2][0] = 2;
      box.facets[2][1] = 3;

      // -x face.
      box.facets[3].resize(2);
      box.facets[3][0] = 3;
      box.facets[3][1] = 0;

      return box;
    }
    else
    {
      POLY_ASSERT(Dimension == 3);

      // Add the new generators to points.
      box.points.push_back(low[0]);
      box.points.push_back(low[1]);
      box.points.push_back(low[2]);

      box.points.push_back(high[0]);
      box.points.push_back(low[1]);
      box.points.push_back(low[2]);

      box.points.push_back(high[0]);
      box.points.push_back(high[1]);
      box.points.push_back(low[2]);

      box.points.push_back(low[0]);
      box.points.push_back(high[1]);
      box.points.push_back(low[2]);

      box.points.push_back(low[0]);
      box.points.push_back(low[1]);
      box.points.push_back(high[2]);

      box.points.push_back(high[0]);
      box.points.push_back(low[1]);
      box.points.push_back(high[2]);

      box.points.push_back(high[0]);
      box.points.push_back(high[1]);
      box.points.push_back(high[2]);

      box.points.push_back(low[0]);
      box.points.push_back(high[1]);
      box.points.push_back(high[2]);

      // Construct the box.
      box.facets.resize(6);

      // -z face.
      box.facets[0].resize(4);
      box.facets[0][0] = box.points.size()/3 - 8;
      box.facets[0][1] = box.points.size()/3 - 7;
      box.facets[0][2] = box.points.size()/3 - 6;
      box.facets[0][3] = box.points.size()/3 - 5;

      // +z face.
      box.facets[1].resize(4);
      box.facets[1][0] = box.points.size()/3 - 4;
      box.facets[1][1] = box.points.size()/3 - 3;
      box.facets[1][2] = box.points.size()/3 - 2;
      box.facets[1][3] = box.points.size()/3 - 1;

      // -x face.
      box.facets[2].resize(4);
      box.facets[2][0] = box.points.size()/3 - 8;
      box.facets[2][1] = box.points.size()/3 - 4;
      box.facets[2][2] = box.points.size()/3 - 1;
      box.facets[2][1] = box.points.size()/3 - 5;

      // +x face.
      box.facets[3].resize(4);
      box.facets[3][0] = box.points.size()/3 - 7;
      box.facets[3][1] = box.points.size()/3 - 3;
      box.facets[3][2] = box.points.size()/3 - 2;
      box.facets[3][3] = box.points.size()/3 - 6;

      // -y face.
      box.facets[4].resize(4);
      box.facets[4][0] = box.points.size()/3 - 8;
      box.facets[4][1] = box.points.size()/3 - 7;
      box.facets[4][2] = box.points.size()/3 - 3;
      box.facets[4][3] = box.points.size()/3 - 4;

      // -y face.
      box.facets[5].resize(4);
      box.facets[5][0] = box.points.size()/3 - 6;
      box.facets[5][1] = box.points.size()/3 - 2;
      box.facets[5][2] = box.points.size()/3 - 1;
      box.facets[5][3] = box.points.size()/3 - 5;

      return box;
    }
  }

  //! Return a normalized set of coordinates, also returning the bounding low/high points.
  std::vector<RealType> computeNormalizedPoints(const std::vector<RealType>& points,
                                                const std::vector<RealType>& PLCpoints,
                                                const bool computeBounds,
                                                RealType* low,
                                                RealType* high) const
  {
    POLY_ASSERT(points.size() % Dimension == 0);
    if (computeBounds) {
      double low1[Dimension], high1[Dimension];
      geometry::computeBoundingBox<Dimension, RealType>(points, true, low1, high1);
      geometry::computeBoundingBox<Dimension, RealType>(PLCpoints, true, low, high);
      for (unsigned j = 0; j != Dimension; ++j) {
        low[j] = std::min(low[j], low1[j]);
        high[j] = std::max(high[j], high1[j]);
        POLY_ASSERT(low[j] < high[j]);
      }
    }
    double boxInv[Dimension];
    for (unsigned j = 0; j != Dimension; ++j) boxInv[j] = 1.0/std::max(1e-30, high[j] - low[j]);
    const unsigned n = points.size();
    std::vector<RealType> result(n);
    for (unsigned i = 0; i != n; ++i) {
      const unsigned j = i % Dimension;
      result[i] = (points[i] - low[j])*boxInv[j];
      POLY_ASSERT(result[i] >= 0.0 and result[i] <= 1.0);
    }
    POLY_ASSERT(result.size() == points.size());
    return result;
  }

  private:

  // Disallowed.
  Tessellator(const Tessellator&);
  Tessellator& operator=(const Tessellator&);

  // Private tessellate to set bounding box and degeneracy spacing
  virtual void tessellate(const std::vector<RealType>& points,
                          const std::vector<RealType>& PLCpoints,
                          const PLC<Dimension, RealType>& geometry,
                          const RealType* low,
                          const RealType* high,
                          const RealType dx,
                          Tessellation<Dimension, RealType>& mesh) const
  {
    error("This Tessellator does not support this tessellate routine");
  }
  
  // Allow the orphanage base class to call any private tessellate(s)
  friend class OrphanageBase<Dimension, RealType>;
};

}

#endif
