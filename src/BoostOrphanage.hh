//------------------------------------------------------------------------
// BoostOrphanage
// 
// Descended class for implementing a cell adoption algorithm. Uses
// Boost.Geometry to construct sub-regions of space on which bounded
// tessellations are built. The sub-regions include the orphaned cell
// pieces and their neighbors and are desigend to satisfy the Voronoi
// principle local to the cell orphans.
//
// The Algorithm:
// Make a sub-tessellation using the generators neighboring the orphan.
// Construct a PLC boundary for the sub-tessellation using the geometry 
// obtained by union-ing the orphan with its neighboring cells. By using 
// a PLC boundary from all the orphan neighbors, we can ensure the union 
// gives a contiguous geometry with no voids.
//------------------------------------------------------------------------
#ifndef __Polytope_BoostOrphanage__
#define __Polytope_BoostOrphanage__

#if HAVE_BOOST

#include <set>
#include <map>
#include <vector>

#include "OrphanageBase.hh"
#include "polytope_tessellator_utilities.hh"
#include "Point.hh"

// We use the Boost.Geometry library to handle polygon intersections and such.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>

namespace polytope
{


template<typename RealType>
class BoostOrphanage: public OrphanageBase<2, RealType>
{
public:

  // Polytope typedefs
  typedef int64_t CoordHash;
  typedef Point2<CoordHash> IntPoint;  

  // Boost typedefs
  typedef boost::geometry::model::polygon<IntPoint, false> BGpolygon;
  typedef boost::geometry::model::ring   <IntPoint, false> BGring;
  typedef boost::geometry::model::multi_polygon<BGpolygon> BGmulti_polygon;

  // Constructor, Destructor
  BoostOrphanage(const Tessellator<2, RealType>* tessellator);
  ~BoostOrphanage();

  // Re-tessellate the area in orphaned cells by modifying the existing cell rings
  void adoptOrphans(const std::vector<RealType>& points,
                    const RealType* low,
                    const RealType* high,
                    const RealType dx,
                    std::vector<BGring>& cellRings,
                    std::vector<BGring>& orphans) const;
};

} //end polytope namespace

#endif
#endif
