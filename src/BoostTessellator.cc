//------------------------------------------------------------------------
// BoostTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <set>
#include <list>
#include <map>
#include <limits>
#include "float.h"

// Handy Boost stuff
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "polytope.hh" // Pulls in POLY_ASSERT

#include "within.hh"
#include "nearestPoint.hh"
#include "intersect.hh"
#include "polytope_plc_canned_geometries.hh"
#include "PLC_Boost_2d.hh"

// Fast predicate for determining colinearity of points.
extern double orient2d(double* pa, double* pb, double* pc);

// The Voronoi tools in Boost.Polygon
#include <boost/polygon/voronoi.hpp>

// We use the Boost.Geometry library to handle polygon intersections and such.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

// Include Boost multiprecision types
#include <boost/multiprecision/cpp_int.hpp>

//------------------------------------------------------------------------------
// We also need segments.
//------------------------------------------------------------------------------
namespace polytope {
template<typename CoordType>
struct PolySegment {
  Point2<CoordType> a, b;
  PolySegment(): a(), b() {}
  PolySegment(CoordType x0, CoordType y0, CoordType x1, CoordType y1): a(x0,y0), b(x1,y1) {}
  bool operator==(const PolySegment& rhs) const { return a == rhs.a and b == rhs.b; }
  bool operator!=(const PolySegment& rhs) const { return !(*this == rhs); }
  bool operator<(const PolySegment& rhs) const {
    return (a < rhs.a                ? true :
            a == rhs.a and b < rhs.b ? true :
            false);
  }
};
}

//------------------------------------------------------------------------
// Map Polytope's point class to Boost.Polygon
//------------------------------------------------------------------------
namespace boost{
namespace polygon{

typedef double fptType;
typedef int CoordHash;
typedef polytope::Point2<CoordHash> IntPoint;
typedef polytope::PolySegment<CoordHash> IntSegment;

template <>
struct geometry_concept<IntPoint> { typedef point_concept type; };
  
template <>
struct point_traits<IntPoint> {
  typedef IntPoint point_type;
  typedef CoordHash coordinate_type;
   
  static inline coordinate_type get(const point_type& point, orientation_2d orient) {
    return (orient == HORIZONTAL) ? point.x : point.y;
  }
};

template <>
struct point_mutable_traits<IntPoint> {
  typedef IntPoint point_type;
  typedef CoordHash coordinate_type;
   
  static inline void set(point_type& point, orientation_2d orient, coordinate_type value) {
    if (orient == HORIZONTAL)
      point.x = value;
    else
      point.y = value;
  }
  static inline point_type construct(coordinate_type x, coordinate_type y) {
    return point_type(x,y);
  }
};

//------------------------------------------------------------------------------
// We also need segments.
//------------------------------------------------------------------------------
template <>
struct geometry_concept<IntSegment> { typedef segment_concept type; };

template <>
struct segment_traits<IntSegment> {
  typedef CoordHash coordinate_type;
  typedef IntPoint point_type;

  static inline point_type get(const IntSegment& segment, direction_1d dir) {
    return dir.to_int() ? segment.b : segment.a;
  }
};

// //------------------------------------------------------------------------
// // Custom comparison operator
// //------------------------------------------------------------------------
// struct polytope_ulp_comparison {
//   enum Result {
//     LESS  = -1,
//     EQUAL =  0,
//     MORE  =  1
//   };
//   Result operator()(fptType a, fptType b, unsigned int maxUlps) const {
//     if (a > b) {
//       return a - b <= maxUlps ? EQUAL : LESS;
//     }
//     return b - a <= maxUlps ? EQUAL : MORE;
//   }
// };

// //------------------------------------------------------------------------
// // Custom floating point converter
// //------------------------------------------------------------------------
// struct polytope_fpt_converter {
//   template <typename T>
//   fptType operator()(const T& that) const {
//     return static_cast<fptType>(that);
//   }

//   template <size_t N>
//   fptType operator()(const typename detail::extended_int<N>& that) const {
//     fptType result = 0.0;
//     for (size_t i = 1; i <= (std::min)((size_t)3, that.size()); ++i) {
//       if (i != 1)
//         result *= static_cast<fptType>(0x100000000ULL);
//       result += that.chunks()[that.size() - i];
//     }
//     return (that.count() < 0) ? -result : result;
//   }
// };

// //------------------------------------------------------------------------
// // Custom voronoi diagram traits
// //------------------------------------------------------------------------
// struct polytope_voronoi_diagram_traits {
//   typedef fptType coordinate_type;
//   typedef voronoi_cell<coordinate_type> cell_type;
//   typedef voronoi_vertex<coordinate_type> vertex_type;
//   typedef voronoi_edge<coordinate_type> edge_type;
//   typedef struct {
//   public:
//     enum {ULPS = 128};
//     bool operator()(const vertex_type &v1, const vertex_type &v2) const {
//       return (ulp_cmp(v1.x(), v2.x(), ULPS) == polytope_ulp_comparison::EQUAL and
//               ulp_cmp(v1.y(), v2.y(), ULPS) == polytope_ulp_comparison::EQUAL);
//     }
//   private:
//      polytope_ulp_comparison ulp_cmp;
//   } vertex_equality_predicate_type;
// };

// //------------------------------------------------------------------------
// // Voronoi coordinate type traits
// //------------------------------------------------------------------------
// struct polytope_voronoi_ctype_traits {
//   typedef CoordHash int_type;
//   typedef detail::extended_int<3> int_x2_type;
//   typedef detail::extended_int<3> uint_x2_type;
//   typedef detail::extended_int<128> big_int_type;
//   typedef fptType fpt_type;
//   typedef fptType efpt_type;
//   typedef polytope_ulp_comparison ulp_cmp_type;
//   typedef polytope_fpt_converter to_fpt_converter_type;
//   typedef polytope_fpt_converter to_efpt_converter_type;
// };

} //end boost namespace
} //end polygon namespace


// // The Boost.Polygon Voronoi diagram
// typedef boost::polygon::voronoi_builder<polytope::DimensionTraits<2, RealType>::CoordHash, 
// 					boost::polygon::polytope_voronoi_ctype_traits> VB;
// typedef boost::polygon::voronoi_diagram<double, 
// 					boost::polygon::polytope_voronoi_diagram_traits> VD;

namespace polytope {

using namespace std;
using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename RealType>
BoostTessellator<RealType>::
BoostTessellator():
  Tessellator<2, RealType>() {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename RealType>
BoostTessellator<RealType>::
~BoostTessellator() {
}

//------------------------------------------------------------------------------
// Compute the QuantizedTessellation
//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
tessellateQuantized(QuantizedTessellation& result) const {

  // Sort the input points by the first element of the generator-index pair.
  // The second element provides the pre-sort generator index. Shame on you,
  // Boost, for making us do this.
  const int numGenerators = result.generators.size();
  sort(result.generators.begin(), result.generators.end());
  
  // Build ourselves the segments representing our bounding box.
  typedef PolySegment<CoordHash> IntSegment;
  vector<IntSegment> bounds(4);
  const CoordHash coordMin = QuantizedTessellation2d<CoordHash, RealType>::coordMin;
  const CoordHash coordMax = QuantizedTessellation2d<CoordHash, RealType>::coordMax;
  bounds[0] = IntSegment(coordMin, coordMin, coordMax, coordMin);
  bounds[1] = IntSegment(coordMax, coordMin, coordMax, coordMax);
  bounds[2] = IntSegment(coordMax, coordMax, coordMin, coordMax);
  bounds[3] = IntSegment(coordMin, coordMax, coordMin, coordMin);
  // sort(bounds.begin(), bounds.end());  // Not sure if this is necessary

  // Invoke the Boost.Voronoi diagram constructor
  VD voronoi;
  construct_voronoi(result.generators.begin(), result.generators.end(),
                    bounds.begin(), bounds.end(),
                    &voronoi);  
  // POLY_ASSERT2(voronoi.num_cells() == numGenerators + 4, numGenerators << " : " << voronoi.num_cells() << "!=" <<  (numGenerators + 4));

  // Read out the Voronoi topology to our intermediate QuantizedTessellation format.
  // Iterate over the edges.  Boost has organized them CCW around each generator.
  // Note here we just extract information about the point generators, and build straight
  // edges where those interface with our segment boundaries.
  result.cellEdges = vector<vector<int> >(numGenerators);
  result.edges.reserve(voronoi.num_edges());
  result.nodes.reserve(voronoi.num_vertices());
  map<IntPoint, int> node2id;
  map<std::pair<int, int>, int> edge2id;
  for (typename VD::const_cell_iterator cellItr = voronoi.cells().begin(); 
       cellItr != voronoi.cells().end(); 
       ++cellItr) {
    if (cellItr->contains_point() and cellItr->source_index() < numGenerators) {
      // POLY_ASSERT2(sortedIndex == cellItr->source_index(),
      //              sortedIndex << " != " << cellItr->source_index());
      // POLY_ASSERT(sortedIndex <  numGenerators + 4);
      POLY_ASSERT(cellItr->source_index() <  numGenerators);
      int cellIndex = result.generators[cellItr->source_index()].index;
      POLY_ASSERT(cellIndex   <  numGenerators);      

      // Start the chain walking the edges of this cell.
      const typename VD::edge_type* edge = cellItr->incident_edge();
      do {
        edge = edge->next();

        // The two vertex pointers for this edge
        // NOTE: If edge is infinite, one of these pointers is null.  This should never happen
        // since we added bounding segments.
        const typename VD::vertex_type* v0 = edge->vertex0();
        const typename VD::vertex_type* v1 = edge->vertex1();
        POLY_ASSERT(v0 and v1);

        // Finite edge.  Add the edge to the cell, and any new nodes.
        // Insert vertex 0.
        const IntPoint p0(v0->x(), v0->y());
        int old_size = node2id.size();
        const int j0 = internal::addKeyToMap(p0, node2id);
        if (j0 == old_size) {
          POLY_ASSERT(j0 == result.nodes.size());
          result.nodes.push_back(p0);
        }

        // Insert vertex 1.
        const IntPoint p1(v1->x(), v1->y());
        old_size = node2id.size();
        const int j1 = internal::addKeyToMap(p1, node2id);
        if (j1 == old_size) {
          POLY_ASSERT(j1 == result.nodes.size());
          result.nodes.push_back(p1);
        }
        
        // Now insert the edge.  Since we use oriented single edges between cells, a bit different than Boost.Polygon.
        const pair<int, int> edge = internal::hashEdge(j0, j1);
        POLY_ASSERT((edge.first == j0 and edge.second == j1) or
                    (edge.first == j1 and edge.second == j0));
        old_size = edge2id.size();
        const int e1 = internal::addKeyToMap(edge, edge2id);
        if (e1 == old_size) {
          POLY_ASSERT(e1 == result.edges.size());
          result.edges.push_back(edge);
        }
        if (edge.first == j0) {
          result.cellEdges[cellIndex].push_back(e1);
        } else {
          result.cellEdges[cellIndex].push_back(~e1);
        }
      } while (edge != cellItr->incident_edge());
    }
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class BoostTessellator<double>;

} //end polytope namespace
