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
// Unbounded tessellation
//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
tessellate(const vector<RealType>& points,
           Tessellation<2, RealType>& mesh) const {
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() > 0);
  POLY_ASSERT(points.size() % 2 == 0);

  // Check for collinear generators
  const bool collinear = geometry::collinear<2, RealType>(points, 1.0e-10);
  
  // Use the appropriate cell node routine
  QuantizedTessellation2d<CoordHash, double> quantmesh(points);
  if (collinear) {
    this->computeCollinearQuantTessellation(quantmesh);
  } else {
    this->computeQuantTessellation(quantmesh);
  }

  // Copy the QuantTessellation to the output.
  quantmesh.fillTessellation(mesh);
}

//------------------------------------------------------------------------------
// Box bounded tessellation
//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
tessellate(const vector<RealType>& points,
           RealType* low,
           RealType* high,
           Tessellation<2, RealType>& mesh) const {
  // // Pre-conditions
  // POLY_ASSERT(low != 0 and high != 0);
  // POLY_ASSERT(low[0] <= high[0] and low[1] <= high[1]);
  
  // // Build a PLC with the bounding box, and then use the PLC method.
  // ReducedPLC<2, RealType> box = plc_box<2, RealType>(low, high);
  // this->tessellate(points, box, mesh);
}

//------------------------------------------------------------------------------
// PLC bounded tessellation
//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const vector<RealType>& plcPoints,
           const PLC<2, RealType>& geometry,
           Tessellation<2, RealType>& mesh) const {
  // ReducedPLC<2, RealType> boundary;
  // boundary.facets = geometry.facets;
  // boundary.holes  = geometry.holes;
  // boundary.points = plcPoints;
  // this->tessellate(points, boundary, mesh);
}

//------------------------------------------------------------------------------
// Reduced PLC bounded tessellation
//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const ReducedPLC<2, RealType>& geometry,
           Tessellation<2, RealType>& mesh) const {
  // // Pre-conditions
  // POLY_ASSERT(mesh.empty());
  // POLY_ASSERT(not points.empty());
  // POLY_ASSERT(not geometry.points.empty());

  // // Initialize quantized coordinate system
  // mCoords.initialize(geometry.points, mDegeneracy);
  // mCoords.points = geometry.points;
  // mCoords.facets = geometry.facets;
  // mCoords.holes  = geometry.holes;

  // // Check that the coordinates don't overflow 32-bit integers. Boost.Voronoi
  // // expects integer input and defaults to 32 bits. This can be extended to
  // // greater precision using custom trait classes but has not been tested
  // // yet.  -DPS 01/19/2015
  // POLY_VERIFY2(mCoords.coordMax() < std::numeric_limits<int32_t>::max(),
  //              "BoostTessellator Error: the specified degeneracy spacing "
  //              << mDegeneracy << " is too small to be represented using "
  //              << "32-bit integers given the scale of the input generators. "
  //              << "Please reduce the degeneracy spacing to tessellate.");

  // // Check for collinear generators
  // const bool collinear = geometry::collinear<2, RealType>(points, 1.0e-10);
  // const unsigned numGenerators = points.size()/2;

  // // Use the appropriate cell node routine
  // vector<vector<unsigned> > cellNodes;
  // map<int, PointType> id2node;
  // vector<unsigned> infNodes;
  // if (collinear) 
  // {
  //   this->computeCellNodesCollinear(points, cellNodes, id2node, infNodes);
  // }
  // else 
  // {
  //   this->computeCellNodes(points, cellNodes, id2node, infNodes);
  // }
  // POLY_ASSERT(cellNodes.size() == numGenerators);

  // vector<ReducedPLC<2, CoordType> > orphans;
  // vector<ReducedPLC<2, CoordType> > cells(numGenerators);

  // // Convert boundary to proper point type to compute intersections
  // ReducedPLC<2, CoordType> boundary;
  // boundary.facets = geometry.facets;
  // boundary.holes  = geometry.holes;
  // boundary.points.resize(geometry.points.size());
  // for (int j = 0; j < geometry.points.size()/2; ++j) {
  //   const PointType pp = BTT::quantize(mCoords, &geometry.points[2*j]);
  //   boundary.points[2*j  ] = pp.x;
  //   boundary.points[2*j+1] = pp.y;
  // }

  // // Compute cell-boundary intersections
  // for (int i = 0; i < numGenerators; ++i) {
  //   const PointType generator = BTT::quantize(mCoords, &points[2*i]);
  //   ReducedPLC<2, CoordType> cell = plcOfCell<CoordType>(cellNodes[i], id2node);

  //   //TODO
  //   //TODO implement the bi-infinite edge fix
  //   //TODO
    
  //   POLY_ASSERT2(not BG::boost_intersects<CoordType>(cell),
  //       	 "Cell " << i << " self-intersects before clipping"
  //       	 << endl << cell);
  //   cells[i] = BG::boost_clip<CoordType>(boundary, cell, generator, orphans);
  //   POLY_ASSERT(not BG::boost_intersects<CoordType>(cells[i]));
  // }

  // if (not orphans.empty()) {
  //   cerr << "Orphans detected. Taking no actions." << endl;
  // }

  // // Input nodes and construct the final mesh topology
  // this->constructBoundedTopology(points, geometry, cells, mesh);
}


//------------------------------------------------------------------------------
// Private routines called by tessellate:
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Compute the QuantizedTessellation
//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
computeQuantTessellation(QuantizedTessellation2d<CoordHash, RealType>& result) const{

  // Sort the input points by the first element of the generator-index pair.
  // The second element provides the pre-sort generator index. Shame on you,
  // Boost, for making us do this.
  const int numGenerators = result.generators.size();
  sort(result.generators.begin(), result.generators.end());
  
  // Build ourselves the segments representing our bounding box.
  typedef PolySegment<CoordHash> IntSegment;
  vector<IntSegment> bounds(4);
  const CoordHash coordMin = std::numeric_limits<CoordHash>::min()/4;
  const CoordHash coordMax = std::numeric_limits<CoordHash>::max()/4;
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
// Specialized compute QuantizedTessellation for collinear generators
//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
computeCollinearQuantTessellation(QuantizedTessellation2d<CoordHash, RealType>& result) const {

  // // Call the 1d routine for projecting a line of points
  // vector<RealPoint> nodes;
  // const RealPoint center = mCoords.center();
  // constructCells1d(points, &center.x, mCoords.infiniteRadius(), cellNodes, nodes);
  // POLY_ASSERT(cellNodes.size() == points.size()/2);
  // POLY_ASSERT(nodes.size() == points.size());

  // // Quantize nodes and assign indices
  // std::set<PointType> uniqueNodes;
  // infNodes.resize(nodes.size());
  // for (unsigned i = 0; i < nodes.size(); ++i) {
  //   const PointType p = BTT::quantize(mCoords, &(nodes[i]).x);
  //   uniqueNodes.insert(p);
  //   POLY_ASSERT(uniqueNodes.size() == i+1);
  //   id2node[i] = p;
  //   infNodes[i] = i;   // All nodes are projected inf nodes
  // }
  // POLY_ASSERT(uniqueNodes.size() == nodes.size());
}
//------------------------------------------------------------------------------

// //------------------------------------------------------------------------------
// template<typename RealType>
// void
// BoostTessellator<RealType>::
// constructBoundedTopology(const vector<RealType>& points,
// 			 const ReducedPLC<2, RealType>& geometry,
// 			 const vector<ReducedPLC<2, CoordType> >& cellRings,
// 			 Tessellation<2, RealType>& mesh) const {
//   // Pre-conditions
//   POLY_ASSERT(mesh.empty());
  
//   // Now build the unique mesh nodes and cell info.
//   const unsigned numGenerators = cellRings.size();
//   // map<PointType, int> point2node;
//   map<IntPoint, int> point2node;
//   map<EdgeHash, int> edgeHash2id;
//   map<int, vector<int> > edgeCells;
//   int i, j, k, iedge;
//   mesh.cells = std::vector<std::vector<int> >(numGenerators);
//   for (i = 0; i != numGenerators; ++i) { 
//     POLY_ASSERT(cellRings[i].facets.size() > 2);
//     for (unsigned ifacet = 0; ifacet < cellRings[i].facets.size(); ++ifacet) {
//       POLY_ASSERT(cellRings[i].facets[ifacet].size() == 2);
//       const int i1 = cellRings[i].facets[ifacet][0];
//       const int i2 = cellRings[i].facets[ifacet][1];
//       POLY_ASSERT(i1 != i2);
//       const PointType p1 = PointType(cellRings[i].points[2*i1],
// 				     cellRings[i].points[2*i1+1]);
//       const PointType p2 = PointType(cellRings[i].points[2*i2],
// 				     cellRings[i].points[2*i2+1]);
//       POLY_ASSERT(p1 != p2);
//       const IntPoint ip1 = mCoords.quantize(&p1.x);
//       const IntPoint ip2 = mCoords.quantize(&p2.x);
//       // j = internal::addKeyToMap(p1, point2node);
//       // k = internal::addKeyToMap(p2, point2node);
//       j = internal::addKeyToMap(ip1, point2node);
//       k = internal::addKeyToMap(ip2, point2node);
//       if (j != k) {
//       POLY_ASSERT(j != k);
//       iedge = internal::addKeyToMap(internal::hashEdge(j,k), edgeHash2id);
//       edgeCells[iedge].push_back(j < k ? i : ~i);
//       POLY_ASSERT2(edgeCells[iedge].size() == 1 or edgeCells[iedge][0]*edgeCells[iedge][1] <= 0,
//                    "BLAGO: " << iedge << " " << j << " " << k << " " << edgeCells[iedge][0] 
//                    << " " << edgeCells[iedge][1]);
//       mesh.cells[i].push_back(j < k ? iedge : ~iedge);
//       }
//     }
//     POLY_ASSERT(mesh.cells[i].size() >= 3);
//   }
//   POLY_ASSERT(edgeCells.size() == edgeHash2id.size());
  
//   // Fill in the mesh nodes.
//   RealPoint node;
//   mesh.nodes = std::vector<RealType>(2*point2node.size());
//   // for (typename std::map<PointType, int>::const_iterator itr = point2node.begin();
//   for (typename std::map<IntPoint, int>::const_iterator itr = point2node.begin();
//        itr != point2node.end(); 
//        ++itr) {
//     // const RealPoint p = BTT::dequantize(mCoords, itr->first);
//     const RealPoint p = mCoords.dequantize(&(itr->first).x);
//     i = itr->second;
//     POLY_ASSERT(i < mesh.nodes.size()/2);
//     // if (not BG::boost_within(p, geometry)) {
//     //   RealPoint result;
//     //   RealType dist = nearestPoint(&node.x, 
//     // 				      geometry.points.size()/2,
//     // 				      &geometry.points[0],
//     // 				      geometry,
//     // 				      &result.x);
//     //   POLY_ASSERT(dist < 1.0e-10);
//     //   node = result;
//     // }
//     node = p;
//     mesh.nodes[2*i]   = node.x;
//     mesh.nodes[2*i+1] = node.y;
//   }
  
//   // Fill in the mesh faces.
//   mesh.faces = std::vector<std::vector<unsigned> >(edgeHash2id.size());
//   for (typename std::map<EdgeHash, int>::const_iterator itr = edgeHash2id.begin();
//        itr != edgeHash2id.end(); 
//        ++itr) {
//     const EdgeHash& ehash = itr->first;
//     i = itr->second;
//     POLY_ASSERT(i < mesh.faces.size());
//     POLY_ASSERT(mesh.faces[i].size() == 0);
//     mesh.faces[i].push_back(ehash.first);
//     mesh.faces[i].push_back(ehash.second);
//   }

//   // Fill in the mesh faceCells.
//   mesh.faceCells = std::vector<std::vector<int> >(edgeHash2id.size());
//   for (i = 0; i != mesh.faces.size(); ++i) {
//     if (not(edgeCells[i].size() == 1 or edgeCells[i].size() == 2)) {
//       const int n1 = mesh.faces[i][0], n2 = mesh.faces[i][1];
//       std::cerr << "Blago! " << i << " "
// 		<< edgeCells[i].size() << " : " << n1 << " " << n2 << " : ("
//                 << mesh.nodes[2*n1] << " " << mesh.nodes[2*n1 + 1] << ") ("
//                 << mesh.nodes[2*n2] << " " << mesh.nodes[2*n2 + 1] << ")" << std::endl;
//       for (j = 0; j != edgeCells[i].size(); ++j) {
//         std::cerr << " --> " << edgeCells[i][j] << " " 
//                   << points[2*edgeCells[i][j]] << " " 
//                   << points[2*edgeCells[i][j]+1] << std::endl;
//       }
//     }
//     POLY_ASSERT(edgeCells[i].size() == 1 or edgeCells[i].size() == 2);
//     mesh.faceCells[i] = edgeCells[i];
//   }

//   // Post-conditions
//   POLY_ASSERT(mesh.faceCells.size() == mesh.faces.size());
//   POLY_ASSERT(mesh.cells.size()     == numGenerators    );
// }
// //------------------------------------------------------------------------------





// //------------------------------------------------------------------------------
// // Private tessellate routines
// //------------------------------------------------------------------------------






// //------------------------------------------------------------------------------
// template<typename RealType>
// void
// BoostTessellator<RealType>::
// tessellate(const std::vector<RealType>& points,
//            const ReducedPLC<2, CoordHash>& intGeometry,
//            const QuantizedCoordinates<2, RealType>& coords,
//            vector<ReducedPLC<2, CoordHash> >& intCells) const {
//   // Pre-conditions
//   POLY_ASSERT(not intGeometry.empty());
//   POLY_ASSERT(not points.empty() > 0);
//   POLY_ASSERT(points.size() % 2 == 0);
//   POLY_ASSERT(not coords.empty());

//   // The Quantized coordinates
//   mCoords = coords;
  
//   const bool collinear = geometry::collinear<2, RealType>(points, mDegeneracy);
//   vector<vector<unsigned> > cellNodes;
//   map<int, PointType> id2node;
//   vector<unsigned> infNodes;
  
//   // Use the appropriate cell node routine
//   if (collinear)
//   {
//     this->computeCellNodesCollinear(points, cellNodes, id2node, infNodes);
//   }
//   else 
//   {
//     this->computeCellNodes(points, cellNodes, id2node, infNodes);
//   }
//   POLY_ASSERT(cellNodes.size() == points.size()/2);

//   // vector<ReducedPLC<2, CoordHash> > dummy;
//   // intCells.resize(numGenerators);
//   // for (int i = 0; i < numGenerators; ++i) {
//   //    const IntPoint intGenerator 
//   // }
// }
// //------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template class BoostTessellator<double>;


} //end polytope namespace
