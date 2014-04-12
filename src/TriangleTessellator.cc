//------------------------------------------------------------------------
// TriangleTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <limits>
#include <numeric>

#include "polytope.hh"    // Pulls in POLY_ASSERT and TriangleTessellator.hh.
#include "PLC_CSG_2d.hh"
#include "PLC_CSG_Boost_2d.hh"
#include "simplifyPLCfacets.hh"
#include "polytope_plc_canned_geometries.hh"

// Since triangle isn't built to work out-of-the-box with C++, we 
// slurp in its source here, bracketing it with the necessary dressing.
#define TRILIBRARY
#define REAL double
#define ANSI_DECLARATORS 
#define CDT_ONLY // Conforming Delaunay triangulations only! 

// Because Hang Si has messed with some of Shewchuk's predicates and
// included them with his own Tetgen library, we need to rename some of 
// the symbols therein to prevent duplicate symbols from confusing the 
// linker. Gross.
#define exactinit triangle_exactinit
#define fast_expansion_sum_zeroelim triangle_fast_expansion_sum_zeroelim
#define scale_expansion_zeroelim triangle_scale_expansion_zeroelim
#define estimate triangle_estimate
#define orient3dadapt triangle_orient3dadapt
#define orient3d triangle_orient3d
#define incircleadapt triangle_incircleadapt
#define incircle triangle_incircle
#include "triangle.c"


// Fast predicate for determining colinearity of points.
extern double orient2d(double* pa, double* pb, double* pc);

namespace polytope {

using namespace std;
using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------
// A collection of helper functions
//------------------------------------------------------------------------

namespace {

//------------------------------------------------------------------------------
// Given an array of 3 integers and 2 unique values, find the other one.
//------------------------------------------------------------------------------
void
findOtherTriIndex(const int* indices,
                  const int a,
                  const int b,
		  int& c) {
  POLY_ASSERT(a == indices[0] or a == indices[1] or a == indices[2]);
  POLY_ASSERT(b == indices[0] or b == indices[1] or b == indices[2]);
  POLY_ASSERT(indices[0] != indices[1] and 
              indices[0] != indices[2] and 
              indices[1] != indices[2]);
  if (a != indices[0] and b != indices[0]) {
    c = indices[0];
  }else {
    c = ((a == indices[1] or b == indices[1]) ? indices[2] : indices[1]);
  }
}


//------------------------------------------------------------------------------
// Given an array of 3 integers and 1 unique value, find the other two.
//------------------------------------------------------------------------------
void
findOtherTriIndices(const int* indices,
                    const int a,
                    int& b,
                    int& c) {
  POLY_ASSERT(a == indices[0] or a == indices[1] or a == indices[2]);
  POLY_ASSERT(indices[0] != indices[1] and 
              indices[0] != indices[2] and 
              indices[1] != indices[2]);
  if (a != indices[0]) {
    b = indices[0];
    c = (a != indices[1] ? indices[1] : indices[2]);
  } else {
    b = indices[1];
    c = indices[2];
  } 
}


//------------------------------------------------------------------------------
// Compute the outward-pointing unit vector from the edge of a triangle with
// nodes p1 and p2. pvert is the third vertex of the triangle
//------------------------------------------------------------------------------
template<typename RealType>
void
computeEdgeUnitVector(RealType* p1, 
		      RealType* p2, 
		      RealType* pvert,
		      RealType* result){
  Point2<RealType> test_point, tricent;
  geometry::computeTriangleCentroid2d(p1, p2, pvert, &tricent.x);
  result[0] = -(p2[1] - p1[1]);
  result[1] =  (p2[0] - p1[0]);
  geometry::unitVector<2, RealType>(result);
  copy(p1, p1+2, &test_point.x);
  test_point.x += result[0];
  test_point.y += result[1];
  if ( orient2d(p1, p2, &tricent.x   )*
       orient2d(p1, p2, &test_point.x) > 0.0){
    result[0] *= -1.0;
    result[1] *= -1.0;
  }
}


//------------------------------------------------------------------------------
// Sort a set of edges around a cell so that sequential edges share nodes.
// We allow for one break in the chain (representing on unbounded surface).
// In such a situation we insert the new edge at the beginning of the chain, and
// return "true" indicating that a new edge was created.
// Note that in order for the logic orienting the returned edge ordering to be
// counter-clockwise to work we assume the edges form a convex area!  This 
// should be OK for unbounded tessellations, but don't pass clipped cells to
// this method.
//------------------------------------------------------------------------------
template<typename RealType>
bool
computeSortedCellEdges(const internal::QuantTessellation<2, RealType>& qmesh,
                       std::vector<std::pair<int, int> >& edges,
                       std::vector<int>& result) {

  typedef Point2<RealType> RealPoint;
  typedef std::pair<int, int> EdgeHash;

  unsigned nedges = edges.size();
  POLY_ASSERT(nedges >= 2);

  // Invert the mapping, from nodes to edges.
  std::map<int, std::set<unsigned> > nodes2edges;
  internal::CounterMap<int> nodeUseCount;
  unsigned i;
  for (i = 0; i != nedges; ++i) {
    nodes2edges[edges[i].first].insert(i);
    nodes2edges[edges[i].second].insert(i);
    ++nodeUseCount[edges[i].first];
    ++nodeUseCount[edges[i].second];
  }

  // // BLAGO
  // cerr << "EDGES:" << endl;
  // for (unsigned i = 0; i != nedges; ++i) {
  //   cerr << "  (" << edges[i].first << " " << edges[i].second << ") " << endl;
  // }
  // cerr << "NODES:" << endl;
  // for (typename std::map<int, std::set<unsigned> >::const_iterator itr = nodes2edges.begin();
  //      itr != nodes2edges.end();
  //      ++itr) {
  //   const int j = itr->first;
  //   const std::set<unsigned>& es = itr->second;
  //   cerr << "   " << j << " use count=" << nodeUseCount[j] << " edges=[ ";
  //   std::copy(es.begin(), es.end(), std::ostream_iterator<unsigned>(std::cerr, " "));
  //   std::cerr << "]" << std::endl;
  // }
  // // BLAGO

  // Look for any edges with one node in the set.  There can be at most
  // two such edges, representing the two ends of the chain.  We introduce a
  // new edge hooking those hanging nodes together, and off we go.
  int lastNode;
  vector<int> hangingNodes;
  for (i = 0; i != nedges; ++i) {
    if (nodeUseCount[edges[i].first] == 1 or
        nodeUseCount[edges[i].second] == 1) {
      POLY_ASSERT((nodeUseCount[edges[i].first] == 1 and nodeUseCount[edges[i].second] == 2) or
                  (nodeUseCount[edges[i].first] == 2 and nodeUseCount[edges[i].second] == 1));
      result.push_back(i);
      nodes2edges[edges[i].first].erase(i);
      nodes2edges[edges[i].second].erase(i);
      lastNode = (nodeUseCount[edges[i].first] == 1 ? edges[i].first : edges[i].second);
      hangingNodes.push_back(lastNode);
    }
  }
  POLY_ASSERT(result.size() == 0 or (hangingNodes.size() == 2 and result.size() == 2));

  // If needed create that new edge and put it in the set.
  // Note in 2D we enforce the convention that we will walk counter-clockwise
  // around cells, so we need to pick the correct direction here to start
  // things off.
  if (hangingNodes.size() == 2) {
    result.insert(result.begin() + 1, edges.size());
    edges.push_back(internal::hashEdge(hangingNodes[0], hangingNodes[1]));
    ++nedges;
    // // BLAGO
    // for (unsigned kk = 0; kk != 3; ++kk) {
    //   cerr << "Edge " << result[kk] << " nodes=(" << edges[result[kk]].first << " " << edges[result[kk]].second << ") pos=[" << qmesh.nodePosition(edges[result[kk]].first) << " " << qmesh.nodePosition(edges[result[kk]].second) << endl;
    // }
    // // BLAGO
    const RealPoint n0 = (qmesh.nodePosition(edges[result[0]].first) +
                          qmesh.nodePosition(edges[result[0]].second))*0.5,
                    n1 = (qmesh.nodePosition(edges[result[1]].first) +
                          qmesh.nodePosition(edges[result[1]].second))*0.5,
                    n2 = (qmesh.nodePosition(edges[result[2]].first) +
                          qmesh.nodePosition(edges[result[2]].second))*0.5;
    if (geometry::triangleVolume2(&n0.x, &n1.x, &n2.x) > 0.0) swap(result[0], result[2]);
    // cerr << "Initial vote : " << result[0] << " " << result[1] << " " << result[2] << endl;
    POLY_ASSERT(result.size() == 3);
  }
  POLY_ASSERT(edges.size() == nedges);

  // Pick a node to start the chain.
  if (hangingNodes.size() == 2) {
    POLY_ASSERT(nodeUseCount[edges[result.back()].first] == 2 or
                nodeUseCount[edges[result.back()].second] == 2);
    lastNode = (nodeUseCount[edges[result.back()].first] == 2 ? 
                edges[result.back()].second :
                edges[result.back()].first);
  } else {
    const RealPoint n0 = qmesh.nodePosition(edges[0].first),
                    n1 = qmesh.nodePosition(edges[0].second),
                    n2 = (qmesh.nodePosition(edges[1].first) + 
                          qmesh.nodePosition(edges[1].second))*0.5;
    lastNode = (geometry::triangleVolume2(&n0.x, &n1.x, &n2.x) > 0.0 ?
                edges[0].first :
                edges[1].second);
  }

  // Walk the remaining edges
  EdgeHash ehash;
  while (result.size() != nedges) {
    POLY_ASSERT(nodes2edges[lastNode].size() > 0);
    result.push_back(*nodes2edges[lastNode].begin());
    ehash = edges[result.back()];
    nodes2edges[ehash.first].erase(result.back());
    nodes2edges[ehash.second].erase(result.back());
    lastNode = (ehash.first == lastNode ? ehash.second : ehash.first);
  }
  POLY_ASSERT(result.size() == nedges);
  
  // Set the orientation for the ordered edges.
  for (i = 0; i != nedges; ++i) {
    int j = result[(i + 1) % nedges];
    if (j < 0) j = ~j;
    if ((edges[result[i]].first == edges[j].first) or 
        (edges[result[i]].first == edges[j].second)) result[i] = ~result[i];
  }

  // // BLAGO
  // cerr << "Raw results: ";
  // for (unsigned i = 0; i != result.size(); ++i) {
  //   int j = result[i];
  //   if (j < 0) j = ~j;
  //   cerr << " [" << result[i] << " (" << edges[j].first << " " << edges[j].second << ")]";
  // }
  // cerr << endl << "Sorted edges: ";
  // for (unsigned i = 0; i != result.size(); ++i) {
  //   int j = result[i];
  //   if (j < 0) {
  //     j = ~j;
  //     cerr << " (" << edges[j].second << " " << edges[j].first << ")";
  //   } else {
  //     cerr << " (" << edges[j].first << " " << edges[j].second << ")";
  //   }
  // }
  // cerr << endl;
  // // BLAGO

  // That's it.
  POLY_BEGIN_CONTRACT_SCOPE;
  {
    POLY_ASSERT(edges.size() == result.size());
    for (int i = 0; i != edges.size(); ++i) {
      const int j = (i + 1) % edges.size();
      const int ii = result[i];
      const int jj = result[j];
      POLY_ASSERT((ii >= 0 ? ii : ~ii) < edges.size());
      POLY_ASSERT((jj >= 0 ? jj : ~jj) < edges.size());
      POLY_ASSERT(((ii >= 0 and jj >= 0) and edges[ii].second == edges[jj].first) or
                  ((ii >= 0 and jj <  0) and edges[ii].second == edges[~jj].second) or
                  ((ii <  0 and jj >= 0) and edges[~ii].first == edges[jj].first) or
                  ((ii <  0 and jj <  0) and edges[~ii].first == edges[~jj].second));
    }
  }
  POLY_END_CONTRACT_SCOPE;
  return !(hangingNodes.empty());
}

//------------------------------------------------------------------------------
// Build a ReducedPLC representation of a 2D cell.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<2, int64_t>
hashReducedPLC(const ReducedPLC<2, RealType>& plc,
	       const RealType* xlow_inner,
	       const RealType* xhigh_inner,
	       const RealType* xlow_outer,
	       const RealType* xhigh_outer,
	       const RealType minTol) {
  typedef geometry::Hasher<2, RealType> HasherType;
  ReducedPLC<2, int64_t> result;
  result.facets = plc.facets;
  result.points.resize(plc.points.size());
  for (unsigned i = 0; i != plc.points.size()/2; ++i) {
    const int64_t ip = HasherType::hashPosition(&plc.points[2*i],
						xlow_inner,
						xhigh_inner,
						xlow_outer,
						xhigh_outer,
						minTol);
    result.points[2*i  ] = HasherType::qxval(ip);
    result.points[2*i+1] = HasherType::qyval(ip);
  }
  return result;
}


//------------------------------------------------------------------------------
// Build a ReducedPLC representation of a 2D cell.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<2, RealType>
unhashReducedPLC(const ReducedPLC<2, int64_t>& plc,
		 const RealType* xlow_inner,
		 const RealType* xhigh_inner,
		 const RealType* xlow_outer,
		 const RealType* xhigh_outer,
		 const RealType minTol) {
  typedef geometry::Hasher<2, RealType> HasherType;
  ReducedPLC<2, RealType> result;
  result.facets = plc.facets;
  result.points.resize(plc.points.size());
  for (unsigned i = 0; i != plc.points.size()/2; ++i) {
    RealType pos[2];
    const uint64_t hashedPosition = plc.points[2*i] + (plc.points[2*i+1] << 31);
    HasherType::unhashPosition(pos,
			       xlow_inner,
			       xhigh_inner,
			       xlow_outer,
			       xhigh_outer,
			       hashedPosition,
			       minTol);
    result.points[2*i  ] = pos[0];
    result.points[2*i+1] = pos[1];
  }
  return result;
}


//------------------------------------------------------------------------------
// Build a ReducedPLC representation of a 2D cell.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<2, RealType>
plcOfCell(const internal::QuantTessellation<2, RealType>& qmesh,
	  const unsigned icell) {
  POLY_ASSERT(icell < qmesh.cells.size());
  typedef Point2<RealType> PointType;
  ReducedPLC<2, RealType> result;
  const unsigned nFaces = qmesh.cells[icell].size();
  result.facets.resize(nFaces, vector<int>(2));
  for (unsigned i = 0; i != nFaces; ++i) {
    const bool flip = qmesh.cells[icell][i] < 0;
    const unsigned iface = flip ? ~qmesh.cells[icell][i] : qmesh.cells[icell][i];
    POLY_ASSERT(iface < qmesh.faces.size());
    POLY_ASSERT(qmesh.faces[iface].size() == 1);
    POLY_ASSERT(qmesh.faces[iface][0] >= 0);
    const unsigned iedge = qmesh.faces[iface][0];
    POLY_ASSERT(iedge < qmesh.edges.size());
    const unsigned ip = flip ? qmesh.edges[iedge].second : qmesh.edges[iedge].first;
    cerr << "Edge : " << iedge << " " << flip << " : (" << qmesh.edges[iedge].first << " " << qmesh.edges[iedge].second << ") : "
         << qmesh.unhashPosition(qmesh.points[qmesh.edges[iedge].first]) << " "
         << qmesh.unhashPosition(qmesh.points[qmesh.edges[iedge].second]) << endl;
    POLY_ASSERT(ip < qmesh.points.size());
    PointType p = qmesh.unhashPosition(qmesh.points[ip]);
    result.points.push_back(p.x);
    result.points.push_back(p.y);
    result.facets[i][0] = i;
    result.facets[i][1] = (i+1)%nFaces;
  }
  POLY_ASSERT(result.points.size()/2 == nFaces);
  return result;
}



} // end anonymous namespace



//------------------------------------------------------------------------------
template<typename RealType>
TriangleTessellator<RealType>::
TriangleTessellator():
  Tessellator<2, RealType>(){
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
TriangleTessellator<RealType>::
~TriangleTessellator() {
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Unbounded tessellation
//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           Tessellation<2, RealType>& mesh) const {
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 2 == 0);
  POLY_ASSERT(points.size() > 2 );
  
  // Generate the internal quantized tessellation.
  internal::QuantTessellation<2, RealType> qmesh;
  vector<RealType> nonGeneratingPoints;
  this->computeUnboundedQuantizedTessellation(points, nonGeneratingPoints, qmesh);

//   // Blago!
//   cerr << "QuantMesh edges:" << endl;
//   for (unsigned i = 0; i != qmesh.edges.size(); ++i) {
//     cerr << "  Edge " << i << ": "
// 	 << qmesh.edges[i].first << ","
// 	 << qmesh.edges[i].second << endl;
//   }
//   cerr << "QuantMesh faces:" << endl;
//   for (unsigned i = 0; i != qmesh.faces.size(); ++i) {
//     cerr << "  Face " << i << ":" << endl << "    ";
//     for (unsigned j = 0; j != qmesh.faces[i].size(); ++j) {
//       if (qmesh.faces[i][j] >= 0)
// 	cerr << qmesh.edges[qmesh.faces[i][j]].first << " ";
//       else
// 	cerr << qmesh.edges[~qmesh.faces[i][j]].second << " ";
//     }
//     cerr << endl;
//   }
//   // Blago!

  // Convert to output tessellation.
  qmesh.tessellation(mesh);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Tessellate within a box
//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
	   RealType* low,
	   RealType* high,
           Tessellation<2, RealType>& mesh) const {
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 2 == 0);
  POLY_ASSERT(points.size() > 2 );
  POLY_ASSERT(low[0] <= high[0] and low[1] <= high[1]);

  // Create a reduced PLC of the bounding box and use the reduced PLC method
  ReducedPLC<2, RealType> box = plc_box<2,RealType>(low, high);
  this->tessellate(points, box, mesh);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Tessellate within a PLC
//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
	   const vector<RealType>& PLCpoints,
	   const PLC<2, RealType>& geometry,
           Tessellation<2, RealType>& mesh) const {
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(not points.empty() and not PLCpoints.empty());
  POLY_ASSERT(points.size() % 2 == 0 and PLCpoints.size() % 2 == 0);
  POLY_ASSERT(points.size() > 2 );
  POLY_ASSERT(not geometry.empty());
  
  // Export to the ReducedPLC method.
  ReducedPLC<2, RealType> boundary;
  boundary.facets = geometry.facets;
  boundary.holes  = geometry.holes;
  boundary.points = PLCpoints;
  this->tessellate(points, boundary, mesh);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Tessellate within a ReducedPLC
//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
	   const ReducedPLC<2, RealType>& geometry,
           Tessellation<2, RealType>& mesh) const {
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 2 == 0);
  POLY_ASSERT(points.size() > 2 );
  POLY_ASSERT(not geometry.empty());

  
  typedef geometry::Hasher<2, RealType> HasherType;
  typedef typename internal::QuantTessellation<2, RealType>::PointHash PointHash;
  typedef typename internal::QuantTessellation<2, RealType>::EdgeHash EdgeHash;
  typedef typename internal::QuantTessellation<2, RealType>::IntPoint IntPoint;
  typedef typename internal::QuantTessellation<2, RealType>::RealPoint RealPoint;

  const unsigned numGenerators = points.size()/2;
  internal::QuantTessellation<2,RealType> qmesh0;
  this->computeUnboundedQuantizedTessellation(points, geometry.points, qmesh0);
  std::cerr << "Computed unbounded quantized tessellation" << qmesh0 << std::endl;

  // Create a new QuantTessellation.  This one will only use the single level of
  // quantization since we know the PLC is within this inner region.
  internal::QuantTessellation<2, RealType> qmesh1;
  qmesh1.generators    = qmesh0.generators;
  qmesh1.low_labframe  = qmesh0.low_labframe;
  qmesh1.high_labframe = qmesh0.high_labframe;
  qmesh1.low_inner     = qmesh0.low_inner;
  qmesh1.high_inner    = qmesh0.high_inner;
  qmesh1.low_outer     = qmesh0.low_inner;
  qmesh1.high_outer    = qmesh0.high_inner;
  qmesh1.degeneracy    = qmesh0.degeneracy;

#if HAVE_BOOST
  ReducedPLC<2, RealType> normalizedGeometry;
  normalizedGeometry.facets = geometry.facets;
  normalizedGeometry.holes  = geometry.holes;
  normalizedGeometry.points = this->computeNormalizedPoints(geometry.points,
							    geometry.points,
							    false,
							    &qmesh0.low_labframe.x,
							    &qmesh0.high_labframe.x);
  const ReducedPLC<2, CoordHash> IntGeometry = hashReducedPLC(normalizedGeometry,
							      const_cast<RealType*>(&qmesh0.low_inner.x),
							      const_cast<RealType*>(&qmesh0.high_inner.x),
							      const_cast<RealType*>(&qmesh0.low_outer.x),
							      const_cast<RealType*>(&qmesh0.high_outer.x),
							      qmesh0.degeneracy);
  std::cerr << "Normalized Geometry:\n" << normalizedGeometry << std::endl;
  std::cerr << "Hashed Geometry:\n" << IntGeometry << std::endl;
#endif

  // Walk the cells in the unbounded tessellation
  for (unsigned icell = 0; icell != numGenerators; ++icell) {
    
    // Intersect cell with boundary
    //
    // Do the clipping integers if Boost.Geometry is available. 
    // Otherwise, reduce to using CSG in floating point.

    // Build a ReducedPLC to represent the cell
    ReducedPLC<2, RealType> cell = plcOfCell(qmesh0, icell);

    std::cerr << "\n------------------------------ Clipping cell " << icell << std::endl;
    std::cerr << "  \nPre-clipped cell:\n" << cell << std::endl;
    
#if HAVE_BOOST
    ReducedPLC<2, CoordHash> IntCell = hashReducedPLC(cell,
						      const_cast<RealType*>(&qmesh0.low_inner.x),
						      const_cast<RealType*>(&qmesh0.high_inner.x),
						      const_cast<RealType*>(&qmesh0.low_outer.x),
						      const_cast<RealType*>(&qmesh0.high_outer.x),
						      qmesh0.degeneracy);
    const CoordHash pid = HasherType::hashPosition(&qmesh0.generators[2*icell],
						   const_cast<RealType*>(&qmesh0.low_inner.x),
						   const_cast<RealType*>(&qmesh0.high_inner.x),
						   const_cast<RealType*>(&qmesh0.low_outer.x),
						   const_cast<RealType*>(&qmesh0.high_outer.x),
						   qmesh0.degeneracy);
    std::cerr << "   Generator: (" << qmesh1.generators[2*icell] << "," << qmesh1.generators[2*icell+1] << ")"
	      << " --> " << pid << " --> " << IntPoint(HasherType::qxval(pid), HasherType::qyval(pid)) << std::endl;
    std::vector<ReducedPLC<2, CoordHash> > orphans;

    std::cerr << "  \nPre-clipped IntCell:\n" << IntCell << std::endl;

    std::cerr << "   Clip... ";
    IntCell = BG::boost_clip(IntGeometry,
			     IntCell,
			     IntPoint(HasherType::qxval(pid), HasherType::qyval(pid)),
			     orphans);
    std::cerr << "Resulting PLC: " << IntCell << std::endl;
    std::cerr << "DONE!" << std::endl;

    // Blago!
    if (not orphans.empty())  cerr << "Orphans detected, but no actions taken" << endl;
    // Blago!

    cell = unhashReducedPLC(IntCell,
			    const_cast<RealType*>(&qmesh0.low_inner.x),
			    const_cast<RealType*>(&qmesh0.high_inner.x),
			    const_cast<RealType*>(&qmesh0.low_outer.x),
			    const_cast<RealType*>(&qmesh0.high_outer.x),
			    qmesh1.degeneracy);

//     std::cerr << "  \nPost-clipped IntCell:\n" << IntCell << std::endl;
//     std::cerr << "  \nPost-clipped cell:\n" << cell << std::endl;


#else
    cell = CSG::csg_intersect(geometry, cell);
    cell = simplifyPLCfacets(cell, 
			     cell.points,
			     &qmesh1.low_inner.x,
			     &qmesh1.high_inner.x,
			     1.0e-5);
#endif
    POLY_ASSERT(cell.facets.size() >= 3);

    // Add cell and its elements to the new tessellation
    vector<int> nodeIDs, edgeIDs, faceIDs;
    qmesh1.cells.push_back(vector<int>());
    for (unsigned i = 0; i != cell.points.size()/2; ++i) {
      nodeIDs.push_back(qmesh1.addNewNode(HasherType::hashPosition(&cell.points[2*i],
								   const_cast<RealType*>(&qmesh1.low_inner.x),
								   const_cast<RealType*>(&qmesh1.high_inner.x),
								   const_cast<RealType*>(&qmesh1.low_outer.x),
								   const_cast<RealType*>(&qmesh1.high_outer.x),
								   qmesh1.degeneracy)));
    }
    for (unsigned iface = 0; iface != cell.facets.size(); ++iface) {
      const unsigned nnodes = cell.facets[iface].size();
      POLY_ASSERT(nnodes == 2);
      vector<int> face;
      for (unsigned i = 0; i != nnodes; ++i) {
	const unsigned j = (i+1) % nnodes;
	const EdgeHash ehash = internal::hashEdge(nodeIDs[cell.facets[iface][i]],
						  nodeIDs[cell.facets[iface][j]]);
	face.push_back(qmesh1.addNewEdge(ehash));
	if (ehash.first == nodeIDs[cell.facets[iface][j]]) face.back() = ~face.back();
      }
      POLY_ASSERT(face.size() == nnodes);
      const unsigned k = qmesh1.faces.size();
      const unsigned i = qmesh1.addNewFace(face);
      qmesh1.cells.back().push_back(i == k ? i : ~i);
    }
    POLY_ASSERT(qmesh1.cells.back().size() == cell.facets.size());
  }

  // Check the validity of the quantized tessellation
  qmesh1.assertValid();

  // Convert to output tessellation
  qmesh1.tessellation(mesh);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Internal method that returns an intermediated quantized representation
// of the unbounded tessellation.
//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeUnboundedQuantizedTessellation(const vector<RealType>& points,
                                      const vector<RealType>& nonGeneratingPoints,
                                      internal::QuantTessellation<2, RealType>& qmesh) const {

  typedef typename internal::QuantTessellation<2, RealType>::PointHash PointHash;
  typedef typename internal::QuantTessellation<2, RealType>::EdgeHash EdgeHash;
  typedef Point2<RealType> RealPoint;

  qmesh.degeneracy = mDegeneracy;

  // Check for collinearity and use the appropriate routine
  bool isCollinear = geometry::collinear<2,RealType>(points, 1.0e-10);

  // Call a special routine to build up the quantized tessellation if the
  // input points are really 1D. This routine is purely geometric and is
  // independent of tessellator.
  if (isCollinear) {
    POLY_ASSERT2(false, "Need to write the special 1D tessellation routine "
		 << "for collinear 2D points for the QuantTessellation.");
    ////computeUnboundedQuantizedTessellationCollinear(points, nonGeneratingPoints, qmesh);
  } 

  // It's a fully-2D problem. Do the tessellator-specific stuff
  else {
    
    // Normalize the input generators
    const unsigned numGenerators = points.size()/2;
    unsigned i, j, k;
//     qmesh.generators = points;
//     geometry::computeBoundingBox<2,RealType>(points, 
// 					     true, 
// 					     &qmesh.low_labframe.x, 
// 					     &qmesh.high_labframe.x);
//     qmesh.low_inner = qmesh.low_labframe;
//     qmesh.high_inner = qmesh.high_labframe;
    qmesh.generators = this->computeNormalizedPoints(points, 
						     nonGeneratingPoints, 
						     true,
						     &qmesh.low_labframe.x, 
						     &qmesh.high_labframe.x);


//     for (i=0; i<numGenerators; ++i)
//       cerr << i << ": "
// 	   << "(" << points[2*i] << "," << points[2*i+1] << ") --> "
// 	   << "(" << qmesh.generators[2*i] << "," << qmesh.generators[2*i+1] << ")   "
// 	   << sqrt((qmesh.generators[2*i  ]-0.5)*(qmesh.generators[2*i  ]-0.5) +
// 		   (qmesh.generators[2*i+1]-0.5)*(qmesh.generators[2*i+1]-0.5)) << endl;
    
    // Compute the triangularization
    triangulateio delaunay;
    computeDelaunay(qmesh.generators, delaunay);
    
    // Find the circumcenters of each triangle, and build the set of triangles
    // associated with each generator.
    qmesh.low_inner  = RealPoint(0.0, 0.0);
    qmesh.high_inner = RealPoint(1.0, 1.0);
    qmesh.low_outer  = RealPoint( numeric_limits<RealType>::max(),
				  numeric_limits<RealType>::max());
    qmesh.high_outer = RealPoint(-numeric_limits<RealType>::max(),
				 -numeric_limits<RealType>::max());
    vector<RealPoint> circumcenters(delaunay.numberoftriangles);
    int p, q, r;
    EdgeHash pq, pr, qr;
    vector<unsigned> triMask(delaunay.numberoftriangles, 0);
    map<EdgeHash, vector<unsigned> > edge2tris;
    map<int, set<unsigned> > gen2tri;
    for (i = 0; i != delaunay.numberoftriangles; ++i) {
      p  = delaunay.trianglelist[3*i  ];
      q  = delaunay.trianglelist[3*i+1];
      r  = delaunay.trianglelist[3*i+2];
//       POLY_ASSERT(p < numGenerators and q < numGenerators and r < numGenerators);
      geometry::computeCircumcenter(&delaunay.pointlist[2*p],
				    &delaunay.pointlist[2*q],
				    &delaunay.pointlist[2*r],
				    &circumcenters[i].x);
      pq = internal::hashEdge(p,q);
      pr = internal::hashEdge(p,r);
      qr = internal::hashEdge(q,r);
      POLY_ASSERT(orient2d(&delaunay.pointlist[2*p],
			   &delaunay.pointlist[2*q],
			   &delaunay.pointlist[2*r]) != 0);
      if (p < numGenerators and q < numGenerators and r < numGenerators) {
	triMask[i] = 1;
	edge2tris[pq].push_back(i);
	edge2tris[pr].push_back(i);
	edge2tris[qr].push_back(i);
	gen2tri[p].insert(i);
	gen2tri[q].insert(i);
	gen2tri[r].insert(i);
	qmesh.low_outer.x  = min(qmesh.low_outer.x , circumcenters[i].x);
	qmesh.low_outer.y  = min(qmesh.low_outer.y , circumcenters[i].y);
	qmesh.high_outer.x = max(qmesh.high_outer.x, circumcenters[i].x);
	qmesh.high_outer.y = max(qmesh.high_outer.y, circumcenters[i].y);
      }
    }
    POLY_ASSERT(circumcenters.size() == delaunay.numberoftriangles);
    POLY_ASSERT(gen2tri.size() == numGenerators);
    POLY_ASSERT(std::accumulate(triMask.begin(), triMask.end(), 0) > 0);
    
    POLY_BEGIN_CONTRACT_SCOPE;
    {
      for (typename map<EdgeHash, vector<unsigned> >::const_iterator itr = edge2tris.begin();
	   itr != edge2tris.end();
	   ++itr) POLY_ASSERT(itr->second.size() == 1 or itr->second.size() == 2);
      for (map<int, set<unsigned> >::const_iterator itr = gen2tri.begin();
	   itr != gen2tri.end();
	   ++itr) POLY_ASSERT(itr->second.size() >= 1);
      POLY_ASSERT(qmesh.low_outer.x  <= qmesh.high_outer.x);
      POLY_ASSERT(qmesh.low_outer.y  <= qmesh.high_outer.y);
      POLY_ASSERT(qmesh.low_inner.x  >= qmesh.low_outer.x );
      POLY_ASSERT(qmesh.low_inner.y  >= qmesh.low_outer.y );
      POLY_ASSERT(qmesh.high_inner.x <= qmesh.high_outer.x);
      POLY_ASSERT(qmesh.high_inner.y <= qmesh.high_outer.y);
    }
    POLY_END_CONTRACT_SCOPE;
  
    // Expand the outer bounding box and choose infinite sphere radius
    qmesh.low_outer.x  = min(qmesh.low_outer.x , qmesh.low_inner.x );
    qmesh.low_outer.y  = min(qmesh.low_outer.y , qmesh.low_inner.y );
    qmesh.high_outer.x = max(qmesh.high_outer.x, qmesh.high_inner.x);
    qmesh.high_outer.y = max(qmesh.high_outer.y, qmesh.high_inner.y);
    RealType rinf = 4.0*max(qmesh.high_outer.x - qmesh.low_outer.x,
			    qmesh.high_outer.y - qmesh.low_outer.y);
    const RealPoint centroid_outer = (qmesh.low_outer + qmesh.high_outer)/2;
    qmesh.low_outer.x  = centroid_outer.x - 1.05*rinf;
    qmesh.low_outer.y  = centroid_outer.y - 1.05*rinf;
    qmesh.high_outer.x = centroid_outer.x + 1.05*rinf;
    qmesh.high_outer.y = centroid_outer.y + 1.05*rinf;

    // Quantize circumcenters and add map them to unique IDs
    map<int, unsigned> tri2id;
    for (i = 0; i != delaunay.numberoftriangles; ++i) {
      if (triMask[i] == 1) {
	tri2id[i] = qmesh.addNewNode(circumcenters[i]);
      }
    }
    POLY_ASSERT(tri2id.size() == std::accumulate(triMask.begin(), triMask.end(), 0));
    //POLY_ASSERT(tri2id.size() == delaunay.numberoftriangles);
  
    // The exterior edges of the triangularization have "unbounded" rays, originating
    // at the circumcenter of the corresponding triangle and passing perpendicular to
    // the edge. Find those surface edges and project unbounded rays through them.
    bool test;
    RealPoint ehat, pinf, tricent;
    map<EdgeHash, unsigned> projEdge2id;
    int i1, i2, ivert;
    qmesh.infNodes = vector<unsigned>();
    for (typename map<EdgeHash, vector<unsigned> >::const_iterator edgeItr = edge2tris.begin();
         edgeItr != edge2tris.end();
         ++edgeItr) {
      const EdgeHash& edge = edgeItr->first;
      const vector<unsigned>& tris = edgeItr->second;
      if (tris.size() == 1) {
        i = tris[0];
        POLY_ASSERT(i < delaunay.numberoftriangles);
        i1 = edge.first;
        i2 = edge.second;
        findOtherTriIndex(&delaunay.trianglelist[3*i], i1, i2, ivert);
        computeEdgeUnitVector(&delaunay.pointlist[2*i1],
			      &delaunay.pointlist[2*i2],
			      &delaunay.pointlist[2*ivert],
			      &ehat.x);
        
        // Compute the intersection of the infinite edge with the inf sphere
        test = geometry::rayCircleIntersection(&circumcenters[i].x,
					       &ehat.x,
					       &centroid_outer.x,
					       rinf,
					       1.0e-10,
					       &pinf.x);
        POLY_ASSERT(test);
  
        // Add the projected point to the quantized tessellation
        k = qmesh.point2id.size();
        j = qmesh.addNewNode(pinf);
        POLY_ASSERT(projEdge2id.find(edge) == projEdge2id.end());
        projEdge2id[edge] = j;
        if (k != qmesh.point2id.size()) qmesh.infNodes.push_back(j);
      }
    }
  
    
//     // Blago!
//     cerr << "Nodes and IDs:" << endl;
//     for (int pp=0; pp<qmesh.points.size(); ++pp) cerr << pp << "  " << qmesh.nodePosition(pp) << endl;
//     // Blago!


    // The faces corresponding to each triangle edge
    qmesh.faces.reserve(edge2tris.size());
    qmesh.cells = vector<vector<int> >(numGenerators);
    int iedge, iface;
    unsigned ii, jj;
    RealType vol;
    RealPoint n0, n1;
    for (map<int, set<unsigned> >::const_iterator genItr = gen2tri.begin();
         genItr != gen2tri.end(); 
         ++genItr) {
      p = genItr->first;
//       std::cerr << "Generator " << p << " at " << points[2*p] << "," << points[2*p+1] << std::endl;
      const set<unsigned>& tris = genItr->second;
      POLY_ASSERT(p < numGenerators);
      vector<EdgeHash> meshEdges;
      for (set<unsigned>::const_iterator triItr = tris.begin();
           triItr != tris.end(); 
  	 ++triItr){
        i = *triItr;
        POLY_ASSERT(i < delaunay.numberoftriangles);
        POLY_ASSERT(tri2id.find(i) != tri2id.end());
        ii = tri2id[i];
        
// 	std::cerr << "  " << ii << " " << qmesh.nodePosition(ii) << std::endl;

        // Get the other indices for this triangle, given one of its vertices p
        findOtherTriIndices(&delaunay.trianglelist[3*i], p, q, r);
        pq = internal::hashEdge(p,q);
        pr = internal::hashEdge(p,r);
        
        // Is pq a surface edge?
        if (edge2tris[pq].size() == 1){
          POLY_ASSERT(edge2tris[pq][0] == i);
          POLY_ASSERT(projEdge2id.find(pq) != projEdge2id.end());
          jj = projEdge2id[pq];
	  POLY_ASSERT(jj != ii);
	  meshEdges.push_back(internal::hashEdge(ii,jj));
        } else {
	  POLY_ASSERT((edge2tris[pq].size() == 2 and edge2tris[pq][0] == i)
		      or edge2tris[pq][1] == i);
          k = (edge2tris[pq][0] == i ? edge2tris[pq][1] : edge2tris[pq][0]);
	  POLY_ASSERT(tri2id.find(k) != tri2id.end());
          jj = tri2id[k];
	  if (jj != ii) meshEdges.push_back(internal::hashEdge(ii,jj));
        }
        
        // Is pr a surface edge?
        if (edge2tris[pr].size() == 1){
          POLY_ASSERT(edge2tris[pr][0] == i);
          POLY_ASSERT(projEdge2id.find(pr) != projEdge2id.end());
          jj = projEdge2id[pr];
	  POLY_ASSERT(ii != jj);
          meshEdges.push_back(internal::hashEdge(ii,jj));
        } else {
	  POLY_ASSERT((edge2tris[pr].size() == 2 and edge2tris[pr][0] == i)
		      or edge2tris[pr][1] == i);
          k = (edge2tris[pr][0] == i ? edge2tris[pr][1] : edge2tris[pr][0]);
	  POLY_ASSERT(tri2id.find(k) != tri2id.end());
          jj = tri2id[k];
	  if (jj != ii) meshEdges.push_back(internal::hashEdge(ii,jj));
        }
      }
  
      // Arrange the edges in the corectly sorted and sign oriented order
      sort(meshEdges.begin(), meshEdges.end());
      meshEdges.erase(unique(meshEdges.begin(), meshEdges.end()), meshEdges.end());
      if (meshEdges.size() > 1) {
        vector<int> edgeOrder;
        const bool infEdge = computeSortedCellEdges(qmesh, meshEdges, edgeOrder);
	
// 	// Blago!
// 	cerr << "    ";
// 	for (int pp=0; pp<meshEdges.size(); ++pp)
// 	  cerr << "(" << meshEdges[pp].first << "," << meshEdges[pp].second << ") ";
// 	cerr << endl;
// 	cerr << "    ";
// 	for (int pp=0; pp<edgeOrder.size(); ++pp)
// 	  cerr << edgeOrder[pp] << " ";
// 	cerr << endl;
// 	// Blago!

        if (meshEdges.size() > 2) {
  	
          // Add the edges and faces to the quantized mesh. (They are equal in 2D.)
          for (vector<int>::const_iterator itr = edgeOrder.begin();
               itr != edgeOrder.end();
               ++itr) {
            const bool flip = (*itr < 0);
            k = (flip ? ~(*itr) : *itr);
            iedge = qmesh.addNewEdge(meshEdges[k]);
            vector<int> face(1, iedge);
            iface = qmesh.addNewFace(face);
            POLY_ASSERT(iedge == iface);
  
            // Determine the orientation of the face with respect to the cell
            n0 = qmesh.nodePosition(meshEdges[k].first);
            n1 = qmesh.nodePosition(meshEdges[k].second);
            vol = geometry::triangleVolume2(&qmesh.generators[2*p], &n1.x, &n0.x);
            POLY_ASSERT(vol != 0.0);
            if (vol > 0.0) {
              qmesh.cells[p].push_back(iface);
            } else {
              qmesh.cells[p].push_back(~iface);
            }
          }
  	
          // Did we create a new infEdge? If so we know it was the second element
          // in the ordered list.
          if (infEdge) {
            j = internal::positiveID(edgeOrder[1]);
            k = qmesh.edge2id[meshEdges[j]];
            qmesh.infEdges.push_back(k);
            qmesh.infFaces.push_back(k);
          }
        }
      }
      
      // How does meshEdges only have one element?
      else {
        cerr << "BLAGO!" << endl 
	     << p << " " << tris.size() << " " << meshEdges.size() << endl;
      }
      POLY_ASSERT(meshEdges.size() > 1);
    }
    
    // All infFaces have been stored in the quantized mesh at this point.
    // Two complications may still exist for an infinite cell:
    //
    // 1. Projected edges may intersect
    // This can occur when a generator on the boundary has two surface triangles
    // that are nearly flat, but the two triangle edges on the surface are not
    // collinear. There is a critical threshold in which Triangle does
    // recognize the edges as collinear, but projecting rays orthogonal to the 
    // edges creates an intersection if the inf sphere is sufficiently large.
    // The internal floating point precision of Triangle creates this error. If
    // Triangle could recognize the edges as not being exactly collinear, it would
    // add a third triangle with circumcenter at the position of the ray intersection.
    // In this instance, the two original triangles are no longer on the surface,
    // and the generator is now internal.
    //
    // 2. Infinite face may intersect domain again
    // If projected edges are nearly collinear, then the inf face connecting their
    // projected nodes could intersect the internal bounding box (or PLC boundary,
    // if it exists). An additional node will have to be projected in this instance
    // in between the previous two. Two infinite faces will be constructed for this
    // cell in this case. The existing infEdge and infFace data in qmesh will need
    // to be modified.
    //
    // This is where to do it at some point...
  
////   for (i = 0; i != numGenerators; ++i) {
////     POLY_ASSERT(cellInfEdges[i].size() == 0 or cellInfEdges[i].size() == 1);
////     if (cellInfEdges[i].size() == 1) {
////     }
////   }
  
    // Post-conditions
    qmesh.assertValid();
  
    // Clean up.
    trifree((VOID*)delaunay.pointlist);
    trifree((VOID*)delaunay.pointmarkerlist);
    trifree((VOID*)delaunay.trianglelist);
    trifree((VOID*)delaunay.edgelist);
    trifree((VOID*)delaunay.edgemarkerlist);
    trifree((VOID*)delaunay.segmentlist);
    trifree((VOID*)delaunay.segmentmarkerlist);
  }
}  

//------------------------------------------------------------------------------
// Private routines called by tessellate:
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeDelaunay(const vector<RealType>& points,
                triangulateio& delaunay) const {
  triangulateio in;
   
  // Find the range of the generator points.
  const unsigned numGenerators = points.size()/2;

//   POLY_ASSERT(!mCoords.empty());
//   RealType low [2] = {mCoords.low [0], mCoords.low [1]};
//   RealType high[2] = {mCoords.high[0], mCoords.high[1]};

  // Determine bounding box for points
  RealType low[2], high[2];
  geometry::computeBoundingBox<2,RealType>(points, true, low, high);
//   POLY_ASSERT(low[0] >= 0.0 and high[0] <= 1.0);
//   POLY_ASSERT(low[1] >= 0.0 and high[1] <= 1.0);
  
  RealType box [2] = {high[0] - low[0], high[1] - low[1]};
  const RealType boxsize = 8.0*max(box[0], box[1]);
  
  const RealType xmin = 0.5*(low[0] + high[0]) - boxsize;
  const RealType xmax = 0.5*(low[0] + high[0]) + boxsize;
  const RealType ymin = 0.5*(low[1] + high[1]) - boxsize;
  const RealType ymax = 0.5*(low[1] + high[1]) + boxsize;

  // Add the generators
  in.numberofpoints = numGenerators + 4;
  in.pointlist = new RealType[2*in.numberofpoints];
  copy(points.begin(), points.end(), in.pointlist);
  in.pointlist[2*numGenerators  ] = xmin;  in.pointlist[2*numGenerators+1] = ymin;
  in.pointlist[2*numGenerators+2] = xmax;  in.pointlist[2*numGenerators+3] = ymin;
  in.pointlist[2*numGenerators+4] = xmax;  in.pointlist[2*numGenerators+5] = ymax;
  in.pointlist[2*numGenerators+6] = xmin;  in.pointlist[2*numGenerators+7] = ymax;
  in.numberofsegments = 0;

  // // Add the generators
  // in.numberofpoints = numGenerators;
  // in.pointlist = new RealType[2*in.numberofpoints];
  // copy(points.begin(), points.end(), in.pointlist);
  // in.numberofsegments = 0;
    
  // No point attributes or markers.
  in.numberofpointattributes = 0;
  in.pointattributelist = 0; 
  in.pointmarkerlist = 0;   // No point markers.
  in.segmentmarkerlist = 0;
  in.numberofholes = 0;
  in.holelist = 0;

  // No regions.
  in.numberofregions = 0;
  in.regionlist = 0;
  
  // Set up the structure for the triangulation.
  delaunay.pointlist = 0;
  delaunay.pointattributelist = 0;
  delaunay.pointmarkerlist = 0;
  delaunay.trianglelist = 0;
  delaunay.triangleattributelist = 0;
  delaunay.neighborlist = 0;
  delaunay.segmentlist = 0;
  delaunay.segmentmarkerlist = 0;
  delaunay.edgelist = 0;
  delaunay.edgemarkerlist = 0;
  delaunay.holelist = 0;

  // Do the triangulation. Switches pass to triangle are:
  // -Q : Quiet (shaddap!), no output on the terminal except errors.
  // -z : Indices are all numbered from zero.
  // -e : Generates edges and places them in out.edgelist.
  // -c : Generates convex hull and places it in out.segmentlist.
  // -p : Uses the given PLC information.
  // triangulate((char*)"Qzec", &in, &delaunay, 0);
  triangulate((char*)"Qz", &in, &delaunay, 0);
  
  // Make sure we got something.
  if (delaunay.numberoftriangles == 0)
    error("TriangleTessellator: Delauney triangulation produced 0 triangles!");
  if (delaunay.numberofpoints != in.numberofpoints) {
    char err[1024];
    snprintf(err, 1024, "TriangleTessellator: Delauney triangulation produced %d triangles\n(%d generating points given)", 
             delaunay.numberofpoints, (int)numGenerators);
    error(err);
  }
  
  // Clean up
  delete [] in.pointlist;
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TriangleTessellator<double>;

}
