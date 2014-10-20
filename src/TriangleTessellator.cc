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
#include "PLC_Boost_2d.hh"
#include "simplifyPLCfacets.hh"
#include "polytope_plc_canned_geometries.hh"

// Since triangle isn't built to work out-of-the-box with C++, we 
// slurp in its source here, bracketing it with the necessary dressing.
#define TRILIBRARY
#define REAL double
#define ANSI_DECLARATORS 
#define CDT_ONLY // Conforming Delaunay triangulations only! 

#define NORMALIZE_GENERATORS false
#define DUMP_UNBOUNDED_DEBUG_MESH true

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
// Compute the outward-pointing unit vector from an edge having vertices p1, p2
// NOTE: This assumes the edges and vertices around a figure are ordered CCW
//------------------------------------------------------------------------------
template<typename RealType>
void
computeEdgeNormal(const RealType* p1, 
		  const RealType* p2,
		  RealType* result) {
  result[0] =  (p2[1] - p1[1]);
  result[1] = -(p2[0] - p1[0]);
  geometry::unitVector<2, RealType>(result);
}


//------------------------------------------------------------------------------
// Sort a set of edges around a face so that sequential edges share nodes.
// We allow for one break in the chain (representing on unbounded surface).
// In such a situation we insert the new edge at the beginning of the chain, and
// return "true" indicating that a new edge was created.
//------------------------------------------------------------------------------
bool
computeSortedEdgeNodes(std::vector<std::pair<int, int> >& edges,
                       std::vector<int>& result) {
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
  if (hangingNodes.size() == 2) {
    result.insert(result.begin() + 1, edges.size());
    POLY_ASSERT(hangingNodes[0] != hangingNodes[1]);
    edges.push_back(internal::hashEdge(hangingNodes[0], hangingNodes[1]));
    ++nedges;
    POLY_ASSERT(result.size() == 3);
  }
  POLY_ASSERT(edges.size() == nedges);

  // Pick a node to start the chain.
  if (hangingNodes.size() == 2) {
    POLY_ASSERT(nodeUseCount[edges[result.back()].first] == 2 or
                nodeUseCount[edges[result.back()].second] == 2);
    lastNode = (nodeUseCount[edges[result.back()].first] == 2 ? 
                edges[result.back()].first :
                edges[result.back()].second);
  } else {
    lastNode = edges[0].first;
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
  
  // Set the orientation for the ordered edges.
  lastNode = (edges[result[0]].first == edges[result[1]].first ? edges[result[0]].first : edges[result[0]].second);
  for (i = 1; i != nedges; ++i) {
    POLY_ASSERT(edges[result[i]].first == lastNode or edges[result[i]].second == lastNode);
    if (edges[result[i]].first == lastNode) {
      lastNode = edges[result[i]].second;
    } else {
      lastNode = edges[result[i]].first;
      result[i] = ~result[i];
    }
  }

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
template<typename RealType, typename IntType>
ReducedPLC<2, IntType>
hashReducedPLC(const ReducedPLC<2, RealType>& plc,
               const internal::QuantTessellation<2, RealType>& qmesh) {
  typedef geometry::Hasher<2, RealType> HasherType;
  ReducedPLC<2, IntType> result;
  result.facets = plc.facets;
  result.holes  = plc.holes;
  result.points.resize(plc.points.size());
  for (unsigned i = 0; i != plc.points.size()/2; ++i) {
    const IntType ip = qmesh.hashPosition(&plc.points[2*i]);
    result.points[2*i  ] = HasherType::qxval(ip);
    result.points[2*i+1] = HasherType::qyval(ip);
  }
  return result;
}


//------------------------------------------------------------------------------
// Build a Real ReducedPLC from a hashed int ReducedPLC.
//------------------------------------------------------------------------------
template<typename RealType, typename IntType>
ReducedPLC<2, RealType>
unhashReducedPLC(const ReducedPLC<2, IntType>& plc,
                 const internal::QuantTessellation<2, RealType>& qmesh) {
  typedef geometry::Hasher<2, RealType> HasherType;
  ReducedPLC<2, RealType> result;
  result.facets = plc.facets;
  result.holes  = plc.holes;
  result.points.resize(plc.points.size());
  for (unsigned i = 0; i != plc.points.size()/2; ++i) {
    const IntType ip = HasherType::hash(plc.points[2*i], plc.points[2*i+1]);
    qmesh.unhashPosition(ip, &result.points[2*i]);
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
    const int iedge = qmesh.faces[iface][0];
    POLY_ASSERT(iedge >= 0);
    const int ip = flip ? qmesh.edges[iedge].first : qmesh.edges[iedge].second;
    POLY_ASSERT(ip >= 0);
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


//------------------------------------------------------------------------------
// Build a ReducedPLC representation of a 2D cell.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<2, int64_t>
plcOfIntCell(const internal::QuantTessellation<2, RealType>& qmesh,
	     const unsigned icell) {
  POLY_ASSERT(icell < qmesh.cells.size());
  typedef typename internal::QuantTessellation<2, RealType>::CoordHash CoordHash;
  typedef typename internal::QuantTessellation<2, RealType>::IntPoint PointType;
  typedef typename geometry::Hasher<2, RealType> HasherType;
  ReducedPLC<2, CoordHash> result;
  const unsigned nFaces = qmesh.cells[icell].size();
  result.facets.resize(nFaces, vector<int>(2));
  for (unsigned i = 0; i != nFaces; ++i) {
    const bool flip = qmesh.cells[icell][i] < 0;
    const unsigned iface = flip ? ~qmesh.cells[icell][i] : qmesh.cells[icell][i];
    POLY_ASSERT(iface < qmesh.faces.size());
    POLY_ASSERT(qmesh.faces[iface].size() == 1);
    const int iedge = qmesh.faces[iface][0];
    POLY_ASSERT(iedge >= 0);
    const int ip = flip ? qmesh.edges[iedge].first : qmesh.edges[iedge].second;
    POLY_ASSERT(ip >= 0);
    POLY_ASSERT(ip < qmesh.points.size());
    PointType p = qmesh.hashedPosition(qmesh.points[ip]);
    result.points.push_back(p.x);
    result.points.push_back(p.y);
    result.facets[i][0] = i;
    result.facets[i][1] = (i+1)%nFaces;
  }
  POLY_ASSERT(result.points.size()/2 == nFaces);
  return result;
}

//------------------------------------------------------------------------------
// Figure out if any of the quantized points of this cell are outside the
// inner bounding box.
//------------------------------------------------------------------------------
template<typename RealType>
bool
hasOuterPoints(const internal::QuantTessellation<2, RealType>& qmesh,
	       const unsigned icell) {
  bool onlyInner = true;
  for (unsigned i = 0; i != qmesh.cells[icell].size(); ++i) {
    const bool flip = qmesh.cells[icell][i] < 0;
    const unsigned iface = flip ? ~qmesh.cells[icell][i] : qmesh.cells[icell][i];
    POLY_ASSERT(iface < qmesh.faces.size());
    POLY_ASSERT(qmesh.faces[iface].size() == 1);
    const int iedge = qmesh.faces[iface][0];
    POLY_ASSERT(iedge >= 0);
    const int ip = flip ? qmesh.edges[iedge].first : qmesh.edges[iedge].second;
    onlyInner *= (qmesh.points[ip] < (1ULL << 63));
  }
  return not onlyInner;
}

//------------------------------------------------------------------------------
// Construct a ReducedPLC from a range of coordinates.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<2,RealType>
plcOfBoundingBox(RealType* low, RealType* high) {
  ReducedPLC<2,RealType> result;
  result.points.resize(8);
  result.points[0] = low [0];  result.points[1] = low [1];
  result.points[2] = high[0];  result.points[3] = low [1];
  result.points[4] = high[0];  result.points[5] = high[1];
  result.points[6] = low [0];  result.points[7] = high[1];
  result.facets.resize(4, std::vector<int>(2));
  for (unsigned i=0; i!=4; ++i) {
    result.facets[i][0] = i;
    result.facets[i][1] = (i+1) % 4;
  }
  return result;
}



} // end anonymous namespace





template<typename RealType>
void
computeUnboundedQuantizedTessellationCollinear(const vector<RealType>& points,
					       const vector<RealType>& nonGeneratingPoints,
					       internal::QuantTessellation<2, RealType>& qmesh) {
  // Typedefs
  typedef Point2<RealType> RealPoint;
  typedef pair<int, int> EdgeHash;

  const unsigned N = points.size()/2;
  RealPoint p1, p2, r1, r2, node, midpt;
  bool test;
  int i, inode, icell1, icell2;

  // The center of the domain
  const RealPoint center = (qmesh.low_labframe + qmesh.high_labframe)/2;
  const RealType rinf = 4.0*std::max(qmesh.high_labframe.x - qmesh.low_labframe.x,
				     qmesh.high_labframe.y - qmesh.low_labframe.y);

  qmesh.low_inner  = center - RealPoint(rinf,rinf);
  qmesh.high_inner = center + RealPoint(rinf,rinf);

  qmesh.low_outer  = center - RealPoint(rinf,rinf)*1.05;
  qmesh.high_outer = center + RealPoint(rinf,rinf)*1.05;

  POLY_ASSERT(qmesh.low_inner.x  <= qmesh.high_inner.x and
	      qmesh.low_inner.y  <= qmesh.high_inner.y and
	      qmesh.low_outer.x  <= qmesh.high_outer.x and
	      qmesh.low_outer.y  <= qmesh.high_outer.y);
  POLY_ASSERT(qmesh.low_inner.x  >= qmesh.low_outer.x  and
	      qmesh.low_inner.y  >= qmesh.low_outer.y );
  POLY_ASSERT(qmesh.high_inner.x <= qmesh.high_outer.x and
	      qmesh.high_inner.y <= qmesh.high_outer.y);
  
  // Order the generators by position 0,...,N
  vector<pair<RealPoint, int> > pointIndexPairs;
  for (i = 0; i != N; ++i){
    pointIndexPairs.push_back(std::make_pair(RealPoint(points[2*i], points[2*i+1]), i));
    
  }
  sort(pointIndexPairs.begin(), pointIndexPairs.end(),
       internal::pairCompareFirst<RealPoint,int> );
  
  // Size the quant tessellation
  qmesh.edges.resize(3*N-1);
  qmesh.faces.resize(3*N-1, vector<int>(1));
  qmesh.cells.resize(N, vector<int>());

  // ------ Nodes and edges for min generator's cell ---------- //

  {
    // The left and right generators
    p1     = pointIndexPairs[0].first;
    p2     = pointIndexPairs[1].first;
    icell1 = pointIndexPairs[0].second;
    icell2 = pointIndexPairs[1].second;
    
    // Midpoint and direction vectors
    midpt = (p1 + p2)/2;
    r1.x = p2.x - p1.x;
    r1.y = p2.y - p1.y;
    geometry::unitVector<2,RealType>(&r1.x);
    r2.x =  r1.y;
    r2.y = -r1.x;
    
    // Extra inf node used to bound the first cell
    r1 *= -1.0;
    test = geometry::rayCircleIntersection(&p1.x, &r1.x, &center.x, 
					   rinf, 1.0e-10, &node.x);
    POLY_ASSERT(test);
    inode = qmesh.addNewNode(node);
    POLY_ASSERT(inode == 0);
    qmesh.infNodes.push_back(inode);
    
    // Node 1: endpt of first interior face
    test = geometry::rayCircleIntersection(&midpt.x, &r2.x, &center.x, 
					   rinf, 1.0e-10, &node.x);
    POLY_ASSERT(test);
    inode = qmesh.addNewNode(node);
    POLY_ASSERT(inode == 1);
    qmesh.infNodes.push_back(inode);  
    
    // Node 2: other endpt of first interior face
    r2 *= -1.0;
    test = geometry::rayCircleIntersection(&midpt.x, &r2.x, &center.x, 
					   rinf, 1.0e-10, &node.x);
    POLY_ASSERT(test);
    inode = qmesh.addNewNode(node);
    POLY_ASSERT(inode == 2);
    qmesh.infNodes.push_back(inode);
    
    // Register the edges and redundant faces
    qmesh.edges[0] = internal::hashEdge(0,1);
    qmesh.edges[1] = internal::hashEdge(1,2);
    qmesh.edges[2] = internal::hashEdge(0,2);
    qmesh.faces[0][0] = 0;
    qmesh.faces[1][0] = 1;
    qmesh.faces[2][0] = 2;

    // All the faces around cell 0
    qmesh.cells[icell1].push_back( 0);
    qmesh.cells[icell1].push_back( 1);
    qmesh.cells[icell1].push_back(~2);
      
    // Start the faces around cell 1
    qmesh.cells[icell2].push_back(~1);
  }

  // ------ The interior generators b/w min and max ------- //

  for (i = 1; i != N-1; ++i){
    // The left and right generators
    p1     = pointIndexPairs[i  ].first;
    p2     = pointIndexPairs[i+1].first;
    icell1 = pointIndexPairs[i  ].second;
    icell2 = pointIndexPairs[i+1].second;
    
    // Midpoint and direction vectors
    midpt = (p1 + p2)/2;
    r1.x = p2.x - p1.x;
    r1.y = p2.y - p1.y;
    geometry::unitVector<2,RealType>(&r1.x);
    r2.x =  r1.y;
    r2.y = -r1.x;
    
    // Node 0: endpt of interior face
    test = geometry::rayCircleIntersection(&midpt.x, &r2.x, &center.x, 
					   rinf, 1.0e-10, &node.x);
    POLY_ASSERT(test);
    inode = qmesh.addNewNode(node);
    POLY_ASSERT(inode == 2*i+1);
    qmesh.infNodes.push_back(inode);  
    
    // Node 1: other endpt of interior face
    r2 *= -1.0;
    test = geometry::rayCircleIntersection(&midpt.x, &r2.x, &center.x, 
					   rinf, 1.0e-10, &node.x);
    POLY_ASSERT(test);
    inode = qmesh.addNewNode(node);
    POLY_ASSERT(inode == 2*i+2);
    qmesh.infNodes.push_back(inode);  

    // The edges around the bottom, right, and top of cell i
    qmesh.edges[3*i  ] = internal::hashEdge(2*i-1,2*i+1);
    qmesh.edges[3*i+1] = internal::hashEdge(2*i+1,2*i+2);
    qmesh.edges[3*i+2] = internal::hashEdge(2*i  ,2*i+2);
    qmesh.faces[3*i  ][0] = 3*i;
    qmesh.faces[3*i+1][0] = 3*i+1;
    qmesh.faces[3*i+2][0] = 3*i+2;

    // The rest of the faces around cell i
    qmesh.cells[icell1].push_back(3*i);
    qmesh.cells[icell1].push_back(3*i+1);
    qmesh.cells[icell1].push_back(~(3*i+2));

    // Start the faces for cell i+1
    qmesh.cells[icell2].push_back(~(3*i+1));
  }
 
  // ------ Nodes and edges for max generator's cell ---------- //

  { 
    // Last and next-to-last generators
    p1     = pointIndexPairs[N-1].first;
    p2     = pointIndexPairs[N-2].first;
    icell1 = pointIndexPairs[N-1].second;
    
    // Direction vector
    r1.x = p1.x - p2.x;
    r1.y = p1.y - p2.y;
    geometry::unitVector<2,RealType>(&r1.x);
    
    // Inf node to bound the last generator
    test = geometry::rayCircleIntersection(&p2.x, &r1.x, &center.x, 
					   rinf, 1.0e-10, &node.x);
    POLY_ASSERT(test);
    inode = qmesh.addNewNode(node);
    POLY_ASSERT(inode == 2*N-1);
    qmesh.infNodes.push_back(inode);  
    
    // The last two inf edges/faces
    qmesh.edges[3*N-3] = internal::hashEdge(2*N-3,2*N-1);
    qmesh.edges[3*N-2] = internal::hashEdge(2*N-2,2*N-1);
    qmesh.faces[3*N-3][0] = 3*N-3;
    qmesh.faces[3*N-2][0] = 3*N-2;

    // Finish the face list for the final cell
    qmesh.cells[icell1].push_back(3*N-3);
    qmesh.cells[icell1].push_back(~(3*N-2));
  }
}








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

  // qmesh0.clipToInnerBoundingBox();

#if DUMP_UNBOUNDED_DEBUG_MESH and HAVE_SILO
  {
    Tessellation<2,RealType> debugMesh;
    qmesh0.tessellation(debugMesh);
    vector<double> px(debugMesh.cells.size());
    vector<double> py(debugMesh.cells.size());
    for (int ii = 0; ii != debugMesh.cells.size(); ++ii) {
      px[ii] = points[2*ii];
      py[ii] = points[2*ii+1];
    }
    map<string, double*> fields, cellFields;
    cellFields["gen_x"] = &px[0];
    cellFields["gen_y"] = &py[0];
    SiloWriter<2,RealType>::write(debugMesh, fields, fields, fields, cellFields, "debugMesh_Unbounded");
  }
#endif


  POLY_BEGIN_CONTRACT_SCOPE;
  {
    for (int i=0; i<numGenerators; ++i) {
      POLY_ASSERT(qmesh0.low_inner.x <= points[2*i  ] and points[2*i  ] <= qmesh0.high_inner.x);
      POLY_ASSERT(qmesh0.low_inner.y <= points[2*i+1] and points[2*i+1] <= qmesh0.high_inner.y);
      POLY_ASSERT(qmesh0.low_outer.x <= points[2*i  ] and points[2*i  ] <= qmesh0.high_outer.x);
      POLY_ASSERT(qmesh0.low_outer.y <= points[2*i+1] and points[2*i+1] <= qmesh0.high_outer.y);
    }
    for (int i=0; i<geometry.points.size()/2; ++i) {
      POLY_ASSERT(qmesh0.low_inner.x <= geometry.points[2*i  ] and geometry.points[2*i  ] <= qmesh0.high_inner.x);
      POLY_ASSERT(qmesh0.low_inner.y <= geometry.points[2*i+1] and geometry.points[2*i+1] <= qmesh0.high_inner.y);
      POLY_ASSERT(qmesh0.low_outer.x <= geometry.points[2*i  ] and geometry.points[2*i  ] <= qmesh0.high_outer.x);
      POLY_ASSERT(qmesh0.low_outer.y <= geometry.points[2*i+1] and geometry.points[2*i+1] <= qmesh0.high_outer.y);
    }
  }
  POLY_END_CONTRACT_SCOPE;
  

  // Create a new QuantTessellation.  This one will only use the single level of
  // quantization since we know the PLC is within this inner region.
  //
  // NOTE: I think this is in error. We can't get by with a single level of quantization
  //       because the unbounded data has two levels. We first have to clip the unbounded
  //       data within the outer quantization to the inner quantization boundary first.
  internal::QuantTessellation<2, RealType> qmesh1;
  qmesh1.generators    = qmesh0.generators;
  qmesh1.low_labframe  = qmesh0.low_labframe;
  qmesh1.high_labframe = qmesh0.high_labframe;
  qmesh1.low_inner     = qmesh0.low_inner;
  qmesh1.high_inner    = qmesh0.high_inner;
  qmesh1.low_outer     = qmesh0.low_outer;
  qmesh1.high_outer    = qmesh0.high_outer;
  // qmesh1.low_outer     = qmesh0.low_inner;
  // qmesh1.high_outer    = qmesh0.high_inner;
  qmesh1.degeneracy    = qmesh0.degeneracy;


#if HAVE_BOOST
  ReducedPLC<2, CoordHash> IntGeometry;
#if NORMALIZED_GENERATORS
  ReducedPLC<2, RealType> normalizedGeometry;
  normalizedGeometry.facets = geometry.facets;
  normalizedGeometry.holes  = geometry.holes;
  normalizedGeometry.points = this->computeNormalizedPoints(geometry.points,
							    geometry.points,
							    false,
							    &qmesh1.low_labframe.x,
							    &qmesh1.high_labframe.x);
  IntGeometry = hashReducedPLC<RealType, CoordHash>(normalizedGeometry, qmesh0);
#else
  IntGeometry = hashReducedPLC<RealType, CoordHash>(geometry, qmesh0);
#endif //NORMALIZED_GENERATORS




  const ReducedPLC<2, RealType> innerBoundingBox = plc_box<2,RealType>(&qmesh0.low_inner.x, 
								       &qmesh0.high_inner.x);
  std::vector<ReducedPLC<2, CoordHash> > allOrphans;
#endif //HAVE_BOOST


  // Walk the cells in the unbounded tessellation
  for (unsigned icell = 0; icell != numGenerators; ++icell) {
    
    // ---------------------------------------------------------
    // Intersect cell with boundary
    // ---------------------------------------------------------
    // Do the clipping integers if Boost.Geometry is available. 
    // Otherwise, reduce to using CSG in floating point.

    // std::cerr << "\n------------------------------ Clipping cell " << icell << std::endl;
    // std::cerr << "  \nPre-clipped cell:\n" << plcOfCell(qmesh0, icell) << std::endl;

    
    ReducedPLC<2, RealType> cell;
#if HAVE_BOOST
    const CoordHash pid = qmesh0.hashPosition(const_cast<RealType*>(&qmesh0.generators[2*icell]));
    const IntPoint IntGenerator = qmesh0.hashedPosition(pid);
    ReducedPLC<2, CoordHash> IntCell;

    // Check if any points for this cell are quantized based on the outer bounding box.
    // Clip these cells (in floating point) to the inner bounding box first.
    bool outerPoints = hasOuterPoints(qmesh0, icell);
    if (outerPoints) {
      cell = plcOfCell(qmesh0, icell);
      std::vector<ReducedPLC<2, RealType> > realOrphans;
      ReducedPLC<2, RealType> clippedRealCell = BG::boost_clip(innerBoundingBox,
							       cell,
							       RealPoint(qmesh0.generators[2*icell],
									 qmesh0.generators[2*icell+1]),
							       realOrphans);      


      // std::cerr << "  \nClipped-to-inner cell:\n" << clippedRealCell << std::endl;


      IntCell = hashReducedPLC<RealType, CoordHash>(clippedRealCell, qmesh0);
      POLY_ASSERT(realOrphans.size() == 0);
    } else {
      IntCell = plcOfIntCell(qmesh0, icell);
    }

    // std::cerr << "   Generator: (" << qmesh1.generators[2*icell] << "," << qmesh1.generators[2*icell+1] << ")"
    // 	      << " --> " << pid << " --> " << IntPoint(HasherType::qxval(pid), HasherType::qyval(pid)) << std::endl;
    // std::cerr << "  \nPre-clipped IntCell:\n" << IntCell << std::endl;
    // std::cerr << "   Clip... ";

    std::vector<ReducedPLC<2, CoordHash> > orphans;
    ReducedPLC<2, CoordHash> ClippedIntCell = BG::boost_clip(IntGeometry, IntCell, IntGenerator, orphans);
    allOrphans.insert(allOrphans.end(), orphans.begin(), orphans.end());
    POLY_ASSERT(ClippedIntCell.facets.size() >= 3);

    // std::cerr << "DONE!" << std::endl;
    // std::cerr << "  \nPost-clipped IntCell:\n" << ClippedIntCell << std::endl;

    // Add cell and its elements to the new tessellation
    vector<int> nodeIDs, edgeIDs, faceIDs;
    qmesh1.cells.push_back(vector<int>());
    for (unsigned i = 0; i != ClippedIntCell.points.size()/2; ++i) {
      nodeIDs.push_back(qmesh1.addNewNode(ClippedIntCell.points[2*i], ClippedIntCell.points[2*i+1]));
    }
    for (unsigned iface = 0; iface != ClippedIntCell.facets.size(); ++iface) {
      const unsigned nnodes = ClippedIntCell.facets[iface].size();
      POLY_ASSERT(nnodes == 2);
      vector<int> face;
      for (unsigned i = 0; i != nnodes; ++i) {
	const unsigned j = (i+1) % nnodes;
	POLY_ASSERT(nodeIDs[ClippedIntCell.facets[iface][i]] != nodeIDs[ClippedIntCell.facets[iface][j]]);
	const EdgeHash ehash = internal::hashEdge(nodeIDs[ClippedIntCell.facets[iface][i]],
						  nodeIDs[ClippedIntCell.facets[iface][j]]);
	face.push_back(qmesh1.addNewEdge(ehash));
	if (ehash.first == nodeIDs[ClippedIntCell.facets[iface][j]]) face.back() = ~face.back();
      }
      POLY_ASSERT(face.size() == nnodes);
      const unsigned k = qmesh1.faces.size();
      const unsigned i = qmesh1.addNewFace(face);
      qmesh1.cells.back().push_back(i == k ? i : ~i);
    }
    POLY_ASSERT(qmesh1.cells.back().size() == ClippedIntCell.facets.size());

#else

    // Build a ReducedPLC to represent the cell
    cell = plcOfCell(qmesh0, icell);
    cell = CSG::csg_intersect(geometry, cell);
    cell = simplifyPLCfacets(cell, 
			     cell.points,
			     &qmesh1.low_inner.x,
			     &qmesh1.high_inner.x,
			     1.0e-5);
    POLY_ASSERT(cell.facets.size() >= 3);

    // Add cell and its elements to the new tessellation
    vector<int> nodeIDs, edgeIDs, faceIDs;
    qmesh1.cells.push_back(vector<int>());
    for (unsigned i = 0; i != cell.points.size()/2; ++i) {
      nodeIDs.push_back(qmesh1.addNewNode(&cell.points[2*i]));
    }
    for (unsigned iface = 0; iface != cell.facets.size(); ++iface) {
      const unsigned nnodes = cell.facets[iface].size();
      POLY_ASSERT(nnodes == 2);
      vector<int> face;
      for (unsigned i = 0; i != nnodes; ++i) {
	const unsigned j = (i+1) % nnodes;
	POLY_ASSERT(nodeIDs[cell.facets[iface][i]] != nodeIDs[cell.facets[iface][j]]);
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
    
#endif //HAVE_BOOST
  }
    
    
    // cerr << "\n\n" << qmesh1 << "\n\n" << endl;


#if DUMP_UNBOUNDED_DEBUG_MESH and HAVE_SILO
  {
    Tessellation<2,RealType> debugMesh;
    qmesh1.tessellation(debugMesh);
    vector<double> px(debugMesh.cells.size());
    vector<double> py(debugMesh.cells.size());
    for (int ii = 0; ii != debugMesh.cells.size(); ++ii) {
      px[ii] = points[2*ii];
      py[ii] = points[2*ii+1];
    }
    map<string, double*> fields, cellFields;
    cellFields["gen_x"] = &px[0];
    cellFields["gen_y"] = &py[0];
    SiloWriter<2,RealType>::write(debugMesh, fields, fields, fields, cellFields, "debugMesh_Bounded");
  }
#endif

  // Blago!
  if (not allOrphans.empty()) {
#if HAVE_BOOST and false
    BoostOrphanage<RealType> orphanage(this);
#else
    cerr << "Orphans detected, but no actions taken" << endl;
#endif
  }
  // Blago!


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


#if NORMALIZE_GENERATORS
  qmesh.generators = this->computeNormalizedPoints(points, 
						   nonGeneratingPoints, 
						   true,
						   &qmesh.low_labframe.x, 
						   &qmesh.high_labframe.x);
#else
  std::vector<RealType> tmp = this->computeNormalizedPoints(points, 
							    nonGeneratingPoints, 
							    true,
							    &qmesh.low_labframe.x, 
							    &qmesh.high_labframe.x);
  qmesh.low_inner     = qmesh.low_labframe;
  qmesh.high_inner    = qmesh.high_labframe;
  qmesh.low_outer     = qmesh.low_labframe;
  qmesh.high_outer    = qmesh.high_labframe;
  qmesh.low_labframe  = RealPoint(0.0, 0.0);
  qmesh.high_labframe = RealPoint(1.0, 1.0);
  qmesh.generators = points;
#endif

  // Check for collinearity and use the appropriate routine
  bool isCollinear = geometry::collinear<2,RealType>(points, 1.0e-10);

  // Call a special routine to build up the quantized tessellation if the
  // input points are really 1D. This routine is purely geometric and is
  // independent of tessellator.
  if (isCollinear) {
    computeUnboundedQuantizedTessellationCollinear(qmesh.generators, nonGeneratingPoints, qmesh);
  } 

  // It's a fully-2D problem. Do the tessellator-specific stuff
  else {
    
    // Normalize the input generators
    const unsigned numGenerators = points.size()/2;
    unsigned i, j, k;

    // Call the underlying Delaunay algorithm and get its connectivity
    vector<RealPoint> circumcenters;
    vector<unsigned> triMask;
    map<EdgeHash, vector<unsigned> > edge2tris;
    map<int, set<unsigned> > gen2tri;
    vector<int> triangleList;
    computeDelaunayConnectivity(qmesh.generators,
				circumcenters,
				triMask,
				edge2tris,
				gen2tri,
				triangleList,
				qmesh.low_inner,
				qmesh.high_inner,
				qmesh.low_outer,
				qmesh.high_outer);
    const unsigned numTriangles = triMask.size();
    POLY_ASSERT(numTriangles > 0);
    POLY_ASSERT(circumcenters.size() ==   numTriangles );
    POLY_ASSERT(triangleList.size()  == 3*numTriangles );
    POLY_ASSERT(gen2tri.size()       ==   numGenerators);

    // Expand the outer bounding box and choose infinite sphere radius
    qmesh.low_outer.x  = min(qmesh.low_outer.x , qmesh.low_inner.x );
    qmesh.low_outer.y  = min(qmesh.low_outer.y , qmesh.low_inner.y );
    qmesh.high_outer.x = max(qmesh.high_outer.x, qmesh.high_inner.x);
    qmesh.high_outer.y = max(qmesh.high_outer.y, qmesh.high_inner.y);
    RealType rinf = 1.5*max(qmesh.high_outer.x - qmesh.low_outer.x,
			    qmesh.high_outer.y - qmesh.low_outer.y);
    const RealPoint centroid_outer = (qmesh.low_outer + qmesh.high_outer)*0.5;
    qmesh.low_outer.x  = centroid_outer.x - 1.05*rinf;
    qmesh.low_outer.y  = centroid_outer.y - 1.05*rinf;
    qmesh.high_outer.x = centroid_outer.x + 1.05*rinf;
    qmesh.high_outer.y = centroid_outer.y + 1.05*rinf;

    // Quantize circumcenters and add map them to unique IDs
    map<int, unsigned> tri2id;
    for (i = 0; i != numTriangles; ++i) {
      if (triMask[i] == 1) {
	tri2id[i] = qmesh.addNewNode(circumcenters[i]);
      }
    }
    POLY_ASSERT(tri2id.size() == std::accumulate(triMask.begin(), triMask.end(), 0));

  
    // //Blago!
    // for (i=0; i<numGenerators; ++i) {
    //   cerr << "Generator " << i << " at " << qmesh.generators[2*i] << " " << qmesh.generators[2*i+1] << endl << "   ";
    //   for (std::set<unsigned>::iterator itr = gen2tri[i].begin(); itr != gen2tri[i].end(); ++itr) {
    // 	if (triMask[*itr] == 1) {
    // 	  cerr << "(" << *itr << "," << tri2id[*itr] << ")  ";
    // 	}else{
    // 	  cerr << "(" << *itr << ")  ";
    // 	}
    //   }
    //   cerr << endl;
    // }
    // //Blago!
    

    // The exterior edges of the triangularization have "unbounded" rays, originating
    // at the circumcenter of the corresponding triangle and passing perpendicular to
    // the edge. Find those surface edges and project unbounded rays through them.
    bool test;
    RealPoint ehat, pinf, tricent;
    map<EdgeHash, unsigned> projEdge2id;
    int p, q, r, i1, i2, ivert;
    qmesh.infNodes = vector<unsigned>();
    for (typename map<EdgeHash, vector<unsigned> >::const_iterator edgeItr = edge2tris.begin();
         edgeItr != edge2tris.end();
         ++edgeItr) {
      const EdgeHash& edge = edgeItr->first;
      const vector<unsigned>& tris = edgeItr->second;
      if (tris.size() == 1) {
        i = tris[0];
        POLY_ASSERT(i < numTriangles);
        i1 = edge.first;
        i2 = edge.second;
        findOtherTriIndex(&triangleList[3*i], i1, i2, ivert);
        computeEdgeUnitVector<RealType>(&qmesh.generators[2*i1],
					&qmesh.generators[2*i2],
					&qmesh.generators[2*ivert],
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
    
    // The faces corresponding to each triangle edge
    qmesh.faces.reserve(edge2tris.size());
    qmesh.cells = vector<vector<int> >(numGenerators);
    int iedge, iface;
    unsigned ii, jj;
    RealType vol;
    RealPoint n0, n1;
    EdgeHash pq, pr;
    vector<vector<EdgeHash> > cellInfEdges(numGenerators);
    for (map<int, set<unsigned> >::const_iterator genItr = gen2tri.begin();
         genItr != gen2tri.end(); 
         ++genItr) {
      p = genItr->first;
      const set<unsigned>& tris = genItr->second;
      POLY_ASSERT(p < numGenerators);
      vector<EdgeHash> meshEdges;
      for (set<unsigned>::const_iterator triItr = tris.begin();
           triItr != tris.end(); 
  	 ++triItr){
        i = *triItr;
        POLY_ASSERT(i < numTriangles);
        POLY_ASSERT(tri2id.find(i) != tri2id.end());
        ii = tri2id[i];
        
        // Get the other indices for this triangle, given one of its vertices p
        findOtherTriIndices(&triangleList[3*i], p, q, r);
	POLY_ASSERT(p != q and q != r);
        pq = internal::hashEdge(p,q);
        pr = internal::hashEdge(p,r);
        
	// //Blago!
	// cerr << p << "," << ii << ": " << q << " " << r << endl;
	// //Blago!

        // Is pq a surface edge?
        if (edge2tris[pq].size() == 1){
          POLY_ASSERT(edge2tris[pq][0] == i);
          POLY_ASSERT(projEdge2id.find(pq) != projEdge2id.end());
          jj = projEdge2id[pq];
	  POLY_ASSERT(jj != ii);
	  // cerr << "---" << ii << " " << jj << " " << k << endl;
	  meshEdges.push_back(internal::hashEdge(ii,jj));
        } else {
	  POLY_ASSERT((edge2tris[pq].size() == 2 and edge2tris[pq][0] == i)
		      or edge2tris[pq][1] == i);
          k = (edge2tris[pq][0] == i ? edge2tris[pq][1] : edge2tris[pq][0]);
	  POLY_ASSERT(tri2id.find(k) != tri2id.end());
          jj = tri2id[k];
	  // cerr << "+++" << ii << " " << jj << " " << k << endl;
	  if (jj != ii) meshEdges.push_back(internal::hashEdge(ii,jj));
        }
        
        // Is pr a surface edge?
        if (edge2tris[pr].size() == 1){
          POLY_ASSERT(edge2tris[pr][0] == i);
          POLY_ASSERT(projEdge2id.find(pr) != projEdge2id.end());
          jj = projEdge2id[pr];
	  POLY_ASSERT(ii != jj);
	  // cerr << "---" << ii << " " << jj << " " << k << endl;
          meshEdges.push_back(internal::hashEdge(ii,jj));
        } else {
	  POLY_ASSERT((edge2tris[pr].size() == 2 and edge2tris[pr][0] == i)
		      or edge2tris[pr][1] == i);
          k = (edge2tris[pr][0] == i ? edge2tris[pr][1] : edge2tris[pr][0]);
	  POLY_ASSERT(tri2id.find(k) != tri2id.end());
          jj = tri2id[k];
	  // cerr << "+++" << ii << " " << jj << " " << k << endl;
	  if (jj != ii) meshEdges.push_back(internal::hashEdge(ii,jj));
        }
      }
  
      // //Blago!
      // if (p == 90) {
      // 	cerr << "   Cell " << p << endl << "      ";
      // 	for (unsigned ii = 0; ii != meshEdges.size(); ++ii) {
      // 	  EdgeHash ed = meshEdges[ii];
      // 	  cerr << "(" << ed.first << "," << ed.second << ")  ";
      // 	}
      // 	cerr << endl;
      // }
      // //Blago!


      // Arrange the edges in the corectly sorted and sign oriented order
      sort(meshEdges.begin(), meshEdges.end());
      meshEdges.erase(unique(meshEdges.begin(), meshEdges.end()), meshEdges.end());
      if (meshEdges.size() > 1) {
        vector<int> edgeOrder;
        const bool infEdge = computeSortedEdgeNodes(meshEdges, edgeOrder);
	
        if (meshEdges.size() > 2) {
  	
	  // // BLAGO!
	  // if (p == 90) {
	  //   for (vector<int>::const_iterator itr = edgeOrder.begin();
	  // 	 itr != edgeOrder.end();
	  // 	 ++itr) {
	  //     const bool flip = (*itr < 0);
	  //     k = (flip ? ~(*itr) : *itr);
	  //     n0 = qmesh.nodePosition(meshEdges[k].first);
	  //     n1 = qmesh.nodePosition(meshEdges[k].second);
	  //     CoordHash phash0 = qmesh.points[meshEdges[k].first];
	  //     CoordHash phash1 = qmesh.points[meshEdges[k].second];
	  //     IntPoint in0 = qmesh.hashedPosition(phash0);
	  //     IntPoint in1 = qmesh.hashedPosition(phash1);
	  //     cerr << meshEdges[k].first << "," << meshEdges[k].second << endl
	  // 	   << "   " << n0 << " " << n1 << endl
	  // 	   << "   " << in0 << " " << in1 << endl
	  // 	   << "   " << phash0 << " " << phash1 << endl;
	  //   }
	  // }
	  // // BLAGO!


	  // Add the edges and faces to the quantized mesh. (They are equal in 2D.)
	  for (vector<int>::const_iterator itr = edgeOrder.begin();
	       itr != edgeOrder.end();
	       ++itr) {
	    const bool flip = (*itr < 0);
	    k = (flip ? ~(*itr) : *itr);
	    iedge = qmesh.addNewEdge(meshEdges[k]);
	    //int edgePair[2] = {meshEdges[k].first, meshEdges[k].second};
	    //vector<int> face(edgePair, edgePair+2);
	    vector<int> face(1, iedge);
	    iface = qmesh.addNewFace(face);
	    POLY_ASSERT(iedge == iface);
	    
	    // Determine the orientation of the face with respect to the cell
	    n0 = qmesh.nodePosition(meshEdges[k].first);
	    n1 = qmesh.nodePosition(meshEdges[k].second);
	    vol = geometry::triangleVolume2(&qmesh.generators[2*p], &n1.x, &n0.x);
	    POLY_ASSERT2(vol != 0.0,
			 "Bounded tessellate error: cell " << p << " with position ("
			 << qmesh.generators[2*p] << "," << qmesh.generators[2*p+1]
			 << ") has zero triangle volume with respect to edge "
			 << n0 << " , " << n1);
	    if (vol > 0.0) {
	      qmesh.cells[p].push_back(iface);
	    } else {
	      qmesh.cells[p].push_back(~iface);
	    }
	  }
	  
	  // Did we create a new infEdge? If so we know it was the second element
	  // in the ordered list.
	  if (0 and infEdge) {
	    k = internal::positiveID(edgeOrder[1]);
	    // iedge = qmesh.edge2id[meshEdges[k]];
	    iedge = qmesh.addNewEdge(meshEdges[k]);
	    qmesh.infEdges.push_back(iedge);
	    cellInfEdges[p].push_back(meshEdges[k]);
	    //unsigned edgePair[2] = {meshEdges[k].first, meshEdges[k].second};
	    //vector<int> face(edgePair, edgePair+2);
	    vector<int> face(1, iedge);
	    iface = qmesh.addNewFace(face);
	    //   	  iface = qmesh.face2id[face];
	    POLY_ASSERT(iface == iedge);
	    qmesh.infFaces.push_back(iface);
	    qmesh.cells[p].push_back(iface);
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
  
  }

  // Post-conditions
  qmesh.assertValid();
}  
//------------------------------------------------------------------------------




template<typename RealType>
void
TriangleTessellator<RealType>::
computeDelaunayConnectivity(const vector<RealType>& points,
			    vector<RealPoint>& circumcenters,
			    vector<unsigned>& triMask,
			    map<EdgeHash, vector<unsigned> >& edge2tris,
			    map<int, set<unsigned> >& gen2tri,
			    vector<int>& triangleList,
			    RealPoint& low_inner,
			    RealPoint& high_inner,
			    RealPoint& low_outer,
			    RealPoint& high_outer) const {
  const unsigned numGenerators = points.size()/2;

  // Compute the triangularization
  triangulateio delaunay;
  computeDelaunay(points, delaunay, low_inner, high_inner);
  low_outer.x  = std::min(low_outer.x , low_inner.x );
  low_outer.y  = std::min(low_outer.y , low_inner.y );
  high_outer.x = std::max(high_outer.x, high_inner.x);
  high_outer.y = std::max(high_outer.y, high_inner.y);
    
  // Find the circumcenters of each triangle, and build the set of triangles
  // associated with each generator.
  int i, p, q, r;
  EdgeHash pq, pr, qr;
  circumcenters.resize(delaunay.numberoftriangles);
  triMask.resize(delaunay.numberoftriangles, 0);
  triangleList.resize(3*delaunay.numberoftriangles);
  for (i = 0; i != delaunay.numberoftriangles; ++i) {
    p = delaunay.trianglelist[3*i  ];
    q = delaunay.trianglelist[3*i+1];
    r = delaunay.trianglelist[3*i+2];
    triangleList[3*i  ] = p;
    triangleList[3*i+1] = q;
    triangleList[3*i+2] = r;
    POLY_ASSERT2(orient2d(&delaunay.pointlist[2*p],
			  &delaunay.pointlist[2*q],
			  &delaunay.pointlist[2*r]) > 0,
		 "TriangleTessellator Error: Delaunay vertices are "
		 << "not in CCW order for triangle " << i);
    geometry::computeCircumcenter(&delaunay.pointlist[2*p],
				  &delaunay.pointlist[2*q],
				  &delaunay.pointlist[2*r],
				  &circumcenters[i].x);
    POLY_ASSERT(p != q and p != r and q != r);
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
      low_outer.x  = min(low_outer.x , circumcenters[i].x);
      low_outer.y  = min(low_outer.y , circumcenters[i].y);
      high_outer.x = max(high_outer.x, circumcenters[i].x);
      high_outer.y = max(high_outer.y, circumcenters[i].y);
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
    POLY_ASSERT(low_outer.x  <= high_outer.x);
    POLY_ASSERT(low_outer.y  <= high_outer.y);
    POLY_ASSERT(low_inner.x  >= low_outer.x );
    POLY_ASSERT(low_inner.y  >= low_outer.y );
    POLY_ASSERT(high_inner.x <= high_outer.x);
    POLY_ASSERT(high_inner.y <= high_outer.y);
  }
  POLY_END_CONTRACT_SCOPE;
      
  // Clean up.
  trifree((VOID*)delaunay.pointlist);
  trifree((VOID*)delaunay.pointmarkerlist);
  trifree((VOID*)delaunay.trianglelist);
  trifree((VOID*)delaunay.edgelist);
  trifree((VOID*)delaunay.edgemarkerlist);
  trifree((VOID*)delaunay.segmentlist);
  trifree((VOID*)delaunay.segmentmarkerlist);
}
      
      
//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeDelaunay(const vector<RealType>& points,
	        triangulateio& delaunay,
	        RealPoint& low,
                RealPoint& high) const {
  triangulateio in;
   
  // Find the range of the generator points.
  const unsigned numGenerators = points.size()/2;

  // Determine bounding box for points
  RealType low1[2], high1[2];
  geometry::computeBoundingBox<2,RealType>(points, true, low1, high1);
  low[0]  = min(low[0] , low1[0] );
  low[1]  = min(low[1] , low1[1] );
  high[0] = max(high[0], high1[0]);
  high[1] = max(high[1], high1[1]);

  RealType box [2] = {high[0] - low[0], high[1] - low[1]};
  const RealType boxsize = 8.0*max(box[0], box[1]);
  
  const RealType xmin = 0.5*(low[0] + high[0]) - boxsize;
  const RealType xmax = 0.5*(low[0] + high[0]) + boxsize;
  const RealType ymin = 0.5*(low[1] + high[1]) - boxsize;
  const RealType ymax = 0.5*(low[1] + high[1]) + boxsize;

  low[0]  = min(low[0] , xmin);
  low[1]  = min(low[1] , ymin);
  high[0] = max(high[0], xmax);
  high[1] = max(high[1], ymax);

  // Add the generators
  in.numberofpoints = numGenerators + 4;
  in.pointlist = new RealType[2*in.numberofpoints];
  copy(points.begin(), points.end(), in.pointlist);
  in.pointlist[2*numGenerators  ] = xmin;  in.pointlist[2*numGenerators+1] = ymin;
  in.pointlist[2*numGenerators+2] = xmax;  in.pointlist[2*numGenerators+3] = ymin;
  in.pointlist[2*numGenerators+4] = xmax;  in.pointlist[2*numGenerators+5] = ymax;
  in.pointlist[2*numGenerators+6] = xmin;  in.pointlist[2*numGenerators+7] = ymax;
  in.numberofsegments = 0;

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
