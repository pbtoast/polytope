//------------------------------------------------------------------------
// TriangleTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <numeric>
#include <set>
#include <list>
#include <map>
#include <limits>

#include "polytope.hh" // Pulls in POLY_ASSERT and TriangleTessellator.hh.
#include "convexHull_2d.hh"
#include "nearestPoint.hh"
#include "within.hh"
#include "intersect.hh"
#include "PLC_Boost_2d.hh"
#include "polytope_plc_canned_geometries.hh"

// Pull in triangle
// Since triangle isn't built to work out-of-the-box with C++, we 
// slurp in its source here, bracketing it with the necessary dressing.
#define TRILIBRARY
#define REAL double
#define VOID void
#define ANSI_DECLARATORS 
#define CDT_ONLY // Conforming Delaunay triangulations only! 

#define ENABLE_INTEGER_INTERSECTIONS false

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

extern "C"
{
#include "triangle.h"
}


// Fast predicate for determining colinearity of points.o
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
// Sort a set of edges around a face so that sequential edges share nodes.
// We account for one break in the chain, representing an unbounded surface.
//------------------------------------------------------------------------------
vector<unsigned>
computeSortedFaceNodes(const vector<pair<int, int> >& edges) {
  typedef pair<int, int> EdgeHash;
  vector<unsigned> result;

  const unsigned nedges = edges.size();
  if (nedges > 1) {

    // Invert the mapping, from nodes to edges.
    std::map<int, std::set<EdgeHash> > nodes2edges;
    internal::CounterMap<int> nodeUseCount;
    unsigned i, j;
    for (i = 0; i != nedges; ++i) {
      nodes2edges[edges[i].first].insert(edges[i]);
      nodes2edges[edges[i].second].insert(edges[i]);
      ++nodeUseCount[edges[i].first];
      ++nodeUseCount[edges[i].second];
    }
    
    // Look for any edges with one node in the set.  There can be at most
    // two such edges, representing the two ends of the chain.  We will put 
    // the edges with those nodes first in the ordering, so that all remaining
    // edges should naturally hook together.
    std::vector<EdgeHash> orderedEdges;
    orderedEdges.reserve(nedges);
    int lastNode;
    bool hangingNodes = false;
    for (i = 0; i != nedges; ++i) {
      if (nodeUseCount[edges[i].first] == 1 or
          nodeUseCount[edges[i].second] == 1) {
        POLY_ASSERT((nodeUseCount[edges[i].first] == 1 and nodeUseCount[edges[i].second] == 2) or
                    (nodeUseCount[edges[i].first] == 2 and nodeUseCount[edges[i].second] == 1));
        orderedEdges.push_back(edges[i]);
        nodes2edges[edges[i].first].erase(edges[i]);
        nodes2edges[edges[i].second].erase(edges[i]);
        lastNode = (nodeUseCount[edges[i].first] == 1 ? edges[i].first : edges[i].second);
        hangingNodes = true;
      }
    }
    POLY_ASSERT(orderedEdges.size() == 0 or orderedEdges.size() == 2);

    // Pick a node to start the chain.
    if (hangingNodes) {
      POLY_ASSERT(nodeUseCount[orderedEdges.back().first] == 2 or
                  nodeUseCount[orderedEdges.back().second] == 2);
      lastNode = (nodeUseCount[orderedEdges.back().first] == 2 ? 
                  orderedEdges.back().first :
                  orderedEdges.back().second);
    } else {
      lastNode = edges[0].first;
    }

    // Walk the remaining edges
    while (orderedEdges.size() != nedges) {
      POLY_ASSERT(nodes2edges[lastNode].size() > 0);
      orderedEdges.push_back(*nodes2edges[lastNode].begin());
      nodes2edges[orderedEdges.back().first].erase(orderedEdges.back());
      nodes2edges[orderedEdges.back().second].erase(orderedEdges.back());
      lastNode = (orderedEdges.back().first == lastNode ? orderedEdges.back().second : orderedEdges.back().first);
    }

    // Read the nodes in order.
    if (hangingNodes) {
      result.push_back(nodeUseCount[orderedEdges[0].first] == 1 ? orderedEdges[0].first : orderedEdges[0].second);
      result.push_back(nodeUseCount[orderedEdges[1].first] == 1 ? orderedEdges[1].first : orderedEdges[1].second);
      i = 1;
    } else {
      i = 0;
    }
    for (; i != nedges; ++i) {
      j = (i + 1) % nedges;
      POLY_ASSERT(orderedEdges[i].first == orderedEdges[j].first or
                  orderedEdges[i].first == orderedEdges[j].second or
                  orderedEdges[i].second == orderedEdges[j].first or
                  orderedEdges[i].second == orderedEdges[j].second);
      result.push_back((orderedEdges[i].first == orderedEdges[j].first or orderedEdges[i].first == orderedEdges[j].second) ? 
                       orderedEdges[i].first : 
                       orderedEdges[i].second);
    }
    POLY_ASSERT2((hangingNodes and result.size() == nedges + 1) or
                 ((not hangingNodes) and result.size() == nedges), result.size());      

  } else {

    // There are either one or no edges, so the solution is pretty simple.
    if (nedges == 1) {
      result.push_back(edges[0].first);
      result.push_back(edges[0].second);
    }
  }
  // That's it.
  return result;
}


//------------------------------------------------------------------------------
// Output an pair of node-index lists as a ReducedPLC for the cell
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<2, RealType>
plcOfCell(const vector<unsigned>& nodeIndices,
          map<int, Point2<RealType> >& nodeMap) {
  typedef Point2<RealType> PointType;
  ReducedPLC<2, RealType> result;
  const unsigned nnodes = nodeIndices.size();
  result.points.resize(2*nnodes);
  result.facets.resize(nnodes, vector<int>(2));
  for (unsigned i = 0; i < nnodes; ++i) {
    result.facets[i][0] = i;
    result.facets[i][1] = (i + 1)%nnodes;
    POLY_ASSERT(nodeMap.find(nodeIndices[i]) != nodeMap.end());
    const PointType pt = nodeMap[nodeIndices[i]];
    result.points[2*i  ] = pt.x;
    result.points[2*i+1] = pt.y;
  }
  return result;
}

//------------------------------------------------------------------------------
// Output an pair of node-index lists as a ReducedPLC for the cell
//------------------------------------------------------------------------------
template<typename RealType, typename IntType>
ReducedPLC<2, IntType>
quantizeGeometry(const QuantizedCoordinates<2, RealType>& coords,
                 const ReducedPLC<2, RealType>& geometry) {
  typedef Point2<RealType> RealPoint;
  typedef Point2<IntType> IntPoint;
  ReducedPLC<2, IntType> result;
  result.facets = geometry.facets;
  result.holes  = geometry.holes;
  result.points.resize(geometry.points.size());
  for (unsigned i = 0; i < geometry.points.size()/2; ++i) {
    const RealPoint rp = RealPoint(geometry.points[2*i  ],
                                   geometry.points[2*i+1]);
    const IntPoint ip = coords.quantize(&rp.x);
    result.points[2*i  ] = ip.x;
    result.points[2*i+1] = ip.y;
  }
  return result;
}

//------------------------------------------------------------------------------
// Output an pair of node-index lists as a ReducedPLC for the cell
//------------------------------------------------------------------------------
template<typename RealType, typename IntType>
ReducedPLC<2, RealType>
dequantizeGeometry(const QuantizedCoordinates<2, RealType>& coords,
                   const ReducedPLC<2, IntType>& geometry) {
  typedef Point2<RealType> RealPoint;
  typedef Point2<IntType> IntPoint;
  ReducedPLC<2, RealType> result;
  result.facets = geometry.facets;
  result.holes  = geometry.holes;
  result.points.resize(geometry.points.size());
  for (unsigned i = 0; i < geometry.points.size()/2; ++i) {
    const RealPoint rp = coords.dequantize(&geometry.points[2*i]);
    result.points[2*i  ] = rp.x;
    result.points[2*i+1] = rp.y;
  }
  return result;
}

//------------------------------------------------------------------------------
// Collapse a facet having points within a given tolerance
//------------------------------------------------------------------------------
template<typename RealType>
void
simplifyGeometry(ReducedPLC<2, RealType>& geometry, const RealType tol) {
  typedef Point2<RealType> PointType;
  ReducedPLC<2, RealType> result;
  const int nFacets = geometry.facets.size();
  for (int i = 0; i < nFacets; ++i) {
    POLY_ASSERT(geometry.facets[i].size() == 2);
    const int i1 = geometry.facets[i][0];
    const int i2 = geometry.facets[i][1];
    POLY_ASSERT(i1 < geometry.points.size()/2 and i2 < geometry.points.size()/2);
    const PointType p1 = PointType(geometry.points[2*i1], geometry.points[2*i1+1]);
    const PointType p2 = PointType(geometry.points[2*i2], geometry.points[2*i2+1]);
    const double dx = double(p2.x - p1.x);
    const double dy = double(p2.y - p1.y);
    const double dist = sqrt(dx*dx + dy*dy);
    if (dist > tol) {
      result.points.push_back(geometry.points[2*i1]);
      result.points.push_back(geometry.points[2*i1+1]);
    }
  }
  POLY_ASSERT(result.points.size()/2 > 0 and result.points.size()/2 <= nFacets);
  result.facets.resize(result.points.size()/2, std::vector<int>(2));
  for (int i = 0; i < result.points.size()/2; ++i) {
    result.facets[i][0] = i;
    result.facets[i][1] = (i+1) % (result.points.size()/2);
  }
  geometry = result;
}

//------------------------------------------------------------------------------
// Exprss the Triangle-generated delaunay struct as a Tessellation
//------------------------------------------------------------------------------
template<typename RealType>
void
constructDelaunayMesh(const triangulateio& delaunay,
                      Tessellation<2, RealType>& dmesh) {  
  typedef pair<int, int> EdgeHash;

  int i;
  dmesh.nodes.resize(2*delaunay.numberofpoints);
  for (i = 0; i < 2*delaunay.numberofpoints; ++i) {
    dmesh.nodes[i] = delaunay.pointlist[i];
  }

  int iedge, p, q, r;
  EdgeHash pq, qr, rp;
  map<EdgeHash, int> edge2id;
  map<int, std::vector<int> > edgeCells;
  dmesh.cells.resize(delaunay.numberoftriangles, std::vector<int>(3));
  for (i = 0; i < delaunay.numberoftriangles; ++i) {
    p = delaunay.trianglelist[3*i  ];
    q = delaunay.trianglelist[3*i+1];
    r = delaunay.trianglelist[3*i+2];
    pq = internal::hashEdge(p,q);
    qr = internal::hashEdge(q,r);
    rp = internal::hashEdge(r,p);

    iedge = internal::addKeyToMap(pq, edge2id);
    edgeCells[iedge].push_back(p < q ? i : ~i);
    dmesh.cells[i][0] = (p < q ? iedge : ~iedge);

    iedge = internal::addKeyToMap(qr, edge2id);
    edgeCells[iedge].push_back(q < r ? i : ~i);
    dmesh.cells[i][2] = (q < r ? iedge : ~iedge);

    iedge = internal::addKeyToMap(rp, edge2id);
    edgeCells[iedge].push_back(r < p ? i : ~i);
    dmesh.cells[i][1] = (r < p ? iedge : ~iedge);
  }

  dmesh.faces.resize(edge2id.size(), vector<unsigned>(2));
  for (typename std::map<EdgeHash, int>::const_iterator itr = edge2id.begin();
       itr != edge2id.end(); 
       ++itr) {
    const EdgeHash& edge = itr->first;
    i = itr->second;
    dmesh.faces[i][0] = edge.first;
    dmesh.faces[i][1] = edge.second;
  }

  dmesh.faceCells.resize(edge2id.size());
  for (i = 0; i < dmesh.faces.size(); ++i) dmesh.faceCells[i] = edgeCells[i];
}

//------------------------------------------------------------------------------
// Output a Triangle-generated delaunay mesh into a silo file for viewing.
//------------------------------------------------------------------------------
template<typename RealType>
void
dumpDelaunay(const triangulateio& delaunay,
             string name,
             const int cycle) {
  Tessellation<2, RealType> dmesh;
  constructDelaunayMesh(delaunay, dmesh);
  map<string, double*> fields;
#if HAVE_SILO
  SiloWriter<2,RealType>::write(dmesh, fields, fields, fields, fields, name, cycle, 0.0);
#endif
}

//------------------------------------------------------------------------------
// Comparison operator for two points. They're equal if one lives inside the
// 3x3 region around the other:    _ _ _
//                                |_|_|_|
//                                |_|_|_|
//                                |_|_|_|
//------------------------------------------------------------------------------
template<typename RealType>
class ThreeByThreeCompare {
public:
  bool operator()(const Point2<RealType> pt1, const Point2<RealType> pt2) {
    return (pt1.x < pt2.x-1 ? true :
            (pt1.x == pt2.x-1 or pt1.x == pt2.x or pt1.x == pt2.x+1) and 
            pt1.y < pt2.y-1 ? true : false);
  }
};

//------------------------------------------------------------------------------
// Comparison operator for two points. They're equal if one lives inside the
// plus-sign-shaped region around the other:        _
//                                                _|_|_
//                                               |_|_|_|
//                                                 |_|
//------------------------------------------------------------------------------
template<typename RealType>
class PlusSignCompare {
public:
  bool operator()(const Point2<RealType> pt1, const Point2<RealType> pt2) {
    return (pt1.x < pt2.x-1 ? true :
            pt1.x == pt2.x-1 and pt1.y < pt2.y ? true :
            pt1.x == pt2.x   and pt1.y < pt2.y-1 ? true : false);
  }
};

//------------------------------------------------------------------------------
// Check the orientation of a cell. If CCW, return True.
//------------------------------------------------------------------------------
template<typename RealType>
bool
checkCellOrientation(const RealType* pt, const ReducedPLC<2, RealType>& cell) {
  double orientation = 0.0;
  double pta[2] = {double(pt[0]), double(pt[1])};
  for (int ifacet = 0; ifacet < cell.facets.size(); ++ifacet) {
    double ptb[2] = {double(cell.points[2*cell.facets[ifacet][0]  ]),
                     double(cell.points[2*cell.facets[ifacet][0]+1])};
    double ptc[2] = {double(cell.points[2*cell.facets[ifacet][1]  ]),
                     double(cell.points[2*cell.facets[ifacet][1]+1])};
    orientation += orient2d(pta, ptb, ptc);
  }
  return (orientation > 0.0);
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
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           Tessellation<2, RealType>& mesh) const {
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(mesh.empty());

  // Initialize the bounding box to quantize mesh nodes and project
  // unbounded mesh edges to infinity.
  mCoords.initialize(points, mDegeneracy);

  bool collinear = geometry::collinear<2,RealType>(points, mDegeneracy);
  vector<vector<unsigned> > cellNodes;
  map<int, IntPoint> id2node;
  vector<unsigned> infNodes;

  // The points are collinear. This is easy:
  if (collinear) {
    this->computeCellNodesCollinear(points, cellNodes, id2node, infNodes);
  }

  // The points aren't collinear. Do the complicated stuff.
  else {
    this->computeCellNodes(points, cellNodes, id2node, infNodes);
  }

  // Copy the quantized nodes to the final tessellation.
  const unsigned numNodes = id2node.size();
  RealPoint node;
  mesh.nodes.resize(2*numNodes);
  mesh.infNodes = infNodes;
  for (typename map<int, IntPoint>::const_iterator itr = id2node.begin();
       itr != id2node.end();
       ++itr) {
    node = mCoords.dequantize(&(itr->second).x);
    POLY_ASSERT(itr->first < numNodes);
    mesh.nodes[2*(itr->first)  ] = node.x;
    mesh.nodes[2*(itr->first)+1] = node.y;
  }

  //Finish constructing the cell-face-node-topology
  constructUnboundedMeshTopology(cellNodes, points, mesh);  
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           RealType* low,
           RealType* high,
           Tessellation<2, RealType>& mesh) const {
  POLY_ASSERT(low[0] <= high[0] and low[1] <= high[1]);
  ReducedPLC<2, RealType> box = plc_box<2, RealType>(low, high);
  this->tessellate(points, box, mesh);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const vector<RealType>& plcPoints,
           const PLC<2, RealType>& geometry,
           Tessellation<2, RealType>& mesh) const {
  ReducedPLC<2, RealType> boundary;
  boundary.facets = geometry.facets;
  boundary.holes  = geometry.holes;
  boundary.points = plcPoints;
  this->tessellate(points, boundary, mesh);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const ReducedPLC<2, RealType>& geometry,
           Tessellation<2, RealType>& mesh) const {
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(not geometry.points.empty());
  POLY_ASSERT(mesh.empty());

  // Initialize the bounding box to quantize mesh nodes and project
  // unbounded mesh edges to infinity.
  mCoords.initialize(geometry.points, mDegeneracy);

  const unsigned numGenerators = points.size()/2;
  bool collinear = geometry::collinear<2,RealType>(points, mDegeneracy);
  vector<vector<unsigned> > cellNodes;
  map<int, IntPoint> id2node;
  vector<unsigned> infNodes;   // dummy arg for bounded tessellations

  // The points are collinear. This is easy:
  if (collinear) {
    this->computeCellNodesCollinear(points, cellNodes, id2node, infNodes);
  }

  // The points aren't collinear. Do the complicated stuff.
  else {
    this->computeCellNodes(points, cellNodes, id2node, infNodes);
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);


  // Construct a PLC of the geometric boundary.
  vector<ReducedPLC<2, CoordHash> > intCells(numGenerators);
#if ENABLE_INTEGER_INTERSECTIONS
  const ReducedPLC<2, CoordHash> intGeometry = quantizeGeometry<RealType, CoordHash>(mCoords, geometry);
  vector<ReducedPLC<2, CoordHash> > orphans;
  for (int i = 0; i < numGenerators; ++i) {
    const IntPoint intGenerator = mCoords.quantize(&points[2*i]);
    const ReducedPLC<2, CoordHash> intCell = plcOfCell<CoordHash>(cellNodes[i], id2node);

    intCells[i] = BG::boost_clip<CoordHash>(intGeometry, intCell, intGenerator, orphans);
    POLY_ASSERT2(not BG::boost_intersects(intCells[i]),
                 "Cell " << i << " intersects itself:\n" 
                 << "\nBefore clipping:\n" << intCell
                 << "\nAfter clipping:\n" << intCells[i]
                 << "\nOuter Geometry:\n" << intGeometry
                 << "\nGenerator:\n" << intGenerator);
  }

  if (orphans.size() > 0) {
    cerr << "Orphans detected." << endl;
    BoostOrphanage<RealType> orphanage(this);
    orphanage.adoptOrphans(points, mCoords, intCells, orphans);
  }
#else
  vector<ReducedPLC<2, RealType> > orphans;
  for (int i = 0; i < numGenerators; ++i) {
    const RealPoint generator = RealPoint(points[2*i], points[2*i+1]);
    const ReducedPLC<2, CoordHash> intCell = plcOfCell<CoordHash>(cellNodes[i], id2node);
    const ReducedPLC<2, RealType> cell = dequantizeGeometry(mCoords, intCell);
    ReducedPLC<2, RealType> clippedCell;
    clippedCell = BG::boost_clip<RealType>(geometry, cell, generator, orphans);
    POLY_ASSERT2(not BG::boost_intersects(clippedCell), clippedCell);
    intCells[i] = quantizeGeometry<RealType, CoordHash>(mCoords, clippedCell);
    simplifyGeometry<CoordHash>(intCells[i], 0);
  }

  if (orphans.size() > 0) {
    cerr << "Orphans detected. Taking no action." << endl;
  }
#endif

  constructBoundedTopology(points, geometry, intCells, mesh);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeDelaunayConnectivity(const vector<RealType>& points,
			    vector<RealPoint>& circumcenters,
			    vector<unsigned>& triMask,
			    map<EdgeHash, vector<unsigned> >& edge2tris,
			    map<int, set<unsigned> >& gen2tri,
			    vector<int>& triangleList) const {
  const unsigned numGenerators = points.size()/2;

  // Compute the triangularization
  triangulateio delaunay;
  computeDelaunay(points, delaunay, true);
    
  // Find the circumcenters of each triangle, and build the set of triangles
  // associated with each generator.
  int i, p, q, r;
  EdgeHash pq, pr, qr;
  circumcenters.resize(delaunay.numberoftriangles);
  triMask.resize(delaunay.numberoftriangles, 0);
  triangleList.resize(3*delaunay.numberoftriangles);
  RealPoint lowc  = RealPoint(mCoords.low [0], mCoords.low [1]);
  RealPoint highc = RealPoint(mCoords.high[0], mCoords.high[1]);
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
      lowc.x  = min(lowc.x , circumcenters[i].x);
      lowc.y  = min(lowc.y , circumcenters[i].y);
      highc.x = max(highc.x, circumcenters[i].x);
      highc.y = max(highc.y, circumcenters[i].y);
    }
  }
  POLY_ASSERT(circumcenters.size() == delaunay.numberoftriangles);
  POLY_ASSERT(gen2tri.size() == numGenerators);
  POLY_ASSERT(std::accumulate(triMask.begin(), triMask.end(), 0) > 0);
  POLY_ASSERT(lowc.x <= highc.x and lowc.y <= highc.y);
  mCoords.expand(&lowc.x, &highc.x);
  
  // // Blago!
  // {
  //   triangulateio dtmp;
  //   computeDelaunay(points, dtmp, false);
  //   dumpDelaunay<RealType>(dtmp, "delaunayDEBUG", 0);
  //   trifree((VOID*)dtmp.pointlist);
  //   trifree((VOID*)dtmp.pointmarkerlist);
  //   trifree((VOID*)dtmp.trianglelist);
  //   trifree((VOID*)dtmp.edgelist);
  //   trifree((VOID*)dtmp.edgemarkerlist);
  //   trifree((VOID*)dtmp.segmentlist);
  //   trifree((VOID*)dtmp.segmentmarkerlist);
  //   computeDelaunay(points, dtmp, true);
  //   dumpDelaunay<RealType>(dtmp, "delaunayDEBUG", 1);
  //   trifree((VOID*)dtmp.pointlist);
  //   trifree((VOID*)dtmp.pointmarkerlist);
  //   trifree((VOID*)dtmp.trianglelist);
  //   trifree((VOID*)dtmp.edgelist);
  //   trifree((VOID*)dtmp.edgemarkerlist);
  //   trifree((VOID*)dtmp.segmentlist);
  //   trifree((VOID*)dtmp.segmentmarkerlist);
  // }
  // // Blago!

  POLY_BEGIN_CONTRACT_SCOPE;
  {
    for (typename map<EdgeHash, vector<unsigned> >::const_iterator itr = edge2tris.begin();
         itr != edge2tris.end();
	 ++itr) POLY_ASSERT(itr->second.size() == 1 or itr->second.size() == 2);
    for (map<int, set<unsigned> >::const_iterator itr = gen2tri.begin();
         itr != gen2tri.end();
	 ++itr) POLY_ASSERT(itr->second.size() >= 1);
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

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeCellNodesCollinear(const std::vector<RealType>& points,
                          vector<vector<unsigned> >& cellNodes,
                          map<int, IntPoint>& id2node,
                          vector<unsigned>& infNodes) const { 

  // A dummy expansion just in case we're doing a distributed mesh construction.
  mCoords.expand(&mCoords.low.x, &mCoords.high.x);

  // Call the 1d routine for projecting a line of points
  vector<RealPoint> nodes;
  constructCells1d(points, &(mCoords.center).x, mCoords.rinf, cellNodes, nodes);
  POLY_ASSERT(cellNodes.size() == points.size()/2);
  POLY_ASSERT(nodes.size() == points.size());

  // Quantize nodes and assign indices
  std::set<IntPoint> uniqueNodes;
  infNodes.resize(nodes.size());
  for (unsigned i = 0; i < nodes.size(); ++i) {
    const IntPoint ip = mCoords.quantize(&(nodes[i]).x);
    uniqueNodes.insert(ip);
    POLY_ASSERT(uniqueNodes.size() == i+1);
    id2node[i] = ip;
    infNodes[i] = i;   // All nodes are projected inf nodes
  }
  POLY_ASSERT(uniqueNodes.size() == nodes.size());
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeCellNodes(const vector<RealType>& points,
                 vector<vector<unsigned> >& cellNodes,
                 map<int, IntPoint>& id2node,
                 vector<unsigned>& infNodes) const { 
  const unsigned numGenerators = points.size()/2;
  int i, j, k;
  // Compute the minimum connectivity set from the Delaunay diagram.
  cellNodes.resize(numGenerators);
  vector<RealPoint> circumcenters;
  vector<unsigned> triMask;
  map<EdgeHash, vector<unsigned> > edge2tri;
  map<int, set<unsigned> > gen2tri;
  vector<int> triangleList;
  computeDelaunayConnectivity(points,
                              circumcenters,
                              triMask,
                              edge2tri,
                              gen2tri,
                              triangleList);
  const unsigned numTriangles = triMask.size();
  POLY_ASSERT(numTriangles > 0);
  POLY_ASSERT(circumcenters.size() ==   numTriangles );
  POLY_ASSERT(triangleList.size()  == 3*numTriangles );
  POLY_ASSERT(gen2tri.size()       ==   numGenerators);
  
  // Assign a unique ID to each triangle and circumcenter.
  // Circumcenters that overlap are collapsed into nodes.
  // map<IntPoint, int, ThreeByThreeCompare<CoordHash> > circ2id;
  map<IntPoint, int> circ2id;
  map<int, unsigned> tri2id;
  for (i = 0; i != numTriangles; ++i){
    if (triMask[i] == 1) {
      const IntPoint ip = mCoords.quantize(&circumcenters[i].x);
      POLY_ASSERT(-mCoords.mCoordMax <= ip.x and ip.x <= mCoords.mCoordMax and
                  -mCoords.mCoordMax <= ip.y and ip.y <= mCoords.mCoordMax);
      k = circ2id.size();
      j = internal::addKeyToMap(ip, circ2id);
      tri2id[i] = j;
      if (k != circ2id.size()) id2node[j] = ip;
    }
  }
  POLY_ASSERT(circ2id.size() == id2node.size());

  // Find all the infinite edges of the unbounded Voronoi. Project nodes at
  // "infinity" to represent unbounded Voronoi cells.
  int i1, i2, ivert;
  map<EdgeHash, unsigned> edge2id;
  for (map<EdgeHash, vector<unsigned> >::const_iterator edgeItr = edge2tri.begin();
       edgeItr != edge2tri.end(); 
       ++edgeItr){
    const EdgeHash& edge = edgeItr->first;
    const vector<unsigned>& tris = edgeItr->second;
    if (tris.size() == 1){
      i = tris[0];
      POLY_ASSERT(i < numTriangles);
      i1 = edge.first;
      i2 = edge.second;
      POLY_ASSERT(i1 != i2);
      findOtherTriIndex(&triangleList[3*i], i1, i2, ivert);
      POLY_ASSERT(i1 < numGenerators and i2 < numGenerators and ivert < numGenerators);

      RealPoint ehat;
      RealPoint p1   (points[2*i1]   , points[2*i1   +1]);
      RealPoint p2   (points[2*i2]   , points[2*i2   +1]);
      RealPoint pvert(points[2*ivert], points[2*ivert+1]);
      computeEdgeUnitVector(&p1.x, &p2.x, &pvert.x, &ehat.x);

      RealPoint pinf = mCoords.projectPoint(&circumcenters[i].x, &ehat.x);
      const IntPoint ip = mCoords.quantize(&pinf.x);
      POLY_ASSERT(-mCoords.mCoordMax <= ip.x and ip.x <= mCoords.mCoordMax and
                  -mCoords.mCoordMax <= ip.y and ip.y <= mCoords.mCoordMax);
      k = circ2id.size();
      j = internal::addKeyToMap(ip, circ2id);
      if (k != circ2id.size())  id2node[j] = ip;
      POLY_ASSERT(edge2id.find(edge) == edge2id.end());
      edge2id[edge] = j;
      infNodes.push_back(j);
    }
  }

  // Walk triangles around each generator and build up the cell edges
  int p, q, r;
  EdgeHash pq, pr;
  unsigned ii, jj;
  for (map<int, set<unsigned> >::const_iterator genItr = gen2tri.begin();
       genItr != gen2tri.end(); 
       ++genItr) {
    p = genItr->first;
    POLY_ASSERT(p < numGenerators);
    const set<unsigned>& tris = genItr->second;
    set<EdgeHash> meshEdges;
    bool boundedCell = true;
    for (set<unsigned>::const_iterator triItr = tris.begin();
         triItr != tris.end(); 
         ++triItr) {
      i = *triItr;
      ii = tri2id[i];
      POLY_ASSERT(i < numTriangles);
      POLY_ASSERT(tri2id.find(i) != tri2id.end());

      // Get the other two indices for thsi triangle and hash their edges wrt. p
      findOtherTriIndices(&triangleList[3*i], p, q, r);
      POLY_ASSERT(p != q and p != r);
      pq = internal::hashEdge(p,q);
      pr = internal::hashEdge(p,r);

      // Is pq a surface edge?
      if (edge2tri[pq].size() == 1) {
        boundedCell *= false;
        POLY_ASSERT(edge2tri[pq][0] == i);
        POLY_ASSERT(edge2id.find(pq) != edge2id.end());
        jj = edge2id[pq];
        if (jj != ii) meshEdges.insert(internal::hashEdge(ii,jj));
      }else{
        POLY_ASSERT((edge2tri[pq].size() == 2 and edge2tri[pq][0] == i)
                    or edge2tri[pq][1] == i);
        k = (edge2tri[pq][0] == i ? edge2tri[pq][1] : edge2tri[pq][0]);
        POLY_ASSERT(tri2id.find(k) != tri2id.end());
        jj = tri2id[k];
        if (jj != ii) meshEdges.insert(internal::hashEdge(ii,jj));
      }

      // Is pr a surface edge?
      if (edge2tri[pr].size() == 1){
        boundedCell *= false;
        POLY_ASSERT(edge2tri[pr][0] == i);
        POLY_ASSERT(edge2id.find(pr) != edge2id.end());
        jj = edge2id[pr];
        if (jj != ii) meshEdges.insert(internal::hashEdge(ii,jj));
      }else{
        POLY_ASSERT((edge2tri[pr].size() == 2 and edge2tri[pr][0] == i)
                     or edge2tri[pr][1] == i);
        k = (edge2tri[pr][0] == i ? edge2tri[pr][1] : edge2tri[pr][0]);
        POLY_ASSERT(tri2id.find(k) != tri2id.end());
        jj = tri2id[k];
        if (jj != ii) meshEdges.insert(internal::hashEdge(ii,jj));
      }
    }

    // Hook together the hashed edges around this point.
    cellNodes[p] = 
      computeSortedFaceNodes(vector<EdgeHash>(meshEdges.begin(), meshEdges.end()));

    // Check if the two inf-edges intersect each other.
    if (not boundedCell) {
      const unsigned last = cellNodes[p].size()-1;
    
      // Edge1
      const unsigned n11 = cellNodes[p][0];
      const unsigned n12 = cellNodes[p][last];
      POLY_ASSERT(n11 != n12);
      POLY_ASSERT(id2node.find(n11) != id2node.end());
      POLY_ASSERT(id2node.find(n12) != id2node.end());
      RealPoint rp11 = mCoords.dequantize(&(id2node[n11]).x);
      RealPoint rp12 = mCoords.dequantize(&(id2node[n12]).x);

      // Edge 2
      const unsigned n21 = cellNodes[p][2];
      const unsigned n22 = cellNodes[p][1];     // Projected
      POLY_ASSERT(n21 != n22);
      POLY_ASSERT(id2node.find(n21) != id2node.end());
      POLY_ASSERT(id2node.find(n22) != id2node.end());
      RealPoint rp21 = mCoords.dequantize(&(id2node[n21]).x);
      RealPoint rp22 = mCoords.dequantize(&(id2node[n22]).x);

      // Compute cell self-intersections
      RealPoint result;
      bool intersects;
      if (n12 == n21) intersects = false;
      else            intersects = geometry::segmentIntersection2D(&rp11.x, &rp12.x,
                                                                   &rp21.x, &rp22.x,
                                                                   &result.x,
                                                                   mDegeneracy);
      
      // Collapse the two projected nodes to the intersection point
      if (intersects) {
        // // Blago!
        // cerr << "Self-intersection detected:" << endl
        //      << "  Cell " << p << endl
	//      << "  " << rp11 << "  " << rp12 << endl
	//      << "  " << rp21 << "  " << rp22 << endl
	//      << "  Intersection pt = " << result << endl;
        // // Blago!
        
        if (geometry::distance<2,RealType>(&rp12.x, &result.x) > mDegeneracy and 
            geometry::distance<2,RealType>(&rp21.x, &result.x) > mDegeneracy) {
          const IntPoint intResult = mCoords.quantize(&result.x);
          id2node[n11] = intResult;
          id2node[n22] = intResult;
        }
      }
    }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
constructBoundedTopology(const vector<RealType>& points,
                         const ReducedPLC<2, RealType>& geometry,
                         const vector<ReducedPLC<2, CoordHash> >& intCells,
                         Tessellation<2, RealType>& mesh) const {
  const unsigned numGenerators = points.size()/2;
  int i, j, k, iedge;
  // std::map<IntPoint, int, ThreeByThreeCompare<CoordHash> > point2node;
  std::map<IntPoint, int> point2node;
  std::map<EdgeHash, int> edgeHash2id;
  std::map<int, std::vector<int> > edgeCells;
  mesh.cells = std::vector<std::vector<int> >(numGenerators);
  for (i = 0; i != numGenerators; ++i) { 

    // POLY_ASSERT(checkCellOrientation<CoordHash>
    //             (&(mCoords.quantize(&points[2*i])).x, intCells[i]));

    const unsigned nfacets = intCells[i].facets.size();
    POLY_ASSERT(nfacets > 2);
    for (unsigned ifacet = 0; ifacet < nfacets; ++ifacet) {
      POLY_ASSERT(intCells[i].facets[ifacet].size() == 2);
      const int i1 = intCells[i].facets[ifacet][0];
      const int i2 = intCells[i].facets[ifacet][1];
      POLY_ASSERT(i1 != i2);
      const IntPoint pX1 = IntPoint(intCells[i].points[2*i1], intCells[i].points[2*i1+1]);
      const IntPoint pX2 = IntPoint(intCells[i].points[2*i2], intCells[i].points[2*i2+1]);
      POLY_ASSERT(pX1 != pX2);
      j = internal::addKeyToMap(pX1, point2node);
      k = internal::addKeyToMap(pX2, point2node);
      POLY_ASSERT(j != k);
      iedge = internal::addKeyToMap(internal::hashEdge(j, k), edgeHash2id);
      edgeCells[iedge].push_back(j < k ? i : ~i);
      POLY_ASSERT2(edgeCells[iedge].size() == 1 or edgeCells[iedge][0]*edgeCells[iedge][1] <= 0,
                   "BLAGO: " << iedge << " " << j << " " << k << " " << edgeCells[iedge][0] 
		   << " " << edgeCells[iedge][1]);
      mesh.cells[i].push_back(j < k ? iedge : ~iedge);
    }
    POLY_ASSERT(mesh.cells[i].size() >= 3);
  }
  POLY_ASSERT(edgeCells.size() == edgeHash2id.size());

  // Fill in the mesh nodes
  RealPoint node;
  mesh.nodes = std::vector<RealType>(2*point2node.size());
  for (typename std::map<IntPoint, int>::const_iterator itr = point2node.begin();
       itr != point2node.end(); 
       ++itr) {
    const IntPoint& p = itr->first;
    i = itr->second;
    POLY_ASSERT(i < mesh.nodes.size()/2);
    node = mCoords.dequantize(&p.x);
    mesh.nodes[2*i]   = node.x;
    mesh.nodes[2*i+1] = node.y;
  }

  // Fill in the mesh faces.
  mesh.faces = std::vector<std::vector<unsigned> >(edgeHash2id.size());
  for (typename std::map<EdgeHash, int>::const_iterator itr = edgeHash2id.begin();
       itr != edgeHash2id.end(); ++itr) {
    const EdgeHash ehash = itr->first;
    i = itr->second;
    POLY_ASSERT(i < mesh.faces.size());
    POLY_ASSERT(mesh.faces[i].size() == 0);
    mesh.faces[i].push_back(ehash.first);
    mesh.faces[i].push_back(ehash.second);
  }

  // Fill in the mesh faceCells.
  mesh.faceCells = std::vector<std::vector<int> >(edgeHash2id.size());
  for (i = 0; i != mesh.faces.size(); ++i) {
    if (not(edgeCells[i].size() == 1 or edgeCells[i].size() == 2)) {
      const int n1 = mesh.faces[i][0], n2 = mesh.faces[i][1];
      std::cerr << "Blago! " << i << " " << edgeCells[i].size() << " : " << n1 << " " << n2 << " : ("
                << mesh.nodes[2*n1] << " " << mesh.nodes[2*n1 + 1] << ") ("
                << mesh.nodes[2*n2] << " " << mesh.nodes[2*n2 + 1] << ")" << std::endl;
      for (j = 0; j != edgeCells[i].size(); ++j) {
        std::cerr << " --> " << edgeCells[i][j] << " " 
                  << points[2*edgeCells[i][j]] << " " 
                  << points[2*edgeCells[i][j]+1] << std::endl;
      }
    }
    POLY_ASSERT(edgeCells[i].size() == 1 or edgeCells[i].size() == 2);
    mesh.faceCells[i] = edgeCells[i];
  }

  // Figure out which nodes are on the boundary
  std::set<unsigned> boundaryNodes;
  for (unsigned iface = 0; iface < mesh.faces.size(); ++iface) {
    POLY_ASSERT(mesh.faceCells[iface].size() == 1 or
                mesh.faceCells[iface].size() == 2 );
    if (mesh.faceCells[iface].size() == 1) {
      boundaryNodes.insert(mesh.faces[iface].begin(), mesh.faces[iface].end());
    }
  }
  POLY_ASSERT(boundaryNodes.size() > 0);
  POLY_ASSERT(boundaryNodes.size() <= mesh.nodes.size()/2);
  
  // Snap a subset of the boundary nodes to the PLC boundary locations
  std::set<unsigned> plcNodes;
  std::set<unsigned> indices;
  const RealType tol = 4.0*mCoords.delta;
  for (unsigned ii = 0; ii < geometry.points.size()/2; ++ii) indices.insert(ii);
  for (std::set<unsigned>::iterator iitr = indices.begin();
       iitr != indices.end();
       ++iitr) {
    unsigned ipoint = *iitr;
    unsigned inode;
    const RealPoint rp = RealPoint(geometry.points[2*ipoint], geometry.points[2*ipoint+1]);
    RealType dist = std::numeric_limits<RealType>::max();
    std::set<unsigned>::iterator nodeItr = boundaryNodes.begin();
    while (nodeItr != boundaryNodes.end() and dist > tol) {
      inode = *nodeItr;
      dist = geometry::distance<2, RealType>(&mesh.nodes[2*(*nodeItr)], &rp.x);
      ++nodeItr;
    }
    if (dist < tol) {
      plcNodes.insert(inode);
      mesh.nodes[2*inode  ] = rp.x;
      mesh.nodes[2*inode+1] = rp.y;
      boundaryNodes.erase(inode);
      indices.erase(ipoint);
    }
  }
  // // Blago!
  // if (not indices.empty()) {
  //   for (std::set<unsigned>::iterator itr = indices.begin();
  //        itr != indices.end();
  //        ++itr) {
  //      RealType dist = std::numeric_limits<RealType>::max();
  //      unsigned inode;
  //      for (unsigned ii = 0; ii < mesh.nodes.size()/2; ++ii) {
  //        RealType dnew = geometry::distance<2,RealType>(&geometry.points[2*(*itr)], &mesh.nodes[2*ii]);
  //        if (dnew < dist) {
  //          inode = ii;
  //          dist = dnew;
  //        }
  //      }
  //      cerr << "PLC point " << *itr << endl
  //           << "  location = (" << geometry.points[2*(*itr)  ] << "," << geometry.points[2*(*itr)+1] << ")" << endl
  //           << "  nearNode = " << inode << endl
  //           << "  node loc = (" << mesh.nodes[2*inode] << "," << mesh.nodes[2*inode+1] << ")" << endl
  //           << "  distance = " << dist << endl
  //           << "  boundaryNode? " << (boundaryNodes.find(inode) != boundaryNodes.end() ? "YES" : "NO") << endl;
  //   }
  // }
  // // Blago!
  POLY_ASSERT(indices.empty());
  POLY_ASSERT2(plcNodes.size() == geometry.points.size()/2,
               plcNodes.size() << " != " << geometry.points.size()/2);

  // Snap all boundary nodes to the nearest floating point along boundary
  for (std::set<unsigned>::iterator nodeItr = boundaryNodes.begin();
       nodeItr != boundaryNodes.end();
       ++nodeItr) {
    RealPoint result;
    RealType dist = nearestPoint(&mesh.nodes[2*(*nodeItr)],
                                 geometry.points.size()/2,
                                 &geometry.points[0],
                                 geometry,
                                 &result.x);
    if (dist >= 4.0*mCoords.delta) {
      std::cerr << "Blago!" << std::endl << "   Boundary node " 
                << *nodeItr << " at "
                << "(" << mesh.nodes[2*(*nodeItr)  ] 
                << "," << mesh.nodes[2*(*nodeItr)+1] << ")"
                << " wants to move distance " << dist << " to "
                << "(" << result[0] << "," << result[1] << ")" << std::endl;
    }
    else {
      mesh.nodes[2*(*nodeItr)  ] = result.x;
      mesh.nodes[2*(*nodeItr)+1] = result.y;
    }
    POLY_ASSERT(dist < 4.0*mCoords.delta);
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeDelaunay(const vector<RealType>& points,
                triangulateio& delaunay,
                const bool addStabilizingPoints) const {
  triangulateio in;
   
  // Find the range of the generator points.
  const unsigned numGenerators = points.size()/2;

  if (addStabilizingPoints) 
  {
    // Determine bounding box for points
    POLY_ASSERT(!mCoords.empty());
    RealPoint low  = RealPoint(mCoords.low [0], mCoords.low [1]);
    RealPoint high = RealPoint(mCoords.high[0], mCoords.high[1]);

    RealPoint box = high - low;
    const RealType boxsize = 3.0*max(box.x, box.y);
 
    const RealType xmin = 0.5*(low[0] + high[0]) - boxsize;
    const RealType xmax = 0.5*(low[0] + high[0]) + boxsize;
    const RealType ymin = 0.5*(low[1] + high[1]) - boxsize;
    const RealType ymax = 0.5*(low[1] + high[1]) + boxsize;

    mCoords.low [0] = min(mCoords.low [0], low [0]);
    mCoords.low [1] = min(mCoords.low [1], low [1]);
    mCoords.high[0] = max(mCoords.high[0], high[0]);
    mCoords.high[1] = max(mCoords.high[1], high[1]);

    // Add the generators
    in.numberofpoints = numGenerators + 4;
    in.pointlist = new RealType[2*in.numberofpoints];
    copy(points.begin(), points.end(), in.pointlist);
    in.pointlist[2*numGenerators  ] = xmin;  in.pointlist[2*numGenerators+1] = ymin;
    in.pointlist[2*numGenerators+2] = xmax;  in.pointlist[2*numGenerators+3] = ymin;
    in.pointlist[2*numGenerators+4] = xmax;  in.pointlist[2*numGenerators+5] = ymax;
    in.pointlist[2*numGenerators+6] = xmin;  in.pointlist[2*numGenerators+7] = ymax;
  } 
  else 
  {
    in.numberofpoints = numGenerators;
    in.pointlist = new RealType[2*in.numberofpoints];
    copy(points.begin(), points.end(), in.pointlist);
  }
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


//------------------------------------------------------------------------------
// Private tessellate routines
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const ReducedPLC<2, CoordHash>& intGeometry,
           const QuantizedCoordinates<2, RealType>& coords,
           vector<ReducedPLC<2, CoordHash> >& intCells) const {
  // Pre-conditions
  POLY_ASSERT(not intGeometry.empty());
  POLY_ASSERT(points.size() > 0);
  POLY_ASSERT(points.size() % 2 == 0);
  POLY_ASSERT(not coords.empty());

  // The Quantized coordinates
  mCoords = coords;
  
  const unsigned numGenerators = points.size()/2;
  const bool collinear = geometry::collinear<2, RealType>(points, mDegeneracy);
  vector<vector<unsigned> > cellNodes;
  map<int, IntPoint> id2node;
  vector<unsigned> infNodes;
  
  // Use the appropriate cell node routine
  if (collinear) 
  {
    this->computeCellNodesCollinear(points, cellNodes, id2node, infNodes);
  }
  else 
  {
    this->computeCellNodes(points, cellNodes, id2node, infNodes);
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);
  
  vector<ReducedPLC<2, CoordHash> > dummy;
  intCells.resize(numGenerators);
  for (int i = 0; i < numGenerators; ++i) {
    const IntPoint intGenerator = coords.quantize(&points[2*i]);
    const ReducedPLC<2, CoordHash> intCell = plcOfCell<CoordHash>(cellNodes[i], id2node);
    intCells[i] = BG::boost_clip<CoordHash>(intGeometry, intCell, intGenerator, dummy);
    POLY_ASSERT(dummy.empty());
    POLY_ASSERT2(not BG::boost_intersects(intCells[i]), intCells[i]);
  }  
}
//------------------------------------------------------------------------------









//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TriangleTessellator<double>;

}
