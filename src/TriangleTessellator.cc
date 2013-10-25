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
#include "float.h"

#include "polytope.hh" // Pulls in POLY_ASSERT and TriangleTessellator.hh.

#include "convexHull_2d.hh"
#include "nearestPoint.hh"
#include "within.hh"
#include "intersect.hh"

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

// //------------------------------------------------------------------------------
// // Compute the intersection of a line segment and a PLC. Returns number of
// // intersections, intersected facet(s) along PLC, and intersection point(s)
// //------------------------------------------------------------------------------
// template<typename RealType>
// inline
// unsigned intersectBoundingBox(const RealType* point1,
// 			      const RealType* point2,
// 			      const unsigned numVertices,
// 			      const RealType* vertices,
// 			      const std::vector<std::vector<int> >& facets,
// 			      std::vector<int>& facetIntersections,
// 			      std::vector<RealType>& result) {
//   unsigned i, j, numIntersections=0;
//   RealType intersectionPoint[2];
//   bool intersects, addPoint;
//   const unsigned numFacets = facets.size();
//   for (unsigned ifacet = 0; ifacet != numFacets; ++ifacet) {
//     POLY_ASSERT(facets[ifacet].size() == 2);
//     i = facets[ifacet][0];
//     j = facets[ifacet][1];
//     POLY_ASSERT(i >= 0 and i < numVertices);
//     POLY_ASSERT(j >= 0 and j < numVertices);
//     intersects = geometry::segmentIntersection2D(point1, point2,
// 						 &vertices[2*i], &vertices[2*j],
// 						 intersectionPoint);
//     addPoint = false;
//     if (intersects) {
//       if (numIntersections == 0) addPoint = true;
//       else {
//         if (intersectionPoint[0] != result.back() - 1 and
//             intersectionPoint[1] != result.back() ) addPoint = true;
//         else addPoint = false;
//       }
//     }
//      
//     if (addPoint) {
//       RealType p[2];
//       p[(ifacet+1)%2] = vertices[2*ifacet + (ifacet+1)%2];
//       p[ ifacet   %2] = intersectionPoint[ifacet%2];
//       ++numIntersections;
//       result.push_back(p[0]);
//       result.push_back(p[1]);
//       facetIntersections.push_back(ifacet);
//     }
//   }
//   POLY_ASSERT(numIntersections == result.size()/2);
//   return numIntersections;
// }
  

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
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(!points.empty());
  POLY_ASSERT(points.size() % 2 == 0);
  POLY_ASSERT(points.size() > 2 );
  
  mCoords.initialize(points);

  this->computeVoronoiUnbounded(points, mesh);
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
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(!points.empty());
  POLY_ASSERT(points.size() % 2 == 0);
  
  // // BLAGO!
  // std::ofstream dumpfile;
  // dumpfile.open("generators.txt");
  // const unsigned n = points.size()/2;
  // std::cerr << "HERE WE GO : " << n << std::endl;
  // for (unsigned i = 0; i != n; ++i) {
  //   dumpfile << points[2*i] << " " << points[2*i+1] << std::endl;
  // }
  // dumpfile.close();
  // // BLAGO!

  // Build a PLC with the bounding box, and then use the PLC method.
  ReducedPLC<2, RealType> box = this->boundingBox(low, high);
  this->tessellate(points, box.points, box, mesh);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const vector<RealType>& PLCpoints,
           const PLC<2, RealType>& geometry,
           Tessellation<2, RealType>& mesh) const {
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(!points.empty() and !PLCpoints.empty());
  POLY_ASSERT(points.size() % 2 == 0 and PLCpoints.size() % 2 == 0);
  POLY_ASSERT(!geometry.facets.empty());
  
  mCoords.initialize(points);
  
  this->computeVoronoiBounded(points, PLCpoints, geometry, mesh);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Private routines called by tessellate:
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeCellNodes(const vector<RealType>& points,
		 vector<RealPoint>& nodeList,
		 vector<vector<unsigned> >& cellNodes,
		 vector<unsigned>& infNodes) const {
  POLY_ASSERT(!points.empty());
  POLY_ASSERT(points.size() != 2);

  const unsigned numGenerators = points.size()/2;
  cellNodes.resize(numGenerators);
  
  // Compute the triangularization
  triangulateio delaunay;
  computeDelaunay(points, delaunay);

  //--------------------------------------------------------
  // Create the Voronoi tessellation from the triangulation.
  //--------------------------------------------------------
  
  // Create the Voronoi nodes from the list of triangles. Each triangle 
  // has 3 nodes p, q, r, and corresponds to a Voronoi node at (X,Y), say.

  // Find the circumcenters of each triangle, and build the set of triangles
  // associated with each generator.
  vector<RealPoint> circumcenters(delaunay.numberoftriangles);
  vector<unsigned> triMask(delaunay.numberoftriangles, 0);
  map<EdgeHash, vector<unsigned> > edge2tri;
  map<int, set<unsigned> > gen2tri;
  int pindex, qindex, rindex, i, j, k, triCount = 0;
  EdgeHash pq, pr, qr;
 
  RealType lowc [2] = {mCoords.low [0], mCoords.low [1]};
  RealType highc[2] = {mCoords.high[0], mCoords.high[1]};
  for (i = 0; i != delaunay.numberoftriangles; ++i) {
    pindex = delaunay.trianglelist[3*i  ];
    qindex = delaunay.trianglelist[3*i+1];
    rindex = delaunay.trianglelist[3*i+2];
    pq = internal::hashEdge(pindex, qindex);
    pr = internal::hashEdge(pindex, rindex);
    qr = internal::hashEdge(qindex, rindex);
    geometry::computeCircumcenter(&delaunay.pointlist[2*pindex],
				  &delaunay.pointlist[2*qindex],
				  &delaunay.pointlist[2*rindex],
				  &circumcenters[i].x);
    if (orient2d(&delaunay.pointlist[2*pindex],
		 &delaunay.pointlist[2*qindex],
		 &delaunay.pointlist[2*rindex]) == 0) {
      cerr << "BLAGO!" << endl
	   << delaunay.pointlist[2*pindex  ] << "," 
	   << delaunay.pointlist[2*pindex+1] << endl
	   << delaunay.pointlist[2*qindex  ] << "," 
	   << delaunay.pointlist[2*qindex+1] << endl
	   << delaunay.pointlist[2*rindex  ] << "," 
	   << delaunay.pointlist[2*rindex+1] << endl;
	}
    POLY_ASSERT(orient2d(&delaunay.pointlist[2*pindex],
                         &delaunay.pointlist[2*qindex],
                         &delaunay.pointlist[2*rindex]) != 0);
    if (pindex < numGenerators and qindex < numGenerators and rindex < numGenerators) {
      triMask[i] = 1;
      gen2tri[pindex].insert(i);
      gen2tri[qindex].insert(i);
      gen2tri[rindex].insert(i);
      edge2tri[pq].push_back(i);
      edge2tri[pr].push_back(i);
      edge2tri[qr].push_back(i);
      lowc [0] = min(lowc [0], circumcenters[i].x);
      lowc [1] = min(lowc [1], circumcenters[i].y);
      highc[0] = max(highc[0], circumcenters[i].x);
      highc[1] = max(highc[1], circumcenters[i].y);
      ++triCount;
    }
  }
  POLY_ASSERT(circumcenters.size() == delaunay.numberoftriangles);
  POLY_ASSERT(triMask.size()       == delaunay.numberoftriangles);
  POLY_ASSERT(std::accumulate(triMask.begin(), triMask.end(), 0) > 0);
  
  POLY_ASSERT(lowc[0] <= highc[0] and lowc[1] <= highc[1]);
  mCoords.expand(lowc, highc);

  // Determine which circumcenters lie inside the inner bounding box
  // Map circumcenters and triangle indices to global id's
  int inside;
  IntPoint ip;
  map<IntPoint, int> circ2id;
  map<int, unsigned> tri2id;
  map<IntPoint, pair<int,int> > nodeData;
  for (i = 0; i != delaunay.numberoftriangles; ++i){
    if (triMask[i] == 1) {
      ip = mCoords.quantize(&circumcenters[i].x);
      inside = (mCoords.isInside(&circumcenters[i].x)) ? 1 : 0;
      k = circ2id.size();
      j = internal::addKeyToMap(ip, circ2id);
      tri2id[i] = j;
      if (k != circ2id.size()) nodeData[ip] = make_pair(j,inside);
    }
  }
  POLY_ASSERT(circ2id.size() == nodeData.size());

  // The exterior edges of the triangularization have "unbounded" rays, originating
  // at the circumcenter of the corresponding triangle and passing perpendicular to
  // the edge
  RealPoint ehat, pinf;
  map<EdgeHash, unsigned> edge2id;
  int i1, i2, ivert;
  infNodes = vector<unsigned>(circ2id.size());
  for (map<EdgeHash, vector<unsigned> >::const_iterator edgeItr = edge2tri.begin();
       edgeItr != edge2tri.end(); ++edgeItr){
    const EdgeHash& edge = edgeItr->first;
    const vector<unsigned>& tris = edgeItr->second;
    if (tris.size() == 1){
      i = tris[0];
      POLY_ASSERT(i < delaunay.numberoftriangles);
      i1 = edge.first;
      i2 = edge.second;
      findOtherTriIndex(&delaunay.trianglelist[3*i], i1, i2, ivert);
      computeEdgeUnitVector(&delaunay.pointlist[2*i1],
			    &delaunay.pointlist[2*i2],
			    &delaunay.pointlist[2*ivert],
			    &ehat.x);
      
      pinf = mCoords.projectPoint(&circumcenters[i].x, &ehat.x);
      IntPoint ip = mCoords.quantize(&pinf.x);
      int inside = 0;
      //pinf = mCoords.projectPoint(&circumcenters[i].x, &ehat.x);
      //IntPoint ip = mCoords.quantize(&pinf.x);
      //int inside = 1;
      
      k = circ2id.size();
      j = internal::addKeyToMap(ip, circ2id);
      if (k != circ2id.size())  nodeData[ip] = make_pair(j,inside);
      POLY_ASSERT(edge2id.find(edge) == edge2id.end());
      edge2id[edge] = j;
      if (k != circ2id.size()) infNodes.push_back(1);

      // //Blago!
      // cerr << endl;
      // cerr << "Boundary edge          = (" << edge.first << "," << edge.second << ")" << endl;
      // cerr << "Neighboring Triangle   = " << i << endl;
      // cerr << "Third triangle vertex  = " << ivert << endl;
      // cerr << "Real projected infNode = " << pinf << endl;
      // cerr << "Int projected infNode  = " << ip << endl;
      // cerr << "j index                = " << j << endl;
      // cerr << "Circumcenter box       = " << cbox[0] << " " << cbox[1] << endl;
      // //Blago!
    }
  }
  POLY_ASSERT(circ2id.size() == nodeData.size());

  // Store the final collection of nodes in an ordered vector
  nodeList.resize(nodeData.size());
  for (typename map<IntPoint, pair<int,int> >::const_iterator itr = nodeData.begin();
       itr != nodeData.end();
       ++itr) {
    i = itr->second.first;
    POLY_ASSERT(i < nodeData.size());
    IntPoint ip  = itr->first;
    RealPoint rp = mCoords.dequantize(&ip.x, (itr->second.second==1));
    nodeList[i] = rp;
  }
  
  // The faces corresponding to each triangle edge
  unsigned ii, jj;
  bool unboundedCell;
  for (map<int, set<unsigned> >::const_iterator genItr = gen2tri.begin();
       genItr != gen2tri.end(); 
       ++genItr) {
    pindex = genItr->first;
    const set<unsigned>& tris = genItr->second;
    POLY_ASSERT(pindex < numGenerators);
    
    set<EdgeHash> meshEdges;
    unboundedCell = false;
    for (set<unsigned>::const_iterator triItr = tris.begin();
         triItr != tris.end(); ++triItr){
      i = *triItr;
      POLY_ASSERT(i < delaunay.numberoftriangles);
      POLY_ASSERT(tri2id.find(i) != tri2id.end());
      ii = tri2id[i];
      
      // Get the other indices for this triangle, given one of its vertices pindex
      findOtherTriIndices(&delaunay.trianglelist[3*i], pindex, qindex, rindex);
      pq = internal::hashEdge(pindex,qindex);
      pr = internal::hashEdge(pindex,rindex);
      
      // Is pq a surface edge?
      if (edge2tri[pq].size() == 1){
        POLY_ASSERT(edge2tri[pq][0] == i);
        POLY_ASSERT(edge2id.find(pq) != edge2id.end());
        jj = edge2id[pq];
        if (jj != ii) meshEdges.insert(internal::hashEdge(ii,jj));
      }else{
         POLY_ASSERT((edge2tri[pq].size() == 2 and edge2tri[pq][0] == i)
                     or edge2tri[pq][1] == i);
        k = (edge2tri[pq][0] == i ? edge2tri[pq][1] : edge2tri[pq][0]);
        jj = tri2id[k];
        if (jj != ii) meshEdges.insert(internal::hashEdge(ii,jj));
      }
      
      // Is pr a surface edge?
      if (edge2tri[pr].size() == 1){
        POLY_ASSERT(edge2tri[pr][0] == i);
        POLY_ASSERT(edge2id.find(pr) != edge2id.end());
        jj = edge2id[pr];
        if (jj != ii) meshEdges.insert(internal::hashEdge(ii,jj));
      }else{
         POLY_ASSERT((edge2tri[pr].size() == 2 and edge2tri[pr][0] == i)
                     or edge2tri[pr][1] == i);
        k = (edge2tri[pr][0] == i ? edge2tri[pr][1] : edge2tri[pr][0]);
        jj = tri2id[k];
        if (jj != ii) meshEdges.insert(internal::hashEdge(ii,jj));
      }

      // If either of pq or pr has one triangle neighbor,
      // this is an unbounded cell
      unboundedCell += (edge2tri[pq].size() == 1 or edge2tri[pr].size() == 1);
    }

    // // Blago!
    // cerr << "\n\nBLAGO!\n" << pindex << endl;
    // for (set<EdgeHash>::iterator itr = meshEdges.begin();  
    //      itr != meshEdges.end(); ++itr )
    //    cerr << "(" << itr->first << "," << itr->second << ")" << endl;
    // // Blago!
    
    cellNodes[pindex] = 
      computeSortedFaceNodes(vector<EdgeHash>(meshEdges.begin(), meshEdges.end()));

    // Check for self intersections in the sorted cell nodes.
    // By convention, the chain of nodes for an unbounded cell
    // begins and ends at projected inf nodes. This is only
    // meant to correct self intersections due to projection
    if (unboundedCell) {
      unsigned last = cellNodes[pindex].size()-1;

      // Edge 1
      const unsigned n11 = cellNodes[pindex][0];     // Projected
      const unsigned n12 = cellNodes[pindex][last];
      POLY_ASSERT(n11 != n12);
      POLY_ASSERT(n11 < nodeList.size());
      POLY_ASSERT(n12 < nodeList.size());
      RealPoint rp11 = nodeList[n11], rp12 = nodeList[n12];
      
      // Edge 2
      const unsigned n21 = cellNodes[pindex][2];
      const unsigned n22 = cellNodes[pindex][1];     // Projected
      POLY_ASSERT(n21 != n22);
      POLY_ASSERT(n21 < nodeList.size());
      POLY_ASSERT(n22 < nodeList.size());
      RealPoint rp21 = nodeList[n21], rp22 = nodeList[n22];

      // Compute cell self-intersections
      RealPoint result;
      bool intersects;
      if (n12 == n21) intersects = false;
      else            intersects = geometry::segmentIntersection2D(&rp11.x, &rp12.x,
                                                                   &rp21.x, &rp22.x,
                                                                   &result.x);

      // Self-intersecting cell ring: replace the offending projected
      // points with the intersection point
      if (intersects) {
	
	// Blago!
	cerr << "Self-intersection detected:" << endl
	     << "  Cell " << pindex << endl
	     << "  " << rp11 << "  " << rp12 << endl
	     << "  " << rp21 << "  " << rp22 << endl
	     << "  Intersection pt = " << result << endl;
	// Blago!

	nodeList[n11] = RealPoint(0.95*result.x + 0.05*rp12.x,
				  0.95*result.y + 0.05*rp12.y);
	nodeList[n22] = RealPoint(0.95*result.x + 0.05*rp21.x,
				  0.95*result.y + 0.05*rp21.y);
      }

      
      // Compute intersection with bounding box
      vector<RealType> results;
      unsigned nintersect = intersect(&rp11.x, &rp22.x, 
                                      mCoords.points.size()/2,
                                      &mCoords.points[0],
                                      mCoords, results);
      
      // Node chain endpoints connect through the inner bounding box:
      // Project an additional point outward
      if (nintersect > 0) {
        RealPoint p = RealPoint(points[2*pindex], points[2*pindex+1]);
        RealPoint r = RealPoint(rp22.x-rp11.x, rp22.y-rp11.y);
        RealPoint rperp = (orient2d(&rp11.x, &rp12.x, &rp22.x) > 0.0) ? 
           RealPoint(-r.y, r.x) : RealPoint(r.y, -r.x);
        geometry::unitVector<2, RealType>(&rperp.x);
        RealPoint newNode = mCoords.projectPoint(&p.x, &rperp.x);
        unsigned newNodeIndex = nodeList.size();

         // Blago!
	cerr << "Bounding box intersection detected:" << endl
	     << "  Cell " << pindex << " at " << p << endl
	     << "  " << rp11 << "  " << rp12 << endl
             << "  " << rp21 << "  " << rp22 << endl
             << "  New node = " << newNode << endl;
        for (unsigned ii = 0; ii != nintersect; ++ii) {
           cerr << "  Intersect Pt " << ii << ": "
                << results[2*ii] << " " << results[2*ii+1] << endl;
        }
        cerr << endl;
	// Blago!

        nodeList.push_back(newNode);
        vector<unsigned>::iterator it = cellNodes[pindex].begin() + 1;
        cellNodes[pindex].insert(it, newNodeIndex);
        // cellNodes[pindex].push_back(newNodeIndex);
      }
    }
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);

  // Post conditions
  POLY_ASSERT(cellNodes.size() == numGenerators); 
  POLY_ASSERT(!nodeList.empty());

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
computeCellRings(const vector<RealType>& points,
		 const vector<RealPoint>& nodeList,
		 vector<vector<unsigned> >& cellNodes,
		 Clipper2d<CoordHash>& clipper,
                 vector<BGring>& cellRings) const {
  const unsigned numGenerators = points.size()/2;
  int i;
  
  // Circumcenters that lie outside the bounding box of the PLC boundary 
  // are quantized based on different criteria to avoid contaminating the 
  // degeneracy spacing of the mesh nodes. We will project these outer 
  // circumcenters to the edges of the bounding box so all nodes follow 
  // the input degeneracy spacing.

  // Quantize the PLCpoints
  vector<RealPoint> boundingBoxPoints(4);
  boundingBoxPoints[0] = RealPoint(mCoords.low [0], mCoords.low [1]);
  boundingBoxPoints[1] = RealPoint(mCoords.high[0], mCoords.low [1]);
  boundingBoxPoints[2] = RealPoint(mCoords.high[0], mCoords.high[1]);
  boundingBoxPoints[3] = RealPoint(mCoords.low [0], mCoords.high[1]);
  
  RealPolygon boundingBox;
  RealRing boundingBoxRing = RealRing(boundingBoxPoints.begin(),
				      boundingBoxPoints.end());
  boost::geometry::correct(boundingBoxRing);
  boost::geometry::assign(boundingBox, boundingBoxRing);

  // Initialize the object to handle cell intersections
  Clipper2d<RealType> boundingBoxClipper(boundingBox);

  // Walk the nodes around each generator and build up the cell ring
  IntPoint X, ip1, ip2;
  RealPoint rp1, rp2;
  vector<RealRing> boundingBoxOrphans;
  vector<BGring> orphans;
  cellRings.resize(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
    // Removed this because inner/outer bounding box intersections
    // are now done with boost.geometry. The ring concept enforces
    // CCW orientation without having to reverse it manually. 09/19/2013

    // // // Check the orientation of the node list and reverse it if it's CW
    // //  {
    // //    POLY_ASSERT(cellNodes[i].size() > 1);
    // //    unsigned i1 = cellNodes[i][0];
    // //    unsigned i2 = cellNodes[i][1];
    // //    POLY_ASSERT(i1 < numNodes and i2 < numNodes);
    // //    rp1 = nodeList[i1];
    // //    rp2 = nodeList[i2];
    // //    if (orient2d(&rp1.x, &rp2.x, (double*)&points[2*i]) < 0) {
    // //      reverse(cellNodes[i].begin(), cellNodes[i].end());
    // //    }
    // //  }

    // Blago!
    bool Blago = false;
    if(Blago) cerr << "---------- Cell " << i << " -----------" << endl;
    // Blago!

    // Create a BG ring and fill in any ring concept stuff that's missing
    vector<RealPoint> realCellNodes;
    for (vector<unsigned>::const_iterator itr = cellNodes[i].begin();
	 itr != cellNodes[i].end(); 
	 ++itr) {
      rp1 = nodeList[*itr];
      realCellNodes.push_back(rp1);
    }
    POLY_ASSERT(realCellNodes.size() == cellNodes[i].size());
    RealRing realCellRing(realCellNodes.begin(), realCellNodes.end());
    boost::geometry::correct(realCellRing);

    // Blago!
    if(Blago) {
    cerr << "\n\nPre-Clip" << endl;
    for (typename RealRing::const_iterator itr = realCellRing.begin();
	 itr != realCellRing.end();
	 ++itr ) cerr << *itr << endl;
    }
    // Blago!

    // Clip real cell against inner bounding box
    /////POLY_ASSERT(!boost::geometry::intersects(realCellRing));
    boundingBoxClipper.clipCell(RealPoint(points[2*i], points[2*i+1]),
				realCellRing,
				boundingBoxOrphans);

    // Blago!
    if(Blago) {
    cerr << "\n\nPost-Clip" << endl;
    for (typename RealRing::const_iterator itr = realCellRing.begin();
	 itr != realCellRing.end();
	 ++itr ) cerr << *itr << endl;
    }
    // Blago!


    // Remove any repeated points
    boost::geometry::unique(realCellRing);
    POLY_ASSERT(!realCellRing.empty());
    POLY_ASSERT(realCellRing.front() == realCellRing.back());
    POLY_ASSERT(boundingBoxOrphans.empty());

    // Convert ring clipped by bounding box to integers
    vector<IntPoint> cellBoundary;
    for (typename RealRing::iterator itr = realCellRing.begin();
	 itr != realCellRing.end(); ++itr) {
      POLY_ASSERT2((*itr).x >= mCoords.low [0] and
		   (*itr).x <= mCoords.high[0],
		   *itr << " " << mCoords.low[0] << " " << mCoords.high[0]);
      POLY_ASSERT2((*itr).y >= mCoords.low [1] and
		   (*itr).y <= mCoords.high[1], 
		   *itr << " " << mCoords.low[1] << " " << mCoords.high[1]);
      cellBoundary.push_back(mCoords.quantize(&(*itr).x));
    }
    boost::geometry::assign(cellRings[i],
			    BGring(cellBoundary.begin(),
				   cellBoundary.end()));
    boost::geometry::correct(cellRings[i]);
    POLY_ASSERT(!cellRings[i].empty());
    POLY_ASSERT(cellRings[i].front() == cellRings[i].back());

    // Compute the boundary intersections
    clipper.clipCell(mCoords.quantize(&points[2*i]), cellRings[i], orphans);
    
    // Remove any repeated points
    boost::geometry::unique(cellRings[i]);
    

    // Blago!
    if(Blago){
    cerr << endl << "Final clipped cell ring " << i << endl;
    for (typename BGring::iterator itr = cellRings[i].begin();
         itr != cellRings[i].end(); ++itr) {
      cerr << mCoords.dequantize(&(*itr).x) << endl;
    }
    }
    // Blago!

  }

  // // Blago!
  // for (i = 0; i != orphans.size(); ++i) {
  //   cerr << endl << "Orphan " << i << endl;
  //   for (typename BGring::iterator itr = orphans[i].begin();
  //        itr != orphans[i].end(); ++itr) {
  //     cerr << (*itr).realx(mLow[0], mDelta) << " "
  //          << (*itr).realy(mLow[1], mDelta) << endl;
  //   }
  // }
  // // Blago!

  // If any orphaned cells exist, run the adoption algorithm
  // and modify the neighboring cell rings
  if (!orphans.empty()) {
    BoostOrphanage<RealType> orphanage(this);
    orphanage.adoptOrphans(points, mCoords, cellRings, orphans);
  }
  
  // Mark any cell ring points coinciding with the bounding PLC
  for (i = 0; i != numGenerators; ++i) {
    if (boost::geometry::intersects(cellRings[i], clipper.mBoundary)) {
      for (typename BGring::iterator itr = cellRings[i].begin();
	   itr != cellRings[i].end()-1; ++itr) {
	for (typename BGring::const_iterator oItr = clipper.mBoundary.outer().begin();
	     oItr != clipper.mBoundary.outer().end();
	     ++oItr) {
	  if (*itr == *oItr)  itr->index = 1;
	}
	vector<BGring>& holes = clipper.mBoundary.inners();
	for (unsigned ihole = 0; ihole != holes.size(); ++ihole) {
	  for (typename BGring::const_iterator iItr = holes[ihole].begin();
	       iItr != holes[ihole].end()-1;
	       ++iItr) {
	    if (*itr == *iItr)  itr->index = 1;
	  }
	}
	// //Blago!
	// if (itr->index == 1) std::cerr << (*itr) << std::endl;
	// //Blago!
      }
    }
  }


  
  // Post-conditions
  POLY_ASSERT(cellRings.size() == numGenerators);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeVoronoiUnbounded(const vector<RealType>& points,
			Tessellation<2, RealType>& mesh) const {
  POLY_ASSERT(!points.empty());
  POLY_ASSERT(points.size() != 2);
  
  const unsigned numGenerators = points.size()/2;
  vector<RealPoint> nodeList;
  vector<vector<unsigned> > cellNodes;
  int i;

  // Check for collinearity and use the appropriate routine
  bool collinear = true;
  if (numGenerators > 2 ) {
    i = 2;
    while (collinear and i != numGenerators) {
      collinear *= geometry::collinear<2,RealType>(&points[0], &points[2], 
                                                   &points[2*i], 1.0e-10);
      ++i;
    }
  }

  // Use the appropriate cell node routine
  if (collinear) {
    computeCellNodesCollinear(points, mCoords, nodeList, cellNodes);
    mesh.infNodes = vector<unsigned>(nodeList.size(), 1);
  } else {
    vector<unsigned> infNodes;
    this->computeCellNodes(points, nodeList, cellNodes, infNodes);
    mesh.infNodes = infNodes;
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);
  POLY_ASSERT(!nodeList.empty());

  // Copy the quantized nodes to the final tessellation.
  RealPoint node;
  const unsigned numNodes = nodeList.size();
  mesh.nodes.resize(2*numNodes);
  for (i = 0; i != nodeList.size(); ++i) {
    mesh.nodes[2*i  ] = nodeList[i].x;
    mesh.nodes[2*i+1] = nodeList[i].y;
  }
  POLY_ASSERT(mesh.infNodes.size() == mesh.nodes.size()/2);
  
  // Finish constructing the cell-face-node-topology
  constructUnboundedMeshTopology(cellNodes, points, mesh);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeVoronoiBounded(const vector<RealType>& points,
		      const vector<RealType>& PLCpoints,
		      const PLC<2, RealType>& geometry,
		      Tessellation<2, RealType>& mesh) const {
  POLY_ASSERT(!points.empty());
  POLY_ASSERT(mesh.empty());
  
  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints = PLCpoints.size()/2;
  vector<RealPoint> nodeList;
  vector<vector<unsigned> > cellNodes;
  int i;
  
  // Check for collinear generators
  bool collinear = true;
  if (numGenerators > 2 ) {
    i = 2;
    while (collinear and i != numGenerators) {
      collinear *= geometry::collinear<2,RealType>(&points[0], &points[2], 
                                                   &points[2*i], 1.0e-10);
      ++i;
    }
  }

  // Use the appropriate cell node routine
  if (collinear) {
    computeCellNodesCollinear(points, mCoords, nodeList, cellNodes);
  } else {
    vector<unsigned> infNodes;
    this->computeCellNodes(points, nodeList, cellNodes, infNodes);
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);
  POLY_ASSERT(!nodeList.empty());
  
  // Quantize the PLCpoints
  vector<IntPoint> IntPLCPoints(numPLCpoints);
  for (i = 0; i < numPLCpoints; ++i){
    IntPLCPoints[i] = mCoords.quantize(&PLCpoints[2*i]);
  }

  // Generate the quantized boundary to handle boost intersections
  BGpolygon boundary;
  constructBoostBoundary(IntPLCPoints, geometry, boundary);

  // Initialize the object to handle cell intersections
  Clipper2d<CoordHash> clipper(boundary);
  
  // Compute bounded cell rings
  vector<BGring> cellRings;
  this->computeCellRings(points, nodeList, cellNodes, clipper, cellRings);
  
  // Input nodes and construct the final mesh topology
  constructBoundedMeshTopology(cellRings, points, PLCpoints, geometry, mCoords, mesh);
}
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

  // Determine bounding box for points
  POLY_ASSERT(!mCoords.empty());
  RealType low [2] = {mCoords.low [0], mCoords.low [1]};
  RealType high[2] = {mCoords.high[0], mCoords.high[1]};

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


//------------------------------------------------------------------------------
// Private tessellate routines
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const vector<CoordHash>& IntPLCpoints,
           const PLC<2, RealType>& geometry,
           const QuantizedCoordinates<2, RealType>& coords,
           vector<vector<vector<CoordHash> > >& IntCells) const {
  // Pre-conditions
  POLY_ASSERT(!geometry.empty());
  POLY_ASSERT(points.size() > 0 and IntPLCpoints.size() > 0);
  POLY_ASSERT(points.size() % 2 == 0);
  POLY_ASSERT(!coords.empty());

  // The Quantized coordinates
  mCoords = coords;

  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints = IntPLCpoints.size()/2;
  vector<RealPoint> nodeList;
  vector<vector<unsigned> > cellNodes;
  int i;
  
  // Check for collinear generators
  bool collinear = true;
  if (numGenerators > 2 ) {
    i = 2;
    while (collinear and i != numGenerators) {
      collinear *= geometry::collinear<2,RealType>(&points[0], &points[2], 
                                                   &points[2*i], 1.0e-10);
      ++i;
    }
  }

  // Use the appropriate cell node routine
  if (collinear) {
    computeCellNodesCollinear(points, mCoords, nodeList, cellNodes);
  }
  else {
    vector<unsigned> infNodes;
    this->computeCellNodes(points, nodeList, cellNodes, infNodes);
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);
  POLY_ASSERT(!nodeList.empty());
  
  // Store the input boundary as a Boost.Geometry polygon
  BGpolygon boundary;
  vector<IntPoint> boundaryPoints;
  for (i = 0; i != numPLCpoints; ++i)
     boundaryPoints.push_back(IntPoint(IntPLCpoints[2*i], IntPLCpoints[2*i+1]));
  constructBoostBoundary(boundaryPoints, geometry, boundary);
  
  // Initialize the object to handle cell intersections
  Clipper2d<CoordHash> clipper(boundary);
  
  // Compute bounded cell rings
  vector<BGring> cellRings;
  this->computeCellRings(points, nodeList, cellNodes, clipper, cellRings);
  
  // Store the rings in a non-Boost way
  IntCells.resize(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
    vector<CoordHash> node(2);
    int index = 0;
    IntCells[i].resize(cellRings[i].size());
    for (typename BGring::const_iterator itr = cellRings[i].begin();
        itr != cellRings[i].end(); ++itr, ++index) {
      node[0] = (*itr).x;
      node[1] = (*itr).y;
      POLY_ASSERT(node.size() == 2);
      IntCells[i][index] = node;
    }
    POLY_ASSERT(IntCells[i].size()     == cellRings[i].size()  );
    POLY_ASSERT(IntCells[i].front()[0] == IntCells[i].back()[0]);
    POLY_ASSERT(IntCells[i].front()[1] == IntCells[i].back()[1]);
  }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template class TriangleTessellator<double>;
}




// //------------------------------------------------------------------------------
// // Predicate to check if either element of a std::pair matches the given value.
// //------------------------------------------------------------------------------
// template<typename T>
// struct MatchEitherPairValue {
//   T mval;
//   MatchEitherPairValue(const T& x): mval(x) {}
//   bool operator()(const std::pair<T, T>& x) const { return (x.first == mval or x.second == mval); }
// };

// //------------------------------------------------------------------------------
// // Get the set of map keys
// //------------------------------------------------------------------------------
// template<typename Key, typename T>
// std::set<Key>
// getMapKeys(std::map<Key, T>& mapIn) {
//   typename std::set<Key> result;
//   for (typename std::map<Key, T>::const_iterator itr = mapIn.begin();
//        itr != mapIn.end(); ++itr){
//     result.insert( itr->first );
//   }
//   return result;
// }



    /*
  // Walk the nodes around each generator and build up the cell ring
  IntPoint X, ip1, ip2;
  RealPoint rp1, rp2;
  unsigned i1, i2, nints;
  vector<BGring> orphans;
  cellRings.resize(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
    // Check the orientation of the node list and reverse it if it's CW
    POLY_ASSERT(cellNodes[i].size() > 1);
    i1 = cellNodes[i][0];
    i2 = cellNodes[i][1];
    POLY_ASSERT(i1 < numNodes and i2 < numNodes);
    ip1 = id2nodes[i1];
    ip2 = id2nodes[i2];
    rp1 = (innerCirc[i1] == 1) ? mCoords.dequantize(&ip1.x) : mOuterCoords.dequantize(&ip1.x);
    rp2 = (innerCirc[i2] == 1) ? mCoords.dequantize(&ip2.x) : mOuterCoords.dequantize(&ip2.x);
    if (orient2d(&rp1.x, &rp2.x, (double*)&points[2*i]) < 0) {
      reverse(cellNodes[i].begin(), cellNodes[i].end());
    }

    // Add first element to end of cell-node list to form BG rings
    cellNodes[i].push_back(cellNodes[i][0]);

    // Blago!
    bool Blago = false;
    if(Blago) cerr << "---------- Cell " << i << " -----------" << endl;
    // Blago!
    
    // Create a BG ring and fill in any ring concept stuff that's missing
    vector<RealPoint> realCellNodes;
    for (vector<unsigned>::const_iterator itr = cellNodes[i].begin();
	 itr != cellNodes[i].end(); 
	 ++itr) {
      ip1 = id2nodes[*itr];
      rp1 = (innerCirc[*itr] == 1) ? 
	mCoords.dequantize(&ip1.x) : 
	mOuterCoords.dequantize(&ip1.x);
      realCellNodes.push_back(rp1);
    }
    POLY_ASSERT(realCellNodes.size() == cellNodes[i].size());
    RealRing realCellRing(realCellNodes.begin(), realCellNodes.end());
    boost::geometry::correct(realCellRing);
    POLY_ASSERT(realCellRing.front() == realCellRing.back());
    //POLY_ASSERT(!boost::geometry::intersects(realCellRing));
    
    // Blago!
    if(Blago) {
    cerr << "\n\nPre-Clip" << endl;
    for (typename RealRing::const_iterator itr = realCellRing.begin();
	 itr != realCellRing.end();
	 ++itr ) cerr << *itr << endl;
    }
    // Blago!

    // Walk node-node pairs and add them according to 4 possible cases
    int numIntersections = 0;
    vector<int> intersectFacets, indices;
    vector<IntPoint> cellBoundary;
    POLY_ASSERT(cellNodes[i].size() > 2);
    for (j = 0; j != cellNodes[i].size()-1; ++j) {
      i1 = cellNodes[i][j  ];
      i2 = cellNodes[i][j+1];
      POLY_ASSERT(i1 != i2);
      POLY_ASSERT(i1 < id2nodes.size() and i2 < id2nodes.size());
      ip1 = id2nodes[i1];
      ip2 = id2nodes[i2];

      // Case 1: Both circumcenters inside bounding box. Add the 2nd point
      if (innerCirc[i1] == 1 and innerCirc[i2] == 1) {
	cellBoundary.push_back(ip2);

	// Blago!
        if(Blago){
	cerr << "Case 1: " 
	     << mCoords.dequantize(&ip1.x) << "  and  " << mCoords.dequantize(&ip2.x) << endl;
        }
	// Blago!

      }
      
      // Case 2: 1st inside, 2nd outside. Find the intersection pt of the 
      // bounding box and the line segment, quantize the pt, and add it.
      // NOTE: Keep track of which facet we exited the bounding box through
      //       and where in the node list we did so.
      else if(innerCirc[i1] == 1 and innerCirc[i2] == 0) {
        ++numIntersections;
        rp1 = mCoords.dequantize(&ip1.x);
        rp2 = mOuterCoords.dequantize(&ip2.x);
	vector<RealType> result;
	vector<int> resultFacets;
	nints = intersectBoundingBox(&rp1.x, &rp2.x, 4, &mCoords.points[0],
				     mCoords.facets, mCoords.delta, 
                                     resultFacets, result);

	// Blago!
        if(Blago){
	cerr << "Case 2: " << rp1.x << " " << rp1.y << "  and  "
	     << rp2.x << " " << rp2.y << endl;
	cerr << "  " << result[0] << "  " << result[1] << endl;
        if (result.size() > 2) cerr << "  " << result[2] << "  " << result[3] << endl;
        }
	// Blago!

        POLY_ASSERT(rp1.x >= mCoords.low[0] and rp1.x <= mCoords.high[0] and
                    rp1.y >= mCoords.low[1] and rp1.y <= mCoords.high[1]);
        POLY_ASSERT(nints == 1 and result.size() == 2 and resultFacets.size() == 1);
	POLY_ASSERT(mCoords.low[0] <= result[0] and result[0] <= mCoords.high[0] and
		    mCoords.low[1] <= result[1] and result[1] <= mCoords.high[1] );
        cellBoundary.push_back(mCoords.quantize(&result[0]));
	intersectFacets.push_back(resultFacets[0]);
        indices.push_back(cellBoundary.size());
      }
      
      // Case 3: 1st outside, 2nd inside. Find the intersection pt of the 
      // bounding box and the line segment, quantize the pt, and add it. Also 
      // add the 2nd point inside the box
      // NOTE: Keep track of which facet we passed through to enter the box
      else if(innerCirc[i1] == 0 and innerCirc[i2] == 1) {
        ++numIntersections;
        rp1 = mOuterCoords.dequantize(&ip1.x);
        rp2 = mCoords.dequantize(&ip2.x);
	vector<RealType> result;
	vector<int> resultFacets;
	nints = intersectBoundingBox(&rp1.x, &rp2.x, 4, &mCoords.points[0],
				     mCoords.facets, mCoords.delta,
                                     resultFacets, result);

	// Blago!
        if(Blago){
	cerr << "Case 3: " << rp1.x << " " << rp1.y << "  and  " 
             << rp2.x << " " << rp2.y << endl;
	cerr << "  " << result[0] << "  " << result[1] << endl;
	}
        // Blago!

        
	POLY_ASSERT2(rp2.x >= mCoords.low[0] and rp2.x <= mCoords.high[0] and
                     rp2.y >= mCoords.low[1] and rp2.y <= mCoords.high[1],
                     "Point " << rp2 << ip2 << " is outside" << endl << mCoords);
	POLY_ASSERT(nints == 1 and result.size() == 2 and resultFacets.size() == 1);
	POLY_ASSERT2(mCoords.low[0] <= result[0] and result[0] <= mCoords.high[0] and
                     mCoords.low[1] <= result[1] and result[1] <= mCoords.high[1],
                     "Intersection point (" << result[0] << "," << result[1] << ") is outside"
                     << endl << mCoords);
	intersectFacets.push_back(resultFacets[0]);
        indices.push_back(-1);
        cellBoundary.push_back(mCoords.quantize(&result[0]));
	cellBoundary.push_back(ip2);
      }
	
      // Case 4: Both outside. Check intersections of the bounding box and the
      // line segment. Quantize and add all intersections
      else {
        rp1 = mOuterCoords.dequantize(&ip1.x);
        rp2 = mOuterCoords.dequantize(&ip2.x);
	vector<RealType> result;
	vector<int> resultFacets;
	nints = intersectBoundingBox(&rp1.x, &rp2.x, 4, &mCoords.points[0],
				     mCoords.facets, mCoords.delta,
                                     resultFacets, result);
	//POLY_ASSERT(nints==0 or nints==2);

	// Blago!
        if(Blago){
	cerr << "Case 4: " << rp1.x << " " << rp1.y << "  and  "
	     << rp2.x << " " << rp2.y << endl;
	}
        // Blago!

	if (nints == 2) {
          numIntersections += nints;
          RealType d1 = geometry::distance<2,RealType>(&result[0],&rp1.x);
	  RealType d2 = geometry::distance<2,RealType>(&result[2],&rp1.x);
	  int enterIndex = (d1 < d2) ? 0 : 1;
	  int exitIndex  = (d1 < d2) ? 1 : 0;
	  POLY_ASSERT(result.size()==4);
	  POLY_ASSERT(mCoords.low[0] <= result[0] and result[0] <= mCoords.high[0] and
		      mCoords.low[1] <= result[1] and result[1] <= mCoords.high[1] and
		      mCoords.low[0] <= result[2] and result[2] <= mCoords.high[0] and
		      mCoords.low[1] <= result[3] and result[3] <= mCoords.high[1]);
          cellBoundary.push_back(mCoords.quantize(&result[2*enterIndex]));
          cellBoundary.push_back(mCoords.quantize(&result[2* exitIndex]));
          intersectFacets.push_back(resultFacets[enterIndex]);
          intersectFacets.push_back(resultFacets[ exitIndex]);
          indices.push_back(-1);
          indices.push_back(cellBoundary.size());
	}
      }
    }

    // // Blago!
    // cerr << "Before adding corners:" << endl;
    // for (int jj = 0; jj != cellBoundary.size(); ++jj){
    //    cerr << cellBoundary[jj].realx(mLow[0],mDelta) << " " 
    //         << cellBoundary[jj].realy(mLow[1],mDelta) << endl;
    // }
    // // Blago!

    // If we exited and re-entered the bounding box while marching through the
    // nodes, we must add all corners of the bounding box between the exit facet
    // and the enter facet, walking CCW. Insert them into the node list.
    if (numIntersections > 0) {
      POLY_ASSERT(numIntersections % 2   == 0               );
      POLY_ASSERT(intersectFacets.size() == numIntersections);
      POLY_ASSERT(indices.size()         == numIntersections);
      int index, start, addCount = 0;
      if (indices[0] < 0) {
         intersectFacets.push_back(intersectFacets[0]);
         indices.push_back(indices[0]);
         start = 1;
      } else {
         start = 0;
      }
      for (j = 0; j != numIntersections/2; ++j) {
	vector<IntPoint> extraBoundaryPoints;
	int exitIndex = 2*j + start;
	int enterIndex = 2*j + start + 1;
	for (k = intersectFacets[exitIndex]; k%4 != intersectFacets[enterIndex]; ++k) {
          index = (k+1)%4;
          extraBoundaryPoints.push_back(mCoords.quantize(&mCoords.points[2*index]));
	}
        POLY_ASSERT(indices[exitIndex] >= 0);
	POLY_ASSERT(indices[exitIndex] + addCount <= cellBoundary.size());
	cellBoundary.insert(cellBoundary.begin() + indices[exitIndex] + addCount,
			    extraBoundaryPoints.begin(),
			    extraBoundaryPoints.end());
	addCount += extraBoundaryPoints.size();
      }
    }
    
    POLY_ASSERT(!cellBoundary.empty());
    cellBoundary.push_back(cellBoundary[0]);
    boost::geometry::assign(cellRings[i], BGring(cellBoundary.begin(),
                                                 cellBoundary.end()));
    boost::geometry::correct(cellRings[i]);
    POLY_ASSERT(!cellRings[i].empty());
    POLY_ASSERT(cellRings[i].front() == cellRings[i].back());
    
    if(Blago){
    // Blago!
    cerr << "Pre-clipped ring:" << endl;
    for (typename BGring::iterator itr = cellRings[i].begin();
         itr != cellRings[i].end(); ++itr) {
      cerr << mCoords.dequantize(&(*itr).x) << endl;
    }
    }
    // Blago!
    */
