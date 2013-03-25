//------------------------------------------------------------------------
// TriangleTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
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

//------------------------------------------------------------------------
// Union a Boost.Geometry ring with a Boost.Geometry multi_polygon.
// The resulting multi_polygon is corrected to ensure it conforms to
// the proper geometric concept.
//------------------------------------------------------------------------
void
createBGUnion( boost::geometry::model::ring<Point2<CoordHash>,false> ring, 
               boost::geometry::model::multi_polygon<
                 boost::geometry::model::polygon<Point2<CoordHash>,false> >& multiPolygon )
{
  boost::geometry::model::multi_polygon<
    boost::geometry::model::polygon<Point2<CoordHash>,false> > temp;
  boost::geometry::union_(multiPolygon, ring, temp);
  boost::geometry::correct(temp);
  multiPolygon=temp;
}


//------------------------------------------------------------------------------
// Store the PLC boundary as a Boost.geometry polygon
//------------------------------------------------------------------------------
template<typename RealType, typename PointType>
void
buildBoostBoundary(const vector<PointType>& IntPLCPoints,
		   const PLC<2,RealType>& geometry,
		   boost::geometry::model::polygon<PointType,false>& boundary)
{
  typedef boost::geometry::model::polygon<PointType,false> BGpolygon;
  int i, j, k;  
  boost::geometry::append( boundary, IntPLCPoints[geometry.facets[0][0]] );
  for (j = 0; j != geometry.facets.size(); ++j){
    POLY_ASSERT(geometry.facets[j].size() == 2);
    i =  geometry.facets[j][1];
    boost::geometry::append( boundary, IntPLCPoints[i]);
  }
  POLY_ASSERT(boundary.outer().size() == geometry.facets.size() + 1);
  POLY_ASSERT(boundary.outer().front() == boundary.outer().back());
  
  // Add any interior holes.
  const unsigned numHoles = geometry.holes.size();
  if (numHoles > 0) {
    typename BGpolygon::inner_container_type& holes = boundary.inners();
    holes.resize(numHoles);
    for (k = 0; k != numHoles; ++k) {
      boost::geometry::append( holes[k], IntPLCPoints[geometry.holes[k][0][0]] );
      for (j = 0; j != geometry.holes[k].size(); ++j) {
	POLY_ASSERT(geometry.holes[k][j].size() == 2);
	i =  geometry.holes[k][j][1];
	boost::geometry::append( holes[k], IntPLCPoints[i]); 
      }
      POLY_ASSERT(holes[k].size() == geometry.holes[k].size() + 1 );
      POLY_ASSERT(holes[k].front() == holes[k].back());
    }
    POLY_ASSERT(boundary.inners().size() == numHoles );
  }
}


//------------------------------------------------------------------------------
// Compute the intersection of a line segment and a PLC. Returns number of
// intersections, intersected facet(s) along PLC, and intersection point(s)
//------------------------------------------------------------------------------
template<typename RealType>
inline
unsigned intersectBoundingBox(const RealType* point1,
			      const RealType* point2,
			      const unsigned numVertices,
			      const RealType* vertices,
			      const std::vector<std::vector<int> >& facets,
			      std::vector<int>& facetIntersections,
			      std::vector<RealType>& result) {
  unsigned i, j, numIntersections=0;
  RealType intersectionPoint[2];
  bool intersects;
  const unsigned numFacets = facets.size();
  for (unsigned ifacet = 0; ifacet != numFacets; ++ifacet) {
    POLY_ASSERT(facets[ifacet].size() == 2);
    i = facets[ifacet][0];
    j = facets[ifacet][1];
    POLY_ASSERT(i >= 0 and i < numVertices);
    POLY_ASSERT(j >= 0 and j < numVertices);
    intersects = geometry::segmentIntersection2D(point1, point2,
						 &vertices[2*i], &vertices[2*j],
						 intersectionPoint);
    if( intersects ){
      RealType p[2];
      p[(ifacet+1)%2] = vertices[2*ifacet + (ifacet+1)%2];
      p[ ifacet   %2] = intersectionPoint[ifacet%2];
      ++numIntersections;
      result.push_back(p[0]);
      result.push_back(p[1]);
//       result.push_back(intersectionPoint[0]);
//       result.push_back(intersectionPoint[1]);
      facetIntersections.push_back(ifacet);
    }
  }
  return numIntersections;
}
  

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
std::vector<unsigned>
computeSortedFaceNodes(const std::vector<std::pair<int, int> >& edges) {
  typedef std::pair<int, int> EdgeHash;
  std::vector<unsigned> result;

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
    
    // Look for any edges with one one node in the set.  There can be at most
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
  Tessellator<2, RealType>() {
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
  POLY_ASSERT(points.size() > 2 );
  
  // Reset the bounding box and quantized mesh spacing
  mLowOuter.assign(  2,  numeric_limits<RealType>::max() );
  mHighOuter.assign( 2, -numeric_limits<RealType>::max() );
  mdxOuter = 0;

  mLowInner.assign(  2,  numeric_limits<RealType>::max() );
  mHighInner.assign( 2, -numeric_limits<RealType>::max() );
  mdxInner = 0;

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
           Tessellation<2, RealType>& mesh) const 
{
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
           Tessellation<2, RealType>& mesh) const 
{
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(!points.empty());
  
  // Reset the bounding box and quantized mesh spacing
  mLowOuter.assign(  2,  numeric_limits<RealType>::max() );
  mHighOuter.assign( 2, -numeric_limits<RealType>::max() );
  mdxOuter = 0;
  
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
		 map<IntPoint, pair<int,int> >& circumcenterMap,
		 map<int, vector<unsigned> >& cellNodes,
		 vector<unsigned>& infNodes) const {

  POLY_ASSERT(!points.empty());
  POLY_ASSERT(points.size() != 2);

  const CoordHash coordMax = (1LL << 26); // numeric_limits<CoordHash>::max() >> 32U;
  const double degeneracy = 1.5e-8;
  
  const unsigned numGenerators = points.size()/2;
  
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
  map<EdgeHash, vector<unsigned> > edge2tri;
  map<int, set<unsigned> > gen2tri;
  int pindex, qindex, rindex, i, j;
  EdgeHash pq, pr, qr;
  // RealType lowc [2] = { numeric_limits<RealType>::max(), 
  //       		numeric_limits<RealType>::max()};
  // RealType highc[2] = {-numeric_limits<RealType>::max(), 
  //       	       -numeric_limits<RealType>::max()};
  RealType lowc [2] = {mLowInner [0], mLowInner [1]};
  RealType highc[2] = {mHighInner[0], mHighInner[1]};
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
    // if (std::abs(orient2d(&delaunay.pointlist[2*pindex],
    //                       &delaunay.pointlist[2*qindex],
    //                       &delaunay.pointlist[2*rindex])) > degeneracy) {
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
    // }
  }
  POLY_ASSERT(circumcenters.size() == delaunay.numberoftriangles);

  // The circumcenters may all lie inside the convex hull of the
  // generators for an unbounded tessellation. Include the generator
  // locations in the high/low search
  geometry::expandBoundingBox<2,RealType>(&delaunay.pointlist[0],
					  2*delaunay.numberofpoints,
					  true, lowc, highc);
  POLY_ASSERT(lowc[0] <= highc[0] and lowc[1] <= highc[1]);
  
  // The bounding box which contains all circumcenters and generators
  RealType cbox[2] = {highc[0] - lowc[0], highc[1] - lowc[1]};
  
  // The bounding circle onto which we project the "infinite" rays of the 
  // unbounded faces of the tessellation.
  const RealType rtmp    = 4.0*max(cbox[0], cbox[1]);
  const RealType ctmp[2] = {0.5*(lowc[0]+highc[0]), 0.5*(lowc[1]+highc[1])};
  
  // We resize mLow and boxsize so that the bounding box
  // contains the "infinite" sphere. mHigh is not really needed.
  mLowOuter [0] = min(mLowOuter [0], ctmp[0]-rtmp);
  mLowOuter [1] = min(mLowOuter [1], ctmp[1]-rtmp);
  mHighOuter[0] = max(mHighOuter[0], ctmp[0]+rtmp);
  mHighOuter[1] = max(mHighOuter[1], ctmp[1]+rtmp);
  
  const RealType rinf     =  0.5*(mHighOuter[0] - mLowOuter[0]);
  const RealType cboxc[2] = {0.5*(mHighOuter[0] + mLowOuter[0]), 
			     0.5*(mHighOuter[1] + mLowOuter[1])};

  const double cboxsize = 2.0*rinf;
  mdxOuter = max(degeneracy, cboxsize/coordMax);

  // // Blago!
  // cerr << "Outer Bounding Box = "
  //      << "(" << mLowOuter[0] << "," << mHighOuter[0] << ")X"
  //      << "(" << mLowOuter[1] << "," << mHighOuter[1] << ")" << endl;
  // cerr << "Outer Mesh Spacing = " << mdxOuter << endl;
  // cerr << "Inner Bounding Box = "
  //      << "(" << mLowInner[0] << "," << mHighInner[0] << ")X"
  //      << "(" << mLowInner[1] << "," << mHighInner[1] << ")" << endl;
  // cerr << "Inner Mesh Spacing = " << mdxInner << endl << endl;
  // // Blago!
  
  // Determine which circumcenters lie inside fine-scale bounding box
  // Map circumcenters and triangle indices to global id's
  int inside;
  IntPoint ip;
  map<IntPoint, int> circ2id;
  map<int, unsigned> tri2id;
  for (i = 0; i != delaunay.numberoftriangles; ++i){
    if (circumcenters[i].x >= mLowInner [0] and
	circumcenters[i].x <= mHighInner[0] and
	circumcenters[i].y >= mLowInner [1] and
	circumcenters[i].y <= mHighInner[1]) {
      inside = 1;
      ip = IntPoint(circumcenters[i].x, circumcenters[i].y,
		    mLowInner[0], mLowInner[1], mdxInner);
    } else {
      inside = 0;
      ip = IntPoint(circumcenters[i].x, circumcenters[i].y,
		    mLowOuter[0], mLowOuter[1], mdxOuter);
    }
    j = internal::addKeyToMap(ip, circ2id);
    tri2id[i] = j;
    if (j == circ2id.size()-1) circumcenterMap[ip] = make_pair(j,inside);
  }
  POLY_ASSERT(circ2id.size() == circumcenterMap.size());

  // The exterior edges of the triangularization have "unbounded" rays, originating
  // at the circumcenter of the corresponding triangle and passing perpendicular to
  // the edge
  bool test;
  RealPoint ehat, pinf;
  map<EdgeHash, unsigned> edge2id;
  int i1, i2, ivert, k;
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
      
      // Get the intersection point along the "infinite" circumcircle
      test = geometry::rayCircleIntersection(&circumcenters[i].x,
                                             &ehat.x,
                                             cboxc,
                                             rinf,
                                             1.0e-6,
                                             &pinf.x);
      POLY_ASSERT(test);
      IntPoint ip(pinf.x, pinf.y, 
		  mLowOuter[0], mLowOuter[1], mdxOuter);
      k = circ2id.size();
      j = internal::addKeyToMap(ip, circ2id);
      if (j == circ2id.size()-1)  circumcenterMap[ip] = make_pair(j,0);
      POLY_ASSERT(edge2id.find(edge) == edge2id.end());
      edge2id[edge] = j;
      if (k != circ2id.size()) infNodes.push_back(1);
    }
  }
  POLY_ASSERT(circ2id.size() == circumcenterMap.size());
  
  // The faces corresponding to each triangle edge
  unsigned ii, jj;
  for (map<int, set<unsigned> >::const_iterator genItr = gen2tri.begin();
       genItr != gen2tri.end(); ++genItr) {
    pindex = genItr->first;
    const set<unsigned>& tris = genItr->second;
    POLY_ASSERT(pindex < numGenerators);
    
    set<EdgeHash> meshEdges;
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
        POLY_ASSERT(jj != ii);
        meshEdges.insert(internal::hashEdge(ii,jj));
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
        POLY_ASSERT(jj != ii);
        meshEdges.insert(internal::hashEdge(ii,jj));
      }else{
         POLY_ASSERT((edge2tri[pr].size() == 2 and edge2tri[pr][0] == i)
                     or edge2tri[pr][1] == i);
        k = (edge2tri[pr][0] == i ? edge2tri[pr][1] : edge2tri[pr][0]);
        jj = tri2id[k];
        if (jj != ii) meshEdges.insert(internal::hashEdge(ii,jj));
      }
    }

    cellNodes[pindex] = 
      computeSortedFaceNodes(vector<EdgeHash>(meshEdges.begin(), meshEdges.end()));
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);

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
computeCellNodesCollinear(const vector<RealType>& points,
			  map<IntPoint, pair<int,int> >& circumcenterMap,
			  map<int, vector<unsigned> >& cellNodes,
			  vector<unsigned>& infNodes) const {  

  const unsigned numGenerators = points.size()/2;
  const CoordHash coordMax = (1LL << 26); // numeric_limits<CoordHash>::max() >> 32U;
  const double degeneracy = 1.5e-8;
  int i;
  
  vector<pair<RealPoint,int> > pointIndexPairs;
  for (i = 0; i != numGenerators; ++i){
    pointIndexPairs.push_back(make_pair(RealPoint(points[2*i], points[2*i+1]), i));
  }
  sort( pointIndexPairs.begin(), pointIndexPairs.end(),
	internal::pairCompareFirst<RealPoint,int> );
  
  // Bounding box for the points
  RealType lowc[2], highc[2];
  geometry::computeBoundingBox<2,RealType>(points, true, lowc, highc);
  POLY_ASSERT(lowc[0] <= highc[0] and lowc[1] <= highc[1]);

  
  // The bounding circle onto which we project the "infinite" rays of the 
  // unbounded faces of the tessellation.
  const RealType cbox[2] = {highc[0] - lowc[0], highc[1] - lowc[1]};
  const RealType rtmp    = 4.0*max(cbox[0], cbox[1]);
  const RealType ctmp[2] = {0.5*(lowc[0]+highc[0]), 
			    0.5*(lowc[1]+highc[1])};

  // We resize mLow and boxsize so that the bounding box
  // contains the "infinite" sphere.
  mLowOuter [0] = min(mLowOuter [0], ctmp[0]-rtmp);
  mLowOuter [1] = min(mLowOuter [1], ctmp[1]-rtmp);
  mHighOuter[0] = max(mHighOuter[0], ctmp[0]+rtmp);  
  mHighOuter[1] = max(mHighOuter[1], ctmp[1]+rtmp);
  
  const RealType rinf     =  0.5*(mHighOuter[0]-mLowOuter[0]);
  const RealType cboxc[2] = {0.5*(mLowOuter[0]+mHighOuter[0]), 
			     0.5*(mLowOuter[1]+mHighOuter[1])};

  const double cboxsize = 2.0*rinf;
  mdxOuter = max(degeneracy, cboxsize/coordMax);

  // // Blago!
  // cerr << "Outer Bounding Box = "
  //      << "(" << mLowOuter[0] << "," << mHighOuter[0] << ")X"
  //      << "(" << mLowOuter[1] << "," << mHighOuter[1] << ")" << endl;
  // cerr << "Outer Mesh Spacing = " << mdxOuter << endl;
  // cerr << "Inner Bounding Box = "
  //      << "(" << mLowInner[0] << "," << mHighInner[0] << ")X"
  //      << "(" << mLowInner[1] << "," << mHighInner[1] << ")" << endl;
  // cerr << "Inner Mesh Spacing = " << mdxInner << endl << endl;
  // // Blago!

  // Number of nodes
  const int nnodes = 2*numGenerators;

  int inside;
  bool test;
  unsigned inode, icell1, icell2;
  RealPoint p1, p2, r1, r2, node, midpt;
  // IntPoint IntNode;
  vector<RealPoint> nodeList(nnodes);
  vector<int> insideList(nnodes);
  infNodes.resize(nnodes);
  
  // ---------------- Nodes and faces for cell 0 ----------------- //

  inode  = 0;
  icell1 = pointIndexPairs[0].second;
  icell2 = pointIndexPairs[1].second;

  // Node position
  p1   = pointIndexPairs[0].first;
  p2   = pointIndexPairs[1].first;
  midpt = RealPoint( 0.5*(p1.x + p2.x),
		     0.5*(p1.y + p2.y) );
  r1.x = p2.x - p1.x;
  r1.y = p2.y - p1.y;
  geometry::unitVector<2,RealType>(&r1.x);
  r2.x =  r1.y;
  r2.y = -r1.x;
  
  // Extra inf node used to bound the first cell
  r1 *= -1.0;
  test = geometry::rayCircleIntersection(&p1.x, &r1.x, cboxc, rinf, 1.0e-10, &node.x);
  POLY_ASSERT(test);
  inside = (node.x >= mLowInner[0] and node.x <= mHighInner[0] and
            node.y >= mLowInner[1] and node.y <= mHighInner[1]) ? 1 : 0;
  nodeList[inode] = node;
  insideList[inode] = inside;

  // IntNode = IntPoint(node.x, node.y, mLowOuter[0], mLowOuter[1], mdxOuter);
  // circumcenterMap[IntNode] = make_pair(inode,0);
  // infNodes.push_back(1);

  // Node 1: endpt of first interior face
  test = geometry::rayCircleIntersection(&midpt.x, &r2.x, cboxc, rinf, 1.0e-10, &node.x);
  POLY_ASSERT(test);
  inside = (node.x >= mLowInner[0] and node.x <= mHighInner[0] and
            node.y >= mLowInner[1] and node.y <= mHighInner[1]) ? 1 : 0;
  nodeList[inode+1] = node;
  insideList[inode+1] = inside;

  // IntNode = IntPoint(node.x, node.y, mLowOuter[0], mLowOuter[1], mdxOuter);
  // circumcenterMap[IntNode] = make_pair(inode+1,0);
  // infNodes.push_back(1);
  
  // Node 2: other endpt of first interior face
  r2 *= -1.0;
  test = geometry::rayCircleIntersection(&midpt.x, &r2.x, cboxc, rinf, 1.0e-10, &node.x);
  POLY_ASSERT(test);
  inside = (node.x >= mLowInner[0] and node.x <= mHighInner[0] and
            node.y >= mLowInner[1] and node.y <= mHighInner[1]) ? 1 : 0;
  nodeList[inode+2] = node;
  insideList[inode+2] = inside;

  // IntNode = IntPoint(node.x, node.y, mLowOuter[0], mLowOuter[1], mdxOuter);
  // circumcenterMap[IntNode] = make_pair(inode+2,0);
  // infNodes.push_back(1);

  // Nodes around cell 0
  cellNodes[icell1].push_back(inode  );
  cellNodes[icell1].push_back(inode+1);
  cellNodes[icell1].push_back(inode+2);

  // Half of the nodes around cell 1
  cellNodes[icell2].push_back(inode+2);
  cellNodes[icell2].push_back(inode+1);
    
  // ------------------ Interior cells ----------------- //

  for (i = 1; i != numGenerators-1; ++i){
    inode  = 2*i+1;
    icell1 = pointIndexPairs[i  ].second;
    icell2 = pointIndexPairs[i+1].second;
    
    p1    = pointIndexPairs[i  ].first;
    p2    = pointIndexPairs[i+1].first;
    midpt = RealPoint( 0.5*(p1.x + p2.x),
                       0.5*(p1.y + p2.y) );
    r1.x = p2.x - p1.x;
    r1.y = p2.y - p1.y;
    geometry::unitVector<2,RealType>(&r1.x);
    r2.x =  r1.y;
    r2.y = -r1.x;
    
    // Node 0: endpt of interior face
    test = geometry::rayCircleIntersection(&midpt.x, &r2.x, cboxc, rinf, 1.0e-10, &node.x);
    POLY_ASSERT(test);
    inside = (node.x >= mLowInner[0] and node.x <= mHighInner[0] and
              node.y >= mLowInner[1] and node.y <= mHighInner[1]) ? 1 : 0;
    nodeList[inode] = node;
    insideList[inode] = inside;

    // IntNode = IntPoint(node.x, node.y, mLowOuter[0], mLowOuter[1], mdxOuter);
    // circumcenterMap[IntNode] = make_pair(inode,0);
    // infNodes.push_back(1);
    
    // Node 1: other endpt of interior face
    r2 *= -1.0;
    test = geometry::rayCircleIntersection(&midpt.x, &r2.x, cboxc, rinf, 1.0e-10, &node.x);
    POLY_ASSERT(test);
    inside = (node.x >= mLowInner[0] and node.x <= mHighInner[0] and
              node.y >= mLowInner[1] and node.y <= mHighInner[1]) ? 1 : 0;
    nodeList[inode+1] = node;
    insideList[inode+1] = inside;
    
    // IntNode = IntPoint(node.x, node.y, mLowOuter[0], mLowOuter[1], mdxOuter);
    // circumcenterMap[IntNode] = make_pair(inode+1,0);
    // infNodes.push_back(1);

    // Other half of the nodes around cell i
    cellNodes[icell1].push_back(inode  );
    cellNodes[icell1].push_back(inode+1);

    // Half of the nodes around cell i+1
    cellNodes[icell2].push_back(inode+1);
    cellNodes[icell2].push_back(inode  );
  }
 
  // ------------- Nodes and faces for final cell ----------------- //
  
  inode  = 2*numGenerators-1;
  icell1 = pointIndexPairs[numGenerators-1].second;
  
  // Node position
  p1   = pointIndexPairs[numGenerators-1].first;
  p2   = pointIndexPairs[numGenerators-2].first;
  r1.x = p1.x - p2.x;
  r1.y = p1.y - p2.y;
  geometry::unitVector<2,RealType>(&r1.x);
  
  test = geometry::rayCircleIntersection(&p2.x, &r1.x, cboxc, rinf, 1.0e-10, &node.x);
  POLY_ASSERT(test);
  inside = (node.x >= mLowInner[0] and node.x <= mHighInner[0] and
            node.y >= mLowInner[1] and node.y <= mHighInner[1]) ? 1 : 0;
  nodeList[inode] = node;
  insideList[inode] = inside;
    
  // IntNode = IntPoint(node.x, node.y, mLowOuter[0], mLowOuter[1], mdxOuter);
  // circumcenterMap[IntNode] = make_pair(inode,0);
  // infNodes.push_back(1);
  
  // Last node for final cell
  cellNodes[icell1].push_back(inode);

  POLY_ASSERT(cellNodes.size() == numGenerators);
  POLY_ASSERT(nodeList.size() == nnodes);
  POLY_ASSERT(insideList.size() == nnodes);

  IntPoint IntNode;
  for (i = 0; i != nnodes; ++i) {
    IntNode = (insideList[i] == 1) ?
       IntPoint(nodeList[i].x, nodeList[i].y, mLowInner[0], mLowInner[1], mdxInner) :
       IntPoint(nodeList[i].x, nodeList[i].y, mLowOuter[0], mLowOuter[1], mdxOuter);
    circumcenterMap[IntNode] = make_pair(i,insideList[i]);
    infNodes[i] = (insideList[i] == 1) ? 0 : 1;
  }
  POLY_ASSERT(circumcenterMap.size() == nnodes);

//   // Blago!
//   for (i = 0; i < cellNodes.size(); ++i){
//     cerr << "Cell " << i << endl;
//     for (int j = 0; j < cellNodes[i].size(); ++j){
//       cerr << "  " << cellNodes[i][j];
//     }
//     cerr << endl;
//   }
//   // Blago!
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeCellRings(const vector<RealType>& points,
		 const vector<RealType>& PLCpoints,
		 const PLC<2, RealType>& geometry,
		 vector<BGring>& cellRings,
		 bool performCellAdoption) const {
  POLY_ASSERT(!points.empty());
  POLY_ASSERT2(!PLCpoints.empty(), "Error: attempting to create a bounded "
               << "tessellation with no bounding points");
  POLY_ASSERT2(!geometry.empty(),  "Error: attempting to create a bounded "
               << "tessellation with no bounding PLC");

  // Find the range of the generator points.
  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints = PLCpoints.size()/2;
  int i, j, k;
  
  mLowOuter [0] = min(mLowOuter [0], mLowInner [0]);
  mLowOuter [1] = min(mLowOuter [1], mLowInner [1]);
  mHighOuter[0] = max(mHighOuter[0], mHighInner[0]);
  mHighOuter[1] = max(mHighOuter[1], mHighInner[1]);

  // Get the node IDs and quantized positions around each generator
  map<IntPoint, pair<int, int> > circumcenterMap;
  map<int, vector<unsigned> > cellNodes;
  vector<unsigned> infNodes;
  if (numGenerators == 2) {
    this->computeCellNodesCollinear(points, circumcenterMap, cellNodes, infNodes);
  }else{
    // Check that the points are not all collinear. If they are, the mesh is 1D
    // degenerate, and we cannot compute the Delaunay.
    bool collinear = true;
    i = 2;
    while (collinear and i != numGenerators){
      collinear *= geometry::collinear<2,RealType>(&points[0], &points[2], &points[2*i], 1.0e-10);
      ++i;
    }

    if (collinear){
      this->computeCellNodesCollinear(points, circumcenterMap, cellNodes, infNodes);
    }else{
      this->computeCellNodes(points, circumcenterMap, cellNodes, infNodes);
    }
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);
  
  // Quantize the PLCpoints
  vector<IntPoint> IntPLCPoints(numPLCpoints);
  for (i = 0; i < numPLCpoints; ++i){
    IntPLCPoints[i] = IntPoint( PLCpoints[2*i], PLCpoints[2*i+1],
				mLowInner[0], mLowInner[1], mdxInner );
  }

  // Generate the quantized boundary to handle boost intersections
  BGpolygon boundary;
  buildBoostBoundary(IntPLCPoints, geometry, boundary);
  
  // Form a PLC for the inner bounding box. Circumcenters that lie outside
  // the bounding box of the PLC boundary are quantized based on different
  // criteria to avoid contaminating the degeneracy spacing of the mesh nodes.
  // We will project these outer circumcenters to the edges of the bounding
  // box so all nodes follow the input degeneracy spacing.
  vector<RealType> bbPoints;
  bbPoints.push_back(mLowInner [0]);  bbPoints.push_back(mLowInner [1]);
  bbPoints.push_back(mHighInner[0]);  bbPoints.push_back(mLowInner [1]);
  bbPoints.push_back(mHighInner[0]);  bbPoints.push_back(mHighInner[1]);
  bbPoints.push_back(mLowInner [0]);  bbPoints.push_back(mHighInner[1]);
  PLC<2,RealType> boundingBox;
  boundingBox.facets.resize(4, vector<int>(2));
  for (i = 0; i != 4; ++i) {
    boundingBox.facets[i][0] = i;
    boundingBox.facets[i][1] = (i+1) % 4;
  }
  
  // Create a reverse look-up map of IDs to circumcenters
  POLY_ASSERT(circumcenterMap.size() > 0);
  map<int, IntPoint> circumcenters;
  vector<int> innerCirc(circumcenterMap.size());
  for (map<IntPoint, pair<int,int> >::const_iterator itr = circumcenterMap.begin();
       itr != circumcenterMap.end(); ++itr) {
    i = itr->second.first;
    circumcenters[i] = itr->first;
    innerCirc[i] = itr->second.second;
  }
  POLY_ASSERT(circumcenters.size() == circumcenterMap.size());  

  // Walk the nodes around each generator and build up the cell ring
  IntPoint X, ip1, ip2;
  RealPoint rp1, rp2;
  unsigned i1, i2, nints;
  bool inside;
  vector<BGring> orphanage;
  cellRings.resize(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
    // Check the orientation of the node list and reverse it if it's CW
    POLY_ASSERT(cellNodes[i].size() > 1);
    i1 = cellNodes[i][0];
    i2 = cellNodes[i][1];
    ip1 = circumcenters[i1];
    ip2 = circumcenters[i2];
    rp1 = (innerCirc[i1] == 1)
      ? RealPoint(ip1.realx(mLowInner[0], mdxInner),
		  ip1.realy(mLowInner[1], mdxInner))
      : RealPoint(ip1.realx(mLowOuter[0], mdxOuter),
		  ip1.realy(mLowOuter[1], mdxOuter));
    rp2 = (innerCirc[i2] == 1)
      ? RealPoint(ip2.realx(mLowInner[0], mdxInner),
		  ip2.realy(mLowInner[1], mdxInner))
      : RealPoint(ip2.realx(mLowOuter[0], mdxOuter),
		  ip2.realy(mLowOuter[1], mdxOuter));
    if (orient2d(&rp1.x, &rp2.x, (double*)&points[2*i]) < 0) {
      reverse(cellNodes[i].begin(), cellNodes[i].end());
    }

    // Add first element to end of cell-node list to form BG rings
    cellNodes[i].push_back(cellNodes[i][0]);

    // // Blago!
    // cerr << "----- Cell " << i << " ------" << endl;
    // // Blago!

    // Walk node-node pairs and add them according to 4 possible cases
    int numIntersections = 0;
    vector<int> intersectFacets, indices;
    vector<IntPoint> cellBoundary;
    POLY_ASSERT(cellNodes[i].size() > 2);
    for (j = 0; j != cellNodes[i].size()-1; ++j) {
      i1 = cellNodes[i][j  ];
      i2 = cellNodes[i][j+1];
      POLY_ASSERT(i1 != i2);
      POLY_ASSERT(i1 < circumcenters.size() and i2 < circumcenters.size());
      ip1 = circumcenters[i1];
      ip2 = circumcenters[i2];

      // Case 1: Both circumcenters inside bounding box. Add the 2nd point
      if (innerCirc[i1] == 1 and innerCirc[i2] == 1) {
	cellBoundary.push_back(ip2);
	// // Blago!
	// cerr << "Case 1: " 
	//      << ip1.realx(mLowInner  [0],mdxInner) << " "
	//      << ip1.realy(mLowInner  [1],mdxInner) << "  and  "
	//      << ip2.realx(mLowInner  [0],mdxInner) << " "
	//      << ip2.realy(mLowInner  [1],mdxInner) << endl;
	// // Blago!
      }
      
      // Case 2: 1st inside, 2nd outside. Find the intersection pt of the 
      // bounding box and the line segment, quantize the pt, and add it.
      // NOTE: Keep track of which facet we exited the bounding box through
      //       and where in the node list we did so.
      else if(innerCirc[i1] == 1 and innerCirc[i2] == 0) {
        ++numIntersections;
	rp1 = RealPoint(ip1.realx(mLowInner[0],mdxInner), 
			ip1.realy(mLowInner[1],mdxInner));
	rp2 = RealPoint(ip2.realx(mLowOuter[0],mdxOuter), 
			ip2.realy(mLowOuter[1],mdxOuter));
	vector<RealType> result;
	vector<int> resultFacets;
	nints = intersectBoundingBox(&rp1.x, &rp2.x, 4, &bbPoints[0],
				     boundingBox.facets, resultFacets, result);
	// // Blago!
	// cerr << "Case 2: " 
	//      << rp1.x << " "
	//      << rp1.y << "  and  "
	//      << rp2.x << " "
	//      << rp2.y << endl;
	// cerr << "  " << result[0] << "  " << result[1] << endl;
	// // Blago!
	POLY_ASSERT(nints == 1 and result.size() == 2 and resultFacets.size() == 1);
	POLY_ASSERT(mLowInner[0] <= result[0] and result[0] <= mHighInner[0] and
		    mLowInner[1] <= result[1] and result[1] <= mHighInner[1] );
        cellBoundary.push_back(IntPoint(result[0], result[1], 
					mLowInner[0], mLowInner[1], mdxInner));
	intersectFacets.push_back(resultFacets[0]);
        indices.push_back(cellBoundary.size());
      }
      
      // Case 3: 1st outside, 2nd inside. Find the intersection pt of the 
      // bounding box and the line segment, quantize the pt, and add it. Also 
      // add the 2nd point inside the box
      // NOTE: Keep track of which facet we passed through to enter the box
      else if(innerCirc[i1] == 0 and innerCirc[i2] == 1) {
        ++numIntersections;
	rp1 = RealPoint(ip1.realx(mLowOuter[0],mdxOuter), 
			ip1.realy(mLowOuter[1],mdxOuter));
	rp2 = RealPoint(ip2.realx(mLowInner[0],mdxInner  ), 
			ip2.realy(mLowInner[1],mdxInner  ));
	vector<RealType> result;
	vector<int> resultFacets;
	nints = intersectBoundingBox(&rp1.x, &rp2.x, 4, &bbPoints[0],
				     boundingBox.facets, resultFacets, result);
	// // Blago!
	// cerr << "Case 3: " 
	//      << rp1.x << " "
	//      << rp1.y << "  and  "
	//      << rp2.x << " "
	//      << rp2.y << endl;
	// cerr << "  " << result[0] << "  " << result[1] << endl;
	// // Blago!
	POLY_ASSERT(nints == 1 and result.size() == 2 and resultFacets.size() == 1);
	POLY_ASSERT(mLowInner[0] <= result[0] and result[0] <= mHighInner[0] and
		    mLowInner[1] <= result[1] and result[1] <= mHighInner[1] );
	intersectFacets.push_back(resultFacets[0]);
        indices.push_back(-1);
        cellBoundary.push_back(IntPoint(result[0], result[1], 
					mLowInner[0], mLowInner[1], mdxInner));
	cellBoundary.push_back(ip2);
      }
	
      // Case 4: Both outside. Check intersections of the bounding box and the
      // line segment. Quantize and add all intersections
      else {
	rp1 = RealPoint(ip1.realx(mLowOuter[0],mdxOuter), 
			ip1.realy(mLowOuter[1],mdxOuter));
	rp2 = RealPoint(ip2.realx(mLowOuter[0],mdxOuter), 
			ip2.realy(mLowOuter[1],mdxOuter));
	vector<RealType> result;
	vector<int> resultFacets;
	nints = intersectBoundingBox(&rp1.x, &rp2.x, 4, &bbPoints[0],
				     boundingBox.facets, resultFacets, result);
	POLY_ASSERT(nints==0 or nints==2);
	// // Blago!
	// cerr << "Case 4: " 
	//      << rp1.x << " "
	//      << rp1.y << "  and  "
	//      << rp2.x << " "
	//      << rp2.y << endl;
	// // Blago!
	if (nints == 2) {
          numIntersections += nints;
          RealType d1 = geometry::distance<2,RealType>(&result[0],&rp1.x);
	  RealType d2 = geometry::distance<2,RealType>(&result[2],&rp1.x);
	  int enterIndex = (d1 < d2) ? 0 : 1;
	  int exitIndex  = (d1 < d2) ? 1 : 0;
	  POLY_ASSERT(result.size()==4);
	  POLY_ASSERT(mLowInner[0] <= result[0] and result[0] <= mHighInner[0] and
		      mLowInner[1] <= result[1] and result[1] <= mHighInner[1] and
		      mLowInner[0] <= result[2] and result[2] <= mHighInner[0] and
		      mLowInner[1] <= result[3] and result[3] <= mHighInner[1]);
	  cellBoundary.push_back(IntPoint(result[2*enterIndex], result[2*enterIndex+1], 
					  mLowInner[0], mLowInner[1], mdxInner));
	  cellBoundary.push_back(IntPoint(result[2*exitIndex], result[2*exitIndex+1], 
					  mLowInner[0], mLowInner[1], mdxInner));
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
    //    cerr << cellBoundary[jj].realx(mLowInner[0],mdxInner) << " " 
    //         << cellBoundary[jj].realy(mLowInner[1],mdxInner) << endl;
    // }
    // // Blago!


    // If we exited and re-entered the bounding box while marching through the
    // nodes, we must add all corners of the bounding box between the exit facet
    // and the enter facet, walking CCW. Insert them into the node list.
    if (numIntersections > 0) {
      POLY_ASSERT(numIntersections % 2 == 0);
      POLY_ASSERT(intersectFacets.size() == numIntersections);
      POLY_ASSERT(indices.size() == numIntersections);
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
	  extraBoundaryPoints.push_back(IntPoint(bbPoints[2*index], bbPoints[2*index+1],
						 mLowInner[0], mLowInner[1], mdxInner));
	}
        POLY_ASSERT(indices[exitIndex] >= 0);
	POLY_ASSERT(indices[exitIndex] + addCount <= cellBoundary.size());
	cellBoundary.insert(cellBoundary.begin() + indices[exitIndex] + addCount,
			    extraBoundaryPoints.begin(),
			    extraBoundaryPoints.end());
	addCount += extraBoundaryPoints.size();
      }
    }
    

    POLY_ASSERT(cellBoundary.size() > 0);
    cellBoundary.push_back(cellBoundary[0]);
    boost::geometry::assign(cellRings[i], BGring(cellBoundary.begin(), cellBoundary.end()));
    boost::geometry::correct(cellRings[i]);
    POLY_ASSERT(cellRings[i].size() > 0);
    
    // // Blago!
    // cerr << "After adding corners:" << endl;
    // for (typename BGring::iterator itr = cellRings[i].begin();
    //      itr != cellRings[i].end(); ++itr) {
    //   cerr << (*itr).realx(mLowInner[0], mdxInner) << " "
    //        << (*itr).realy(mLowInner[1], mdxInner) << endl;
    // }
    // // Blago!

    // Intersect with the boundary to get the bounded cell.
    // Since for complex boundaries this may return more than one polygon, we find
    // the one that contains the generator.
    vector<BGring> cellIntersections;
    boost::geometry::intersection(boundary, cellRings[i], cellIntersections);
    if (cellIntersections.size() == 0) {
      cerr << points[2*i] << " " << points[2*i+1] << endl 
           << boost::geometry::dsv(cellRings[i]) << endl
           << boost::geometry::dsv(boundary) << endl;
    }
    POLY_ASSERT(cellIntersections.size() > 0);
    if (cellIntersections.size() == 1) {
      cellRings[i] = cellIntersections[0];
    } else {
      X = IntPoint(points[2*i], points[2*i+1], 
		   mLowInner[0], mLowInner[1], mdxInner);
      k = cellIntersections.size();
      for (j = 0; j != cellIntersections.size(); ++j) {
        inside = boost::geometry::within(X, cellIntersections[j]);
        if( inside )  k = j;
	else          orphanage.push_back( cellIntersections[j] );
      }
      // If none of the cell intersections contain the point, check if
      // it's on the boundary of the intersection itself
      if (k == cellIntersections.size()) {
         for (j = 0; j != cellIntersections.size(); ++j) {
            bool onBoundary = false;
            for (typename BGring::const_iterator itr = cellIntersections[j].begin();
                 itr != cellIntersections[j].end(); ++itr) {
               onBoundary += (X == *itr);
            }
            if (onBoundary) k = j;
         }
      }
      POLY_ASSERT(k < cellIntersections.size());
      cellRings[i] = cellIntersections[k];
    }
  }
  
  // If any orphaned cells exist, run the adoption algorithm
  // and modify the neighboring cell rings
  if (performCellAdoption and orphanage.size() > 0){
    this->adoptionAlgorithm(points, cellRings, orphanage);
  }
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

  // Make sure we're not modifying an existing tessellation.
  POLY_ASSERT(mesh.empty());

  int i, j;
  const unsigned numGenerators = points.size()/2;
  map<IntPoint, pair<int, int> > circumcenterMap;
  map<int, vector<unsigned> > cellNodes;
  vector<unsigned> infNodes;
  if (numGenerators == 2) {
    this->computeCellNodesCollinear(points, circumcenterMap, cellNodes, infNodes);
  }else{
    // Check that the points are not all collinear. If they are, the mesh is 1D
    // degenerate, and we cannot compute the Delaunay.
    bool collinear = true;
    i = 2;
    while (collinear and i != numGenerators){
      collinear *= geometry::collinear<2,RealType>(&points[0], &points[2], &points[2*i], 1.0e-10);
      ++i;
    }

    if (collinear){
      this->computeCellNodesCollinear(points, circumcenterMap, cellNodes, infNodes);
    }else{
      this->computeCellNodes(points, circumcenterMap, cellNodes, infNodes);
    }
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);
  mesh.infNodes = infNodes;

  // // Blago!
  // for (i = 0; i != cellNodes.size(); ++i){
  //   cerr << "Cell " << i << endl;
  //   for (j = 0; j != cellNodes[i].size(); ++j){
  //     cerr << "  " << cellNodes[i][j];
  //   }
  //   cerr << endl;
  // }
  // // Blago!


  // Copy the quantized nodes to the final tessellation.
  int inside;
  const unsigned numNodes = circumcenterMap.size();
  mesh.nodes.resize(2*numNodes);
  for (map<IntPoint, pair<int,int> >::const_iterator itr = circumcenterMap.begin();
       itr != circumcenterMap.end(); ++itr) {
    i = itr->second.first;
    inside = itr->second.second;
    POLY_ASSERT(i >= 0 and i < numNodes);
    POLY_ASSERT(inside == 0 or inside == 1);
    if (inside == 1) {
      mesh.nodes[2*i  ] = itr->first.realx(mLowInner[0], mdxInner);
      mesh.nodes[2*i+1] = itr->first.realy(mLowInner[1], mdxInner);
    } else {
      mesh.nodes[2*i  ] = itr->first.realx(mLowOuter[0], mdxOuter);
      mesh.nodes[2*i+1] = itr->first.realy(mLowOuter[1], mdxOuter);
    } 
  }
  
  unsigned i1, i2, iface;
  EdgeHash face;
  vector<unsigned> faceVec(2);
  map<EdgeHash, int> face2id;
  internal::CounterMap<unsigned> faceCounter;
  mesh.cells.resize(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
    // Check the orientation of the nodes around the generator. Reverse the order
    // if they're ordered clockwise. It is enough to check the order of the first
    // two nodes in the cellNodes list, since its node list is sequential
    POLY_ASSERT(cellNodes[i].size() > 1);
    if( orient2d(&mesh.nodes[2*cellNodes[i][0]], 
		 &mesh.nodes[2*cellNodes[i][1]], 
		 (double*)&points[2*i]) < 0 ){
      reverse(cellNodes[i].begin(), cellNodes[i].end());
    }
    
    // The ordered mesh nodes around a given generator
    POLY_ASSERT(cellNodes[i].size() > 2);
    for (j = 0; j != cellNodes[i].size(); ++j){
      i1 = cellNodes[i][j];
      i2 = cellNodes[i][(j+1) % cellNodes[i].size()];
      face  = internal::hashEdge(i1,i2);
      iface = internal::addKeyToMap(face,face2id);
      ++faceCounter[iface];
      
      // If you're looking at this face for the first time, add its
      // nodes, classify it as an interior/inf face, and resize faceCells
      POLY_ASSERT( faceCounter[iface] > 0 );
      if( faceCounter[iface] == 1 ){
        faceVec[0] = face.first; faceVec[1] = face.second;
        mesh.faces.push_back(faceVec);
        mesh.faceCells.resize(iface+1);
	mesh.infFaces.resize(iface+1);
//         if ( mesh.infNodes[i1] == 1 and mesh.infNodes[i2] == 1 ){
//           mesh.infFaces.push_back(1);
//         }else{
//           mesh.infFaces.push_back(0);
//         }
      }
      
      
      // Store the cell-face info based on face orientation
      if( orient2d(&mesh.nodes[2*face.first], 
                   &mesh.nodes[2*face.second], 
                   (double*)&points[2*i]) > 0 ){
        mesh.cells[i].push_back(iface);
        mesh.faceCells[iface].push_back(i);
      }else{
        mesh.cells[i].push_back(~iface);
        mesh.faceCells[iface].push_back(~i);
      }

      // Store the infFace info
      if (faceCounter[iface] == 1 and 
	  mesh.infNodes[i1]  == 1 and 
	  mesh.infNodes[i2]  == 1){
	mesh.infFaces[iface] = 1;
      }else{
	mesh.infFaces[iface] = 0;
      }
    }
  }
  POLY_ASSERT(mesh.faceCells.size() == mesh.faces.size());
  POLY_ASSERT(mesh.infFaces.size() == mesh.faces.size());
  POLY_ASSERT(mesh.infNodes.size() == mesh.nodes.size()/2);
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

  const CoordHash coordMax = (1LL << 26); // numeric_limits<CoordHash>::max() >> 32U;
  const double degeneracy = 1.5e-8;

  // Reset the inner bounding box and mesh spacing
  mLowInner.assign ( 2,  numeric_limits<RealType>::max() );
  mHighInner.assign( 2, -numeric_limits<RealType>::max() );
  mdxInner = 0;
  
  RealType low[2], high[2];
  geometry::computeBoundingBox<2,RealType>(PLCpoints, true, low, high);
  POLY_ASSERT(low[0] < high[0] and low[1] < high[1]);  
  
  const RealType box[2] = {high[0] - low[0],  high[1] - low[1]};
  const RealType cen[2] = {0.5*(high[0] + low[0]),
                           0.5*(high[1] + low[1])};
  const double boxsize  = 2.0*max(box[0], box[1]);

  mLowInner [0] = min(mLowInner [0], cen[0] - boxsize);
  mHighInner[0] = max(mHighInner[0], cen[0] + boxsize);
  mLowInner [1] = min(mLowInner [1], cen[1] - boxsize);  
  mHighInner[1] = max(mHighInner[1], cen[1] + boxsize);

  // const RealType box[2] = {mHighInner[0] - mLowInner[0],
  //       		   mHighInner[1] - mLowInner[1]};
  // const double boxsize  = max(box[0], box[1]);
  mdxInner = max(degeneracy, boxsize/coordMax);

  // Find the range of the generator points.
  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints = PLCpoints.size()/2;
  int i, j, k;

  // Compute bounded cell rings
  vector<BGring> cellRings;
  this->computeCellRings(points, PLCpoints, geometry, cellRings, true);
  
  // Now build the unique mesh nodes and cell info.
  map<IntPoint, int> point2node;
  map<EdgeHash, int> edgeHash2id;
  map<int, vector<int> > edgeCells;
  int iedge;
  mesh.clear();
  POLY_ASSERT(mesh.empty());
  mesh.cells = vector<vector<int> >(numGenerators);
  for (i = 0; i != numGenerators; ++i) { 
    for (typename BGring::const_iterator itr = cellRings[i].begin();
         itr != cellRings[i].end()-1; ++itr) {
      const IntPoint& pX1 = *itr;
      const IntPoint& pX2 = *(itr + 1);
      POLY_ASSERT(*itr != *(itr + 1));
      j = internal::addKeyToMap(pX1, point2node);
      k = internal::addKeyToMap(pX2, point2node);
      POLY_ASSERT(j != k);
      iedge = internal::addKeyToMap(internal::hashEdge(j, k), edgeHash2id);
      edgeCells[iedge].push_back(j < k ? i : ~i);
      mesh.cells[i].push_back(j < k ? iedge : ~iedge);
    }
    POLY_ASSERT(mesh.cells[i].size() >= 3);
  }
  POLY_ASSERT(edgeCells.size() == edgeHash2id.size());
  
  // Fill in the mesh nodes.
  RealType node[2];
  mesh.nodes = vector<RealType>(2*point2node.size());
  for (typename map<IntPoint, int>::const_iterator itr = point2node.begin();
       itr != point2node.end(); ++itr) {
    const IntPoint& p = itr->first;
    i = itr->second;
    POLY_ASSERT(i < mesh.nodes.size()/2);
    node[0] = p.realx(mLowInner[0],mdxInner);
    node[1] = p.realy(mLowInner[1],mdxInner);
    
    // Check if nodes are inside boundary (either bounding box or PLC, if defined)
    //bool inside = boost::geometry::within(RealPoint(node[0],node[1]), realBoundary);
    bool inside = within(node, numPLCpoints, &PLCpoints[0], geometry);
    
    if( !inside ){
      RealType result[2];
      RealType dist = nearestPoint( node, numPLCpoints, &PLCpoints[0], geometry, result );
      // Check the node has not moved more than 2.5 quantized mesh spacings. NOTE: this is
      // not a sharp estimate. Theoreticallly, the distance ought to be at most sqrt(2)*cdx, 
      // but nodes will fail this strict of a test.
      POLY_ASSERT2( dist < 5.0*mdxInner,
                    dist << " " << 2.5*mdxInner << " : "
                    << "(" << node[0]   << " " << node[1]   << ") "
                    << "(" << result[0] << " " << result[1] << ")\n" << geometry );
      node[0] = result[0];
      node[1] = result[1];
    }
     
    mesh.nodes[2*i]   = node[0];
    mesh.nodes[2*i+1] = node[1];

    //cerr << "Node " << i << ": (" << mesh.nodes[2*i] << "," << mesh.nodes[2*i+1] << ")" << endl;
  }
  
  // Fill in the mesh faces.
  mesh.faces = vector<vector<unsigned> >(edgeHash2id.size());
  for (typename map<EdgeHash, int>::const_iterator itr = edgeHash2id.begin();
       itr != edgeHash2id.end(); ++itr) {
    const EdgeHash& ehash = itr->first;
    i = itr->second;
    POLY_ASSERT(i < mesh.faces.size());
    POLY_ASSERT(mesh.faces[i].size() == 0);
    mesh.faces[i].push_back(ehash.first);
    mesh.faces[i].push_back(ehash.second);
  }

  // Fill in the mesh faceCells.
  mesh.faceCells = vector<vector<int> >(edgeHash2id.size());
  for (i = 0; i != mesh.faces.size(); ++i) {
    if (not(edgeCells[i].size() == 1 or edgeCells[i].size() == 2)) {
      const int n1 = mesh.faces[i][0], n2 = mesh.faces[i][1];
      cerr << "Blago! " << i << " " << edgeCells[i].size() << " : " << n1 << " " << n2 << " : ("
           << mesh.nodes[2*n1] << " " << mesh.nodes[2*n1 + 1] << ") ("
           << mesh.nodes[2*n2] << " " << mesh.nodes[2*n2 + 1] << ")" << endl;
      for (j = 0; j != edgeCells[i].size(); ++j) cerr << " --> " << edgeCells[i][j] << " " << points[2*edgeCells[i][j]] << " " << points[2*edgeCells[i][j]+1] << endl;
    }
    POLY_ASSERT(edgeCells[i].size() == 1 or edgeCells[i].size() == 2);
    mesh.faceCells[i] = edgeCells[i];
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeDelaunay(const vector<RealType>& points,
                triangulateio& delaunay) const 
{
  triangulateio in;
   
  // Find the range of the generator points.
  const unsigned numGenerators = points.size()/2;

  // Add the generators
  in.numberofpoints = numGenerators;
  in.pointlist = new RealType[2*in.numberofpoints];
  copy(points.begin(), points.end(), in.pointlist);
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
  triangulate((char*)"Qzec", &in, &delaunay, 0);
  
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
// The algorithm for assigning regions of the domain orphaned by our boundary
// clipping strategy to neighboring cells. Details:
//
// Make sub-tessellations using the generators neighboring the orphaned pieces. 
// Construct a PLC boundary for the sub-tessellation by using the geometry 
// obtained by union-ing the orphan with its neighboring cells. The way we
// compute cell neighbors should ensure that the union gives a contiguous 
// geometry with no holes
//
// TODO: Lift this out of here and make standalone class
//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
adoptionAlgorithm(const vector<RealType>& points,
		  vector<BGring>& cellRings,
		  vector<BGring>& orphanage) const {
  POLY_ASSERT(orphanage.size() > 0);
  int i;
  const int numGenerators = points.size()/2;
  
  // Construct map from node points to neighboring cells
  map<IntPoint, set<int> > point2neighbors;
  for (i = 0; i != cellRings.size(); ++i){
    for (typename BGring::const_iterator itr = cellRings[i].begin();
	 itr != cellRings[i].end() - 1; ++itr){
      point2neighbors[*itr].insert(i);
    }
  }
  
  // First aglomerate all orphaned pieces that neighbor one another by searching orphan pairs
  // TODO: Figure out a more efficient way to do this operation
  BGmulti_polygon orphanUnion;
  for (i = 0; i != orphanage.size(); ++i) {
    createBGUnion(orphanage[i], orphanUnion);
  }
  
  
  for (i = 0; i < orphanUnion.size(); ++i) {
    BGring orphan = orphanUnion[i].outer();
    set<int> orphanNeighbors;
    for (typename BGring::const_iterator pointItr = orphan.begin();
	 pointItr != orphan.end()-1; ++pointItr) {
      map<IntPoint, set<int> >::iterator it = point2neighbors.find(*pointItr);
      if (it != point2neighbors.end()){
	set<int> neighborSet = it->second;
	orphanNeighbors.insert(neighborSet.begin(), neighborSet.end());
      }
    }
    POLY_ASSERT( orphanNeighbors.size() > 0 );
            
    // If the orphan only has a single neighbor, we can skip a lot of work.
    // No need to tessellate - simply union the orphan with its neighbor cell.
    vector<BGring> subCellRings;
    if (orphanNeighbors.size() > 1){
      
      // Compute the sub-tessellation from orphan's neighboring points. Union the
      // orphan and its immediate neighbors to get the tessellation boundary
      vector<RealType> subpoints;
      BGmulti_polygon neighborCells;
      createBGUnion(orphan,neighborCells);          
      for (set<int>::const_iterator nbItr = orphanNeighbors.begin();
	   nbItr != orphanNeighbors.end(); ++nbItr){
	subpoints.push_back( points[2*(*nbItr)  ] );
	subpoints.push_back( points[2*(*nbItr)+1] );
	createBGUnion(cellRings[*nbItr],neighborCells);
      }
      POLY_ASSERT2( neighborCells.size() > 0, "Union produced empty set!" );
      if ( neighborCells.size() > 1){
	cerr << "Blago!" << endl;
	for (i = 0; i != neighborCells.size(); ++i){
	  cerr << "Sub-polygon " << i << " in the union has bounding ring" << endl;
	  for (typename BGring::const_iterator itr = neighborCells[i].outer().begin();
	       itr != neighborCells[i].outer().end(); ++itr){
	    cerr << (*itr)
		 << "(" << (*itr).realx(mLowInner[0],mdxInner) 
		 << "," << (*itr).realy(mLowInner[1],mdxInner) << ")" << endl;
	  }
	  POLY_ASSERT(0);
	}
      }
	
      BGring boundaryRing = neighborCells[0].outer();
      POLY_ASSERT(boundaryRing.front() == boundaryRing.back());
      POLY_ASSERT(boundaryRing.size() > 2);
      
      // TODO: Make sure union-ing rings that share a common face results in 
      //       a closed boundary, has no repeated nodes, etc. etc.
      
      // Extract the boundary points from the union
      //
      // TODO: Check whether converting the PLC points back to doubles to compute
      //       the sub-tessellation gives a valid full tessellation after the
      //       cell adoption loop concludes
      vector<RealType> subPLCpoints;
      int nSides = 0;
      for (typename BGring::const_iterator itr = boundaryRing.begin();
	   itr != boundaryRing.end() - 1; ++itr, ++nSides) {
	subPLCpoints.push_back( (*itr).realx(mLowInner[0],mdxInner) );
	subPLCpoints.push_back( (*itr).realy(mLowInner[1],mdxInner) );
      }
      
      // Form the bounding PLC
      PLC<2, RealType> subPLC;
      subPLC.facets.resize(nSides, vector<int>(2) );
      for (unsigned ii = 0; ii < nSides; ++ii) {
	subPLC.facets[ii][0] = ii;
	subPLC.facets[ii][1] = (ii+1) % nSides;
      }
      
      this->computeCellRings(subpoints, subPLCpoints, subPLC, subCellRings, false);

//       // Blago!
//       cerr << endl << "Subpoints:" << endl;
//       for (int ii = 0; ii != subpoints.size(); ++ii){
// 	cerr << "  (" << subpoints[2*ii] << "," << subpoints[2*ii+1] << ")" << endl;
//       }
//       cerr << endl << "SubPLCpoints:" << endl;
//       for (int ii = 0; ii != subPLCpoints.size(); ++ii){
// 	cerr << "  (" << subPLCpoints[2*ii] << "," << subPLCpoints[2*ii+1] << ")" << endl;
//       }
//       cerr << endl << "SubCellRings:" << endl;
//       for (int ii = 0; ii != subCellRings.size(); ++ii) {
// 	cerr << "  Cell " << ii << endl;
// 	for (typename BGring::const_iterator itr = subCellRings[i].begin();
// 	     itr != subCellRings[i].end(); ++itr){
// 	  cerr << (*itr).realx(mLowInner[0],mdxInner) << " " 
// 	       << (*itr).realy(mLowInner[1],mdxInner) << " ";
// 	}
//       }
//       // Blago!

    }
    
    // We're only concerned with the cells in the sub-tessellation whose generators
    // are immediate neighbors of the orphaned chunk. These are the only cells which can
    // "adopt" the orphan based on a local Voronoi principle
    for (set<int>::const_iterator nbItr = orphanNeighbors.begin();
	 nbItr != orphanNeighbors.end(); ++nbItr){
      set<int>::iterator it = orphanNeighbors.find(*nbItr);
      POLY_ASSERT( it != orphanNeighbors.end() );
      int subIndex = distance(orphanNeighbors.begin(), it);
      POLY_ASSERT( subIndex < orphanNeighbors.size() );
      int thisIndex = *nbItr;
      POLY_ASSERT( thisIndex < numGenerators );
      BGring thisRing;
      
      if (orphanNeighbors.size() > 1){
	thisRing = subCellRings[subIndex];	    
	  
	// Simplify the resulting ring. Removes points that are within some minimum
	// distance to their neighbors. Setting distance = 1 merges ring elements
	// that are within one quantized mesh spacing. This essentially removes
	// repeated cell nodes having length-zero cell faces. An unfortunate
	// consequence of using a third-party lib to do all our unions/intersections
	BGring simplifiedRing;
	boost::geometry::simplify(thisRing, simplifiedRing, 2);
	thisRing = simplifiedRing;
	POLY_ASSERT(thisRing.size() > 2);
	POLY_ASSERT(thisRing.front() == thisRing.back());
      }	
	
      // If the orphan has only a single neighbor, just compute its union with
      // that neighbor's cell ring from the full tessellation
      else{
	thisRing = orphan;
      }
      
      // Union this new cell ring with the original cell ring from the full tessellation
      vector<BGring> unionRing;
      boost::geometry::union_( thisRing, cellRings[thisIndex], unionRing );
      if(unionRing.size() > 1){
	cerr << "Blago!" << endl << "Cell " << thisIndex
	     << " has more than one cell ring:" << endl;
	for( i=0; i<unionRing.size(); ++i){
	  cerr << endl << "Ring " << i << ":" << endl;
	  for (typename BGring::const_iterator itr = unionRing[i].begin();
	       itr != unionRing[i].end(); ++itr){
	    cerr << (*itr).realx(mLowInner[0],mdxInner) << " " 
		 << (*itr).realy(mLowInner[1],mdxInner) << " "
		 << (*itr) << endl;
	  }
	}
      }
      POLY_ASSERT(unionRing.size() == 1);
      thisRing = unionRing[0];
      
      // Simplify the final ring. 
      BGring simplifiedRing;
      boost::geometry::simplify(thisRing, simplifiedRing, 1);
      thisRing = simplifiedRing;
      POLY_ASSERT(thisRing.size() > 2);
      POLY_ASSERT(thisRing.front() == thisRing.back());

//       // Union may produce cell rings containing faces broken into multiple
//       // pieces by redundant hanging nodes. This manifests as three (or more)
//       // collinear nodes. We clip them as a final connection step
//       BGring finalRing;
//       bool collinear;
//       for (typename BGring::const_iterator itr = thisRing.begin()+1;
// 	   itr != thisRing.end()-1; ++itr) {
// 	collinear = geometry::collinear<2,CoordHash>(&(*(itr-1)).x,
// 						     &(*(itr  )).x,
// 						     &(*(itr+1)).x,
// 						     1);
// 	if(!collinear) boost::geometry::append(finalRing,*itr);
//       }
//       // Check the beginning/end point
//       collinear = geometry::collinear<2,CoordHash>(&(*(thisRing.end()-1  )).x,
// 						   &(*(thisRing.end()    )).x,
// 						   &(*(thisRing.begin()+1)).x,
// 						   1);
//       if(!collinear) boost::geometry::append(finalRing,*(thisRing.begin()));
//       POLY_ASSERT(finalRing.front() == finalRing.back());
//       thisRing = finalRing;
//       // TODO: Ensure the first and last entries of finalRing are the same
      
      cellRings[thisIndex] = thisRing;
    }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TriangleTessellator<double>;
}




// //------------------------------------------------------------------------------
// // An implementation of the map specialized for testing true/false.
// // If tested with a key that does not exist, it is initialized as false.
// //------------------------------------------------------------------------------
// template<typename Key, 
//          typename Comparator = std::less<Key> >
// class BoolMap: public std::map<Key, bool> {
// public:
//   BoolMap(): std::map<Key, bool>() {}
//   virtual ~BoolMap() {}
//   bool operator[](const Key& key) const {
//     typename std::map<Key, bool>::const_iterator itr = this->find(key);
//     if (itr == this->end()) return false;
//     return itr->second;
//   }
// };

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
