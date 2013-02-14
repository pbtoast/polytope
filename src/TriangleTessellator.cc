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

namespace {

///------------------------------------------------------------------------
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
//------------------------------------------------------------------------

//------------------------------------------------------------------------------
// An implementation of the map specialized for testing true/false.
// If tested with a key that does not exist, it is initialized as false.
//------------------------------------------------------------------------------
template<typename Key, 
         typename Comparator = std::less<Key> >
class BoolMap: public std::map<Key, bool> {
public:
  BoolMap(): std::map<Key, bool>() {}
  virtual ~BoolMap() {}
  bool operator[](const Key& key) const {
    typename std::map<Key, bool>::const_iterator itr = this->find(key);
    if (itr == this->end()) return false;
    return itr->second;
  }
};

// //------------------------------------------------------------------------------
// // Update the vector of thingies to unique indices.
// //------------------------------------------------------------------------------
// template<typename T>
// void
// addToVector(const unsigned i, std::vector<T>& vec) {
//    if( i > vec.size() )
//   const typename std::map<Key, int>::const_iterator itr = key2id.find(key);
//   int result;
//   if (itr == key2id.end()) {
//     result = key2id.size();
//     key2id[key] = result;
//   } else {
//     result = itr->second;
//   }
//   return result;
// }

//------------------------------------------------------------------------------
// Predicate to check if either element of a std::pair matches the given value.
//------------------------------------------------------------------------------
template<typename T>
struct MatchEitherPairValue {
  T mval;
  MatchEitherPairValue(const T& x): mval(x) {}
  bool operator()(const std::pair<T, T>& x) const { return (x.first == mval or x.second == mval); }
};

//------------------------------------------------------------------------------
// Get the set of map keys
//------------------------------------------------------------------------------
template<typename Key, typename T>
std::set<Key>
getMapKeys(std::map<Key, T>& mapIn) {
  typename std::set<Key> result;
  for (typename std::map<Key, T>::const_iterator itr = mapIn.begin();
       itr != mapIn.end(); ++itr){
    result.insert( itr->first );
  }
  return result;
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
      // orderedEdges.push_back(edges[0]);
      // nodes2edges[lastNode].erase(edges[0]);
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
// Counts the number of times we recursively call tessellate
//------------------------------------------------------------------------------
class RecursionCounter {
public:
   inline RecursionCounter() { ++count; }
   inline ~RecursionCounter() { --count; }
   inline operator int() const { return count; }
private:
   static int count;
};
int RecursionCounter::count = 0;


} // end anonymous namespace


//------------------------------------------------------------------------------
template<typename RealType>
TriangleTessellator<RealType>::
TriangleTessellator():
  Tessellator<2, RealType>() {
   mLow.assign(  2,  numeric_limits<RealType>::max() );
   mHigh.assign( 2, -numeric_limits<RealType>::max() );
   mdx = 0;
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
  
  POLY_ASSERT(!points.empty());

  // Make sure we're not modifying an existing tessellation.
  POLY_ASSERT(mesh.empty());

  const CoordHash coordMax = (1LL << 30); // numeric_limits<CoordHash>::max() >> 32U;
  const double degeneracy = 1.0e-12;
  
  // Use an empty PLCpoints vector for the time being
  std::vector<RealType> PLCpoints;
  triangulateio delaunay;
  computeDelaunay(points, PLCpoints, delaunay);

  const unsigned numGenerators = points.size()/2;
  
  //--------------------------------------------------------
  // Create the Voronoi tessellation from the triangulation.
  //--------------------------------------------------------

  // Create the Voronoi nodes from the list of triangles. Each triangle 
  // has 3 nodes p, q, r, and corresponds to a Voronoi node at (X,Y), say.
  // The Voronoi node is located at the center of the triangle, though things
  // get a little squirrely at boundaries.  On boundary edges we create
  // a vertex at the edge center.

  // Find the circumcenters of each triangle, and build the set of triangles
  // associated with each generator.
  // RealType  clow[2] = { numeric_limits<RealType>::max(),  numeric_limits<RealType>::max()};
  // RealType chigh[2] = {-numeric_limits<RealType>::max(), -numeric_limits<RealType>::max()};
  vector<RealPoint> circumcenters(delaunay.numberoftriangles);
  map<EdgeHash, vector<unsigned> > edge2tri;
  map<int, set<unsigned> > gen2tri;
  int pindex, qindex, rindex, i, j;
  EdgeHash pq, pr, qr;
  for (i = 0; i != delaunay.numberoftriangles; ++i) {
    pindex = delaunay.trianglelist[3*i  ];
    qindex = delaunay.trianglelist[3*i+1];
    rindex = delaunay.trianglelist[3*i+2];
    pq = internal::hashEdge(pindex, qindex);
    pr = internal::hashEdge(pindex, rindex);
    qr = internal::hashEdge(qindex, rindex);
    geometry::computeCircumcenter2d(&delaunay.pointlist[2*pindex],
                                    &delaunay.pointlist[2*qindex],
                                    &delaunay.pointlist[2*rindex],
                                    &circumcenters[i].x);
    gen2tri[pindex].insert(i);
    gen2tri[qindex].insert(i);
    gen2tri[rindex].insert(i);
    edge2tri[pq].push_back(i);
    edge2tri[pr].push_back(i);
    edge2tri[qr].push_back(i);
    // cerr << delaunay.pointlist[2*pindex]   << " " 
    //      << delaunay.pointlist[2*pindex+1] << " " 
    //      << delaunay.pointlist[2*qindex]   << " " 
    //      << delaunay.pointlist[2*qindex+1] << " " 
    //      << delaunay.pointlist[2*rindex]   << " " 
    //      << delaunay.pointlist[2*rindex+1] << " "
    //      << circumcenters[i].x             << " " 
    //      << circumcenters[i].y             << endl; 
    mLow[0]  = min(mLow[0], circumcenters[i].x);
    mLow[1]  = min(mLow[1], circumcenters[i].y);
    mHigh[0] = max(mHigh[0], circumcenters[i].x);
    mHigh[1] = max(mHigh[1], circumcenters[i].y);
  }

  // The circumcenters may all lie inside the convex hull of the
  // generators for an unbounded tessellation.
  for (i = 0; i != delaunay.numberofpoints; ++i){
     mLow[0]  = min(mLow[0] , delaunay.pointlist[2*i  ]);
     mLow[1]  = min(mLow[1] , delaunay.pointlist[2*i+1]);
     mHigh[0] = max(mHigh[0], delaunay.pointlist[2*i  ]);
     mHigh[1] = max(mHigh[1], delaunay.pointlist[2*i+1]);
  }
  POLY_ASSERT(circumcenters.size() == delaunay.numberoftriangles);
  POLY_ASSERT(mLow[0] < mHigh[0] and mLow[1] < mHigh[1]);
  
  // The bounding box which contains PLC, and all circumcenters and generators
  RealType cbox[2] = {mHigh[0] - mLow[0], mHigh[1] - mLow[1]};
  
  // The bounding circle onto which we project the "infinite" rays of the 
  // unbounded faces of the tessellation.
  const RealType rinf = 4.0*max(cbox[0], cbox[1]);
  const RealType cboxc[2] = {0.5*(mLow[0]+mHigh[0]), 0.5*(mLow[1]+mHigh[1])};

  // We resize mLow and boxsize so that the bounding box
  // contains the "infinite" sphere. mHigh is not really needed.
  mLow [0] = cboxc[0]-rinf;  mLow [1] = cboxc[1]-rinf;
  mHigh[0] = cboxc[0]+rinf;  mHigh[1] = cboxc[1]+rinf;
  const double cboxsize = 2.0*rinf;
  mdx = max(degeneracy, cboxsize/coordMax);
  
  // Map circumcenters and triangle indices to global id's
  map<IntPoint, int> circ2id;
  map<int, unsigned> tri2id;
  for (i = 0; i != delaunay.numberoftriangles; ++i){
    IntPoint ip(circumcenters[i].x, circumcenters[i].y,
                mLow[0], mLow[1], mdx);
    j = internal::addKeyToMap(ip, circ2id);
    tri2id[i] = j;
  }

  cerr << "Bounding box:"
       << "(" << mLow[0] << "," << mHigh[0] << ")x"
       << "(" << mLow[1] << "," << mHigh[1] << ")" << endl;
  cerr << "Box size    = " << cboxsize << endl;
  cerr << "spacing     = " << mdx << endl;
  
  // The exterior edges of the triangularization have "unbounded" rays, originating
  // at the circumcenter of the corresponding triangle and passing perpendicular to
  // the edge
    
  bool test;
  RealPoint ehat, test_point, tricent, pinf;
  map<EdgeHash, unsigned> edge2id;
  unsigned k, i1, i2;
  mesh.infNodes = vector<unsigned>(circ2id.size());
  for (map<EdgeHash, vector<unsigned> >::const_iterator edgeItr = edge2tri.begin();
       edgeItr != edge2tri.end(); ++edgeItr){
    const EdgeHash& edge = edgeItr->first;
    const vector<unsigned>& tris = edgeItr->second;
    if (tris.size() == 1){
      i = tris[0];
      POLY_ASSERT(i < delaunay.numberoftriangles);
      i1 = edge.first;
      i2 = edge.second;
       
      // Compute the triangle centroid (guaranteed pt inside the triangle)
      pindex = delaunay.trianglelist[3*i  ];
      qindex = delaunay.trianglelist[3*i+1];
      rindex = delaunay.trianglelist[3*i+2];
      geometry::computeTriangleCentroid2d(&delaunay.pointlist[2*pindex],
                                          &delaunay.pointlist[2*qindex],
                                          &delaunay.pointlist[2*rindex],
                                          &tricent.x);
       
      // Unit vector pointing normal to edge and through circumcenter
      ehat.x = -(delaunay.pointlist[2*i1+1] - delaunay.pointlist[2*i2+1]);
      ehat.y =  (delaunay.pointlist[2*i1  ] - delaunay.pointlist[2*i2  ]);
      geometry::unitVector<2, RealType>(&ehat.x);
       
      // Make sure the unit vector points outward
      copy(&delaunay.pointlist[2*i1], &delaunay.pointlist[2*i1] + 2, &test_point.x);
      test_point += ehat;
      if (orient2d(&delaunay.pointlist[2*i1], &delaunay.pointlist[2*i1+1], &tricent.x)*
          orient2d(&delaunay.pointlist[2*i2], &delaunay.pointlist[2*i2+1], &test_point.x) > 0.0){
         ehat *= -1.0;
      }
       
      // Get the intersection point along the "infinite" circumcircle
      test = geometry::rayCircleIntersection(&circumcenters[i].x,
                                             &ehat.x,
                                             cboxc,
                                             rinf,
                                             1.0e-10,
                                             &pinf.x);
      POLY_ASSERT(test);
      IntPoint ip(pinf.x, pinf.y, mLow[0], mLow[1], mdx);
      k = circ2id.size();
      j = internal::addKeyToMap(ip, circ2id);
      POLY_ASSERT(edge2id.find(edge) == edge2id.end());
      edge2id[edge] = j;
      if (k != circ2id.size()) mesh.infNodes.push_back(1);
    }
  }
  
  // Copy the quantized nodes to the final tessellation.
  const unsigned numNodes = circ2id.size();
  mesh.nodes.resize(2*numNodes);
  for (map<IntPoint, int>::const_iterator itr = circ2id.begin();
       itr != circ2id.end(); ++itr) {
    POLY_ASSERT(itr->second >= 0 and itr->second < numNodes);
    i = itr->second;
    mesh.nodes[2*i  ] = itr->first.realx(mLow[0], mdx);
    mesh.nodes[2*i+1] = itr->first.realy(mLow[1], mdx);
  }
    
  // The faces corresponding to each triangle edge
  unsigned ii, jj, iface;
  EdgeHash face;
  map<EdgeHash, int> face2id;
  unsigned numFaces = edge2tri.size();
  double orientation;
  internal::CounterMap<unsigned> faceCounter;
  vector<unsigned> faceVec(2);
  //mesh.faces = vector<vector<unsigned> >(numFaces);
  //mesh.faceCells = vector<vector<int> >(numFaces);
  mesh.cells = vector<vector<int> >(numGenerators);
  mesh.infFaces = vector<vector<unsigned> >(numGenerators);
  //vector<unsigned> mesh.infFaces = vector<unsigned>;
  for (map<int, set<unsigned> >::const_iterator genItr = gen2tri.begin();
       genItr != gen2tri.end(); ++genItr) {
    pindex = genItr->first;
    const set<unsigned>& tris = genItr->second;
    POLY_ASSERT(pindex < numGenerators);
    
    set<EdgeHash> meshEdges;
    for (std::set<unsigned>::const_iterator triItr = tris.begin();
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
        POLY_ASSERT(edge2tri[pq].size() == 2 and 
                    edge2tri[pq][0] == i or edge2tri[pq][1] == i);
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
        POLY_ASSERT(edge2tri[pr].size() == 2 and 
                    edge2tri[pr][0] == i or edge2tri[pr][1] == i);
        k = (edge2tri[pr][0] == i ? edge2tri[pr][1] : edge2tri[pr][0]);
        jj = tri2id[k];
        if (jj != ii) meshEdges.insert(internal::hashEdge(ii,jj));
      }
    }
    
    // Get the face sorted nodes
    const vector<unsigned> faceNodes = 
       computeSortedFaceNodes(vector<EdgeHash>(meshEdges.begin(), meshEdges.end()));
      
    // The ordered mesh nodes around a given generator
    POLY_ASSERT(faceNodes.size() > 2);
    for (i = 0; i != faceNodes.size(); ++i){
      i1 = faceNodes[i];
      i2 = faceNodes[(i+1) % faceNodes.size()];
      face  = internal::hashEdge(i1,i2);
        
      orientation = orient2d(&mesh.nodes[2*face.first], 
                             &mesh.nodes[2*face.second], 
                             &delaunay.pointlist[2*pindex]);
      iface = internal::addKeyToMap(face,face2id);
      ++faceCounter[iface];
        
      if ( i1 >= delaunay.numberoftriangles and i2 >= delaunay.numberoftriangles ){
        mesh.infFaces[pindex].push_back(face.first);
        mesh.infFaces[pindex].push_back(face.second);
      }

      POLY_ASSERT( faceCounter[iface] > 0 );
      if( faceCounter[iface] == 1 ){
        faceVec[0] = face.first; faceVec[1] = face.second;
        mesh.faces.push_back(faceVec);
        mesh.faceCells.resize(iface+1);
      }
        
      if( orientation > 0 ){
        mesh.cells[pindex].push_back(iface);
        mesh.faceCells[iface].push_back(pindex);
      }else{
        mesh.cells[pindex].push_back(~iface);
        mesh.faceCells[iface].push_back(~int(pindex));
      }
    }
  }
  
  // Blago!
  cerr << "Nodes:" << endl;
  for (i = 0; i < numNodes; ++i){
     cerr << "  (" << mesh.nodes[2*i] << "," << mesh.nodes[2*i+1] << ")" << endl;
  }
  cerr << "Faces:" << endl;
  for (i = 0; i < mesh.faces.size(); ++i){
     cerr << "  Face " << i << ":" << endl;
     for (j = 0; j < mesh.faces[i].size(); ++j){
        cerr << "    " << mesh.faces[i][j] << endl;
     }
  }
  cerr << "FaceCells:" << endl;
  for (i = 0; i < mesh.faceCells.size(); ++i){
     cerr << "  Face " << i << ":" << endl;
     for (j = 0; j < mesh.faceCells[i].size(); ++j){
        cerr << "    " << mesh.faceCells[i][j] << endl;
     }
  }
  cerr << "Cells:" << endl;
  for (i = 0; i < mesh.cells.size(); ++i){
     cerr << "  Cell " << i << ":" << endl;
     for (j = 0; j < mesh.cells[i].size(); ++j){
        cerr << "    " << mesh.cells[i][j] << endl;
     }
  }
  cerr << "InfNodes:" << endl;
  for (i = 0; i < mesh.infNodes.size(); ++i){
     cerr << "  " << mesh.infNodes[i] << endl;
  }
  cerr << "Cell-to-Nodes:" << endl;
  std::vector<std::set<unsigned> > cellToNodes = mesh.computeCellToNodes();
  for (i = 0; i < cellToNodes.size(); ++i){
     cerr << "  Cell " << i << ":" << endl;
     for (std::set<unsigned>::const_iterator itr = cellToNodes[i].begin();
          itr != cellToNodes[i].end(); ++itr ){
        cerr << "    " << *itr << endl;
     }
  }
  // Blago!


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
  RecursionCounter recursionDepth;
  POLY_ASSERT(!points.empty());

  // Make sure we're not modifying an existing tessellation.
  POLY_ASSERT(mesh.empty());

  const CoordHash coordMax = (1LL << 30); // numeric_limits<CoordHash>::max() >> 32U;
  const double degeneracy = 1.0e-12;
  
  // // Find the range of the generator points.
  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints = PLCpoints.size()/2;
  int i, j, k;
  
  RealType vertices[2*numPLCpoints];
  for(i=0; i<numPLCpoints; ++i) vertices[i] = PLCpoints[i];
  // for (i = 0; i != numGenerators; ++i) {
  //   low[0] = min(low[0], points[2*i]);
  //   low[1] = min(low[1], points[2*i+1]);
  //   high[0] = max(high[0], points[2*i]);
  //   high[1] = max(high[1], points[2*i+1]);
  // }
  
  for (i = 0; i != numPLCpoints; ++i) {
    mLow [0] = min(mLow [0], PLCpoints[2*i  ]);
    mLow [1] = min(mLow [1], PLCpoints[2*i+1]);
    mHigh[0] = max(mHigh[0], PLCpoints[2*i  ]);
    mHigh[1] = max(mHigh[1], PLCpoints[2*i+1]);
  }
  POLY_ASSERT(mLow[0] < mHigh[0] and mLow[1] < mHigh[1]);

  // Start by creating an unbounded tessellation
  tessellate(points, mesh);
  
  // RealType box[2] = {high[0] - low[0], 
  //                    high[1] - low[1]};
  // const double boxsize = 2.0*max(box[0], box[1]);

  // // No point attributes or markers.
  // in.numberofpointattributes = 0;
  // in.pointattributelist = 0; 
  // in.pointmarkerlist = 0; // No point markers.
  
  // // Define input points, including our false external generators.
  // if( geometry.empty() ){
  //    in.numberofpoints = numGenerators;
  //    in.pointlist = new RealType[2*in.numberofpoints];
  //    copy(points.begin(), points.end(), in.pointlist);
  //    in.numberofsegments = 0;
  // }else{
  //    in.numberofpoints = numGenerators + 4;
  //    in.pointlist = new RealType[2*in.numberofpoints];
  //    copy(points.begin(), points.end(), in.pointlist);
  //    const double xmin = 0.5*(low[0] + high[0]) - boxsize;
  //    const double ymin = 0.5*(low[1] + high[1]) - boxsize;
  //    const double xmax = 0.5*(low[0] + high[0]) + boxsize;
  //    const double ymax = 0.5*(low[1] + high[1]) + boxsize;
  //    in.pointlist[2*numGenerators  ] = xmin; in.pointlist[2*numGenerators+1] = ymin;
  //    in.pointlist[2*numGenerators+2] = xmax; in.pointlist[2*numGenerators+3] = ymin;
  //    in.pointlist[2*numGenerators+4] = xmax; in.pointlist[2*numGenerators+5] = ymax;
  //    in.pointlist[2*numGenerators+6] = xmin; in.pointlist[2*numGenerators+7] = ymax;
  //    // cerr << "Chose bounding range : (" << xmin << " " << ymin << ") (" << xmax << " " << ymax << ")" << endl;
  
  //    // Fill in Triangle's boundary info.  We use our imposed outer box of fake
  //    // generators.
  //    in.numberofsegments = 4;
  //    in.segmentlist = new int[2*in.numberofsegments];
  //    j = 0;
  //    for (i = 0; i != 4; ++i) {
  //       in.segmentlist[j++] = numGenerators + i;
  //       in.segmentlist[j++] = numGenerators + ((i + 1) % 4);
  //    }
  // }
  // in.segmentmarkerlist = 0;
  // in.numberofholes = 0;
  // in.holelist = 0;

  // // No regions.
  // in.numberofregions = 0;
  // in.regionlist = 0;

  // // Set up the structure for the triangulation.
  // delaunay.pointlist = 0;
  // delaunay.pointattributelist = 0;
  // delaunay.pointmarkerlist = 0;
  // delaunay.trianglelist = 0;
  // delaunay.triangleattributelist = 0;
  // delaunay.neighborlist = 0;
  // delaunay.segmentlist = 0;
  // delaunay.segmentmarkerlist = 0;
  // delaunay.edgelist = 0;
  // delaunay.edgemarkerlist = 0;
  // delaunay.holelist = 0;

  // // Do the triangulation. Switches pass to triangle are:
  // // -Q : Quiet (shaddap!), no output on the terminal except errors.
  // // -z : Indices are all numbered from zero.
  // // -e : Generates edges and places them in out.edgelist.
  // // -c : Generates convex hull and places it in out.segmentlist.
  // // -p : Uses the given PLC information.
  // if (geometry.empty()){
  //   triangulate((char*)"Qzec", &in, &delaunay, 0);
  // }else{
  //   triangulate((char*)"Qzep", &in, &delaunay, 0);
  // }
  // //triangulate((char*)"Qzep", &in, &delaunay, 0);

  // // Make sure we got something.
  // if (delaunay.numberoftriangles == 0)
  //   error("TriangleTessellator: Delauney triangulation produced 0 triangles!");
  // if (delaunay.numberofpoints != in.numberofpoints) {
  //   char err[1024];
  //   snprintf(err, 1024, "TriangleTessellator: Delauney triangulation produced %d triangles\n(%d generating points given)", 
  //            delaunay.numberofpoints, (int)numGenerators);
  //   error(err);
  // }

  // //--------------------------------------------------------
  // // Create the Voronoi tessellation from the triangulation.
  // //--------------------------------------------------------

  // // Create the Voronoi nodes from the list of triangles. Each triangle 
  // // has 3 nodes p, q, r, and corresponds to a Voronoi node at (X,Y), say.
  // // The Voronoi node is located at the center of the triangle, though things
  // // get a little squirrely at boundaries.  On boundary edges we create
  // // a vertex at the edge center.

  // // Find the circumcenters of each triangle, and build the set of triangles
  // // associated with each generator.
  // RealType  clow[2] = { numeric_limits<RealType>::max(),  numeric_limits<RealType>::max()};
  // RealType chigh[2] = {-numeric_limits<RealType>::max(), -numeric_limits<RealType>::max()};
  // internal::CounterMap<EdgeHash> edgeCounter;
  // vector<RealPoint> circumcenters(delaunay.numberoftriangles);
  // map<EdgeHash, vector<unsigned> > edge2tri;
  // map<int, set<unsigned> > gen2tri;
  // int k, pindex, qindex, rindex, iedge;
  // EdgeHash pq, pr, qr;
  // for (i = 0; i != delaunay.numberoftriangles; ++i) {
  //   pindex = delaunay.trianglelist[3*i  ];
  //   qindex = delaunay.trianglelist[3*i+1];
  //   rindex = delaunay.trianglelist[3*i+2];
  //   pq = internal::hashEdge(pindex, qindex);
  //   pr = internal::hashEdge(pindex, rindex);
  //   qr = internal::hashEdge(qindex, rindex);
  //   geometry::computeCircumcenter2d(&delaunay.pointlist[2*pindex],
  //                                   &delaunay.pointlist[2*qindex],
  //                                   &delaunay.pointlist[2*rindex],
  //                                   &circumcenters[i].x);
  //   if (pindex < numGenerators and qindex < numGenerators and rindex < numGenerators) {
  //     ++edgeCounter[pq];
  //     ++edgeCounter[pr];
  //     ++edgeCounter[qr];
  //   }
  //   gen2tri[pindex].insert(i);
  //   gen2tri[qindex].insert(i);
  //   gen2tri[rindex].insert(i);
  //   edge2tri[pq].push_back(i);
  //   edge2tri[pr].push_back(i);
  //   edge2tri[qr].push_back(i);
  //   cerr << delaunay.pointlist[2*pindex]   << " " 
  //        << delaunay.pointlist[2*pindex+1] << " " 
  //        << delaunay.pointlist[2*qindex]   << " " 
  //        << delaunay.pointlist[2*qindex+1] << " " 
  //        << delaunay.pointlist[2*rindex]   << " " 
  //        << delaunay.pointlist[2*rindex+1] << " "
  //        << circumcenters[i].x             << " " 
  //        << circumcenters[i].y             << endl; 
  //   clow[0] = min(clow[0], circumcenters[i].x);
  //   clow[1] = min(clow[1], circumcenters[i].y);
  //   chigh[0] = max(chigh[0], circumcenters[i].x);
  //   chigh[1] = max(chigh[1], circumcenters[i].y);
  // }
  
  // POLY_ASSERT(circumcenters.size() == delaunay.numberoftriangles);
  // POLY_ASSERT(clow[0] < chigh[0] and clow[1] < chigh[1]);
  // RealType cbox[2] = {chigh[0] - clow[0], 
  //                     chigh[1] - clow[1]};
  // const double cboxsize = 2.0*max(cbox[0], cbox[1]);
  // const double cdx = max(degeneracy, cboxsize/coordMax);

  // if( recursionDepth == 1 ){
  //    mLow.resize(2);     mHigh.resize(2);
  //    mLow[0] = clow[0];  mHigh[0] = chigh[0];
  //    mLow[1] = clow[1];  mHigh[1] = chigh[1];
  //    mdx = cdx;
  // }
  
  // // Build the polygon representing our boundaries.
  // BGpolygon boundary;
  // // realBGpolygon realBoundary;
  // // Copy the PLC provided boundary information into a Boost.Geometry polygon.
  // vector<IntPoint> boundaryPoints;
  // boundaryPoints.reserve(geometry.facets.size() + 1);
  // i = geometry.facets[0][0];
  // boundaryPoints.push_back(IntPoint(PLCpoints[2*i], PLCpoints[2*i+1],
  //                                   mLow[0], mLow[1], mdx));
  // // boost::geometry::append( realBoundary, realBGpoint(PLCpoints[2*i], PLCpoints[2*i+1]));

  // for (j = 0; j != geometry.facets.size(); ++j) {
  //    POLY_ASSERT(geometry.facets[j].size() == 2);
  //    i =  geometry.facets[j][1];
  //    boundaryPoints.push_back(IntPoint(PLCpoints[2*i], PLCpoints[2*i+1],
  //                                      mLow[0], mLow[1], mdx));
  //    // boost::geometry::append( realBoundary, realBGpoint(PLCpoints[2*i], PLCpoints[2*i+1]));
  // }
  // POLY_ASSERT(boundaryPoints.size() == geometry.facets.size() + 1);
  // POLY_ASSERT(boundaryPoints.front() == boundaryPoints.back());
  // boost::geometry::assign(boundary, BGring(boundaryPoints.begin(), boundaryPoints.end()));

  // // Add any interior holes.
  // const unsigned numHoles = geometry.holes.size();
  // if (numHoles > 0) {
  //    typename BGpolygon::inner_container_type& holes = boundary.inners();
  //    holes.resize(numHoles);
  //    // typename realBGpolygon::inner_container_type& realHoles = realBoundary.inners();
  //    // realHoles.resize(numHoles);
  //    for (k = 0; k != numHoles; ++k) {
  //       boundaryPoints = vector<IntPoint>();
  //       boundaryPoints.reserve(geometry.holes[k].size() + 1);
  //       i = geometry.holes[k][0][0];
  //       boundaryPoints.push_back(IntPoint(PLCpoints[2*i], PLCpoints[2*i+1],
  //                                         mLow[0], mLow[1], mdx));
  //       // boost::geometry::append( realHoles[k], realBGpoint(PLCpoints[2*i], PLCpoints[2*i+1]));
  //       for (j = 0; j != geometry.holes[k].size(); ++j) {
  //          POLY_ASSERT(geometry.holes[k][j].size() == 2);
  //          i =  geometry.holes[k][j][1];
  //          boundaryPoints.push_back(IntPoint(PLCpoints[2*i], PLCpoints[2*i+1],
  //                                            mLow[0], mLow[1], mdx));
  //          // boost::geometry::append( realHoles[k], realBGpoint(PLCpoints[2*i], PLCpoints[2*i+1]));
  //       }
  //       POLY_ASSERT(boundaryPoints.size() == geometry.holes[k].size() + 1);
  //       POLY_ASSERT(boundaryPoints.front() == boundaryPoints.back());
  //       // boost::geometry::assign(holes[k], BGring(boundaryPoints.begin(), boundaryPoints.end()));
  //    }
  // }

  // Quantize the PLCpoints
  std::vector<IntPoint> IntPLCPoints(numPLCpoints);
  for (i = 0; i < numPLCpoints; ++i){
    IntPLCPoints[i] = IntPoint( PLCpoints[2*i], PLCpoints[2*i+1],
				mLow[0], mLow[1], mdx );
  }

  // Generate the quantized boundary to handle boost intersections
  BGpolygon boundary;
  buildBoostBoundary(IntPLCPoints, geometry, boundary);
  
  // Walk each generator and build up it's unique nodes and faces.
  //mesh.cells.resize(numGenerators);
  IntPoint X, IntNode;
  bool inside;
  map<IntPoint, int> point2node;
  map<IntPoint, set<int> > point2neighbors;
  map<EdgeHash, int> edgeHash2id;
  map<int, vector<int> > edgeCells;
  map<int, BGring> cellRings;
  map<int, vector<BGring> > orphanage;
  for (i = 0; i != numGenerators; ++i) {
    vector<IntPoint> cellBoundary;
    for (vector<int>::const_iterator faceItr = mesh.cells[i].begin();
         faceItr != mesh.cells[i].end(); ++faceItr){
      const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
      POLY_ASSERT(iface < mesh.faceCells.size());
      POLY_ASSERT(mesh.faces[iface].size() == 2);
      const unsigned inode1 = *faceItr < 0 ? mesh.faces[iface][1] : mesh.faces[iface][0];
      const unsigned inode2 = *faceItr < 0 ? mesh.faces[iface][0] : mesh.faces[iface][1];
      IntNode = IntPoint(mesh.nodes[2*inode1  ],
                         mesh.nodes[2*inode1+1],
                         mLow[0], mLow[1], mdx);
      point2neighbors[IntNode].insert(i);
      cellBoundary.push_back(IntNode);
      cerr << "Adding (" << mesh.nodes[2*inode1] << "," << mesh.nodes[2*inode1+1] << ") " 
           << IntNode << endl;
      if( mesh.infNodes[inode1]==1 and mesh.infNodes[inode2]==1 ){
         // Check that segment connectig node1 and node2 doesn't intersect inner
         // bounding radius.
         //    If it does: get an intermediate point at the outer "infinite" radius
         //                in between node1 and node2, quantize it, and add it
         //                to the cell ring
      }
    }
    cellBoundary.push_back( cellBoundary[0] );  // Close the ring
    boost::geometry::assign(cellRings[i], BGring(cellBoundary.begin(), cellBoundary.end()));
    boost::geometry::correct(cellRings[i]);
  // }

  // vector<set<unsigned> > cellToNodes = mesh.computeCellToNodes();
  // POLY_ASSERT(cellToNodes.size() == numGenerators);
  // for (i = 0; i != numGenerators; ++i) {

  //    // Add the circumcenters as points for the cell.
  //    set<IntPoint> cellPointSet;
  //    for (set<unsigned>::const_iterator nodeItr = cellToNodes[i].begin();
  //         nodeItr != cellToNodes[i].end(); ++nodeItr) {
  //       cellPointSet.insert(IntPoint(mesh.nodes[2*(*nodeItr)], mesh.nodes[2*(*nodeItr)+1],
  //                                    clow[0], clow[1], cdx));
  //    }
  //    POLY_ASSERT2(cellPointSet.size() >= 3, cellPointSet.size());

  //    // // Build the convex hull of the cell points.
  //    // vector<double> cellPointCoords;
  //    // for (j = 0; j != cellPoints.size(); ++j) {
  //    //   cellPointCoords.push_back(cellPoints[j].x);
  //    //   cellPointCoords.push_back(cellPoints[j].y);
  //    // }
  //    // POLY_ASSERT(cellPointCoords.size() == 2*cellPoints.size());
  //    // PLC<2, double> hull = convexHull_2d<double>(cellPointCoords, low, dx);
  //    // POLY_ASSERT(hull.facets.size() >= 3);
  //    // POLY_ASSERT(hull.facets[0][0] < cellPoints.size());
  //    // vector<RealPoint> ringPoints;
  //    // ringPoints.push_back(cellPoints[hull.facets[0][0]]);
  //    // for (j = 0; j != hull.facets.size(); ++j) {
  //    //   POLY_ASSERT(hull.facets[j].size() == 2);
  //    //   POLY_ASSERT(hull.facets[j][1] < cellPoints.size());
  //    //   ringPoints.push_back(cellPoints[hull.facets[j][1]]);
  //    // }
  //    // cellRings[i] = BGring(ringPoints.begin(), ringPoints.end());

  //    // Build the convex hull of the cell points.
  //    BGmulti_point mpoints(cellPointSet.begin(), cellPointSet.end());
  //    // for (typename set<IntPoint>::const_iterator itr = cellPointSet.begin();
  //    //      itr != cellPointSet.end();
  //    //      ++itr) mpoints.push_back(RealPoint(itr->realx(mLow[0], mdx),
  //    //                                         itr->realy(mLow[1], mdx)));
  //    boost::geometry::convex_hull(mpoints, cellRings[i]);

    // Intersect with the boundary to get the bounded cell.
    // Since for complex boundaries this may return more than one polygon, we find
    // the one that contains the generator.
    vector<BGring> cellIntersections;
    boost::geometry::intersection(boundary, cellRings[i], cellIntersections);
    if (cellIntersections.size() == 0) {
      cerr << points[2*i] << " " << points[2*i+1] << endl 
           << boost::geometry::dsv(cellRings[i]) << endl
         //<< boost::geometry::dsv(mpoints) << endl
           << boost::geometry::dsv(boundary) << endl;
    }
    POLY_ASSERT(cellIntersections.size() > 0);
    if (cellIntersections.size() == 1) {
      cellRings[i] = cellIntersections[0];
    } else {
      X = IntPoint(points[2*i], points[2*i+1], mLow[0], mLow[1], mdx);
      k = cellIntersections.size();
      for (j = 0; j != cellIntersections.size(); ++j) {
        inside = boost::geometry::within(X, cellIntersections[j]);
        if( inside )  k = j;
        else          orphanage[i].push_back( cellIntersections[j] );
      }
      POLY_ASSERT(k < cellIntersections.size());
      cellRings[i] = cellIntersections[k];
    }
    
    // Construct the map from node points to neighboring cells
    // for (typename BGring::const_iterator itr = cellRings[i].begin();
    //       itr != cellRings[i].end() - 1; ++itr) {
    //     // const IntPoint& pX1 = *itr;
    //     // const IntPoint& pX2 = *(itr + 1);
    //     // POLY_ASSERT(*itr != *(itr + 1));
    //     // j = internal::addKeyToMap(pX1, point2node);
    //     // k = internal::addKeyToMap(pX2, point2node);
    //     // POLY_ASSERT(j != k);
    //     // iedge = internal::addKeyToMap(internal::hashEdge(j, k), edgeHash2id);
    //     // edgeCells[iedge].push_back(j < k ? i : ~i);
    //     // mesh.cells[i].push_back(j < k ? iedge : ~iedge);

    //     std::map<IntPoint, set<int> >::iterator p2nItr = point2neighbors.find(*itr);
    //     point2neighbors[*itr].insert(i);
    //  }
    // POLY_ASSERT(mesh.cells[i].size() >= 3);
  }
  // POLY_ASSERT(edgeCells.size() == edgeHash2id.size());

  // Build the map from cell to set of neighboring cells  
  // map<int, set<int> > neighbors;
  // for (i = 0; i != numGenerators; ++i){
  //    for (typename BGring::const_iterator itr = cellRings[i].begin();
  //         itr != cellRings[i].end() - 1; ++itr) {
  //       neighbors[i].insert( point2neighbors[*itr].begin(), point2neighbors[*itr].end() );
  //    }
  //    neighbors[i].erase(i);
  // }
  
  
  // // Blago!
  // for (i = 0; i < numGenerators; ++i){
  //    cerr << "Cell " << i << " has generator neighbors";
  //    for (std::set<int>::const_iterator nbItr = neighbors[i].begin();
  //         nbItr != neighbors[i].end(); ++nbItr){
  //       cerr << " " << *nbItr;
  //    }
  //    cerr << endl << orphanage.size() << endl;;
  // }
  // // Blago!
  
  
  //*********************** Begin Adoption Algorithm ************************
  // Intersecting cells with the boundary has created orphaned cell pieces. ("Won't 
  // somebody please think of the children!?") Make sub-tessellations using the 
  // generators neighboring the orphaned pieces. Construct a PLC boundary for the 
  // sub-tessellation by using the geometry obtained by union-ing the orphan 
  // with its neighboring cells. The way we compute cell neighbors should ensure
  // that the union gives a contiguous geometry with no holes
  if( orphanage.size() > 0 and recursionDepth == 1 ){
     cerr << orphanage.size() << " orphaned cells" << endl;
     for (std::map<int, std::vector<BGring> >::const_iterator orphanItr = orphanage.begin();
          orphanItr != orphanage.end(); ++orphanItr){
        int parent = orphanItr->first;
        cerr << "Parent cell: " << parent << ":" << endl;
        for (unsigned iorphan = 0; iorphan != orphanItr->second.size(); ++iorphan){
           BGring orphan = orphanItr->second[iorphan];

           // Build the neighboring cells of the orphaned chunk
           std::set<int> orphanNeighbors;
           //std::set<int> orphanNeighborhood;
           for (typename BGring::const_iterator pointItr = orphan.begin();
                pointItr != orphan.end() - 1; ++pointItr) {
              std::map<IntPoint, std::set<int> >::iterator it = point2neighbors.find( *pointItr );
              if (it != point2neighbors.end()){
                 std::set<int> neighborSet = it->second;
                 for (std::set<int>::const_iterator setItr = neighborSet.begin();
                      setItr != neighborSet.end(); ++setItr){
                    orphanNeighbors.insert(*setItr);
                    //orphanNeighborhood.insert(*setItr);
                    //orphanNeighborhood.insert(neighbors[*setItr].begin(), neighbors[*setItr].end());
                 }
              }
           }
           POLY_ASSERT( orphanNeighbors.size() > 0 );
           //POLY_ASSERT( orphanNeighborhood.size() > 0 );
        
        
           // // Blago!
           // cerr << "Orphaned piece has neighbors";
           // for( std::set<int>::const_iterator iii = orphanNeighbors.begin();
           //      iii != orphanNeighbors.end(); ++iii){
           //   cerr << " " << *iii;
           // }
           // cerr << endl << "and neighborhood";
           // for( std::set<int>::const_iterator iii = orphanNeighborhood.begin();
           //      iii != orphanNeighborhood.end(); ++iii){
           //   cerr << " " << *iii;
           // }
           // cerr << endl;
           // // Blago!
        
        
           // If the orphan only has a single neighbor, we can skip a lot of work.
           // No need to tessellate - simply union the orphan with its neighbor cell.
           Tessellation<2, RealType> submesh;
           if (orphanNeighbors.size() > 1){
          
              // Compute the sub-tessellation from orphan's neighboring points. Union the
              // orphan and its immediate neighbors to get the tessellation boundary
              std::vector<RealType> subpoints;
          
              BGmulti_polygon neighborCells;
              createBGUnion(orphan,neighborCells);          
              for (std::set<int>::const_iterator nbItr = orphanNeighbors.begin();
                   nbItr != orphanNeighbors.end(); ++nbItr){
                 subpoints.push_back( points[2*(*nbItr)  ] );
                 subpoints.push_back( points[2*(*nbItr)+1] );
                 createBGUnion(cellRings[*nbItr],neighborCells);
              }
              POLY_ASSERT2( neighborCells.size() > 0, "Union produced empty set!" );
              if (neighborCells.size() > 1){
                 cerr << "Blago!" << endl;
                 for (i = 0; i != neighborCells.size(); ++i){
                    cerr << "Sub-polygon " << i << " in the union has bounding ring" << endl;
                    for (typename BGring::const_iterator itr = neighborCells[i].outer().begin();
                         itr != neighborCells[i].outer().end(); ++itr){
                       cerr << (*itr)
                            << "(" << (*itr).realx(mLow[0],mdx) 
                            << "," << (*itr).realy(mLow[1],mdx) << ")" << endl;
                    }
                    POLY_ASSERT(0);
                 }
              }
          
              BGring boundaryRing = neighborCells[0].outer();
          
              // TODO: Make sure union-ing rings that share a common face results in 
              //       a closed boundary, has no repeated nodes, etc. etc.
                    
              // Extract the boundary points from the union
              //
              // TODO: Check whether converting the PLC points back to doubles to compute
              //       the sub-tessellation gives a valid full tessellation after the
              //       cell adoption loop concludes
              std::vector<RealType> subPLCpoints;
              int nSides = 0;
              for (typename BGring::const_iterator itr = boundaryRing.begin();
                   itr != boundaryRing.end() - 1; ++itr, ++nSides) {
                 subPLCpoints.push_back( (*itr).realx(mLow[0],mdx) );
                 subPLCpoints.push_back( (*itr).realy(mLow[1],mdx) );
              }

              // Form the bounding PLC
              PLC<2, RealType> subPLC;
              subPLC.facets.resize(nSides, std::vector<int>(2) );
              for (i = 0; i < nSides; ++i) {
                 subPLC.facets[i][0] = i;
                 subPLC.facets[i][1] = (i+1) % nSides;
              }

              tessellate(subpoints,subPLCpoints,subPLC,submesh);
           }
        
           // We're only concerned with the cells in the sub-tessellation whose generators
           // are immediate neighbors of the orphaned chunk. These are the only cells which can
           // adopt the orphan based on the Voronoi principle of ownership based on "closeness"
           for (std::set<int>::const_iterator nbItr = orphanNeighbors.begin();
                nbItr != orphanNeighbors.end(); ++nbItr){
              std::set<int>::iterator it = orphanNeighbors.find(*nbItr);
              POLY_ASSERT( it != orphanNeighbors.end() );
              int subIndex = std::distance(orphanNeighbors.begin(), it);
              POLY_ASSERT( subIndex < orphanNeighbors.size() );
              int thisIndex = *it;
              POLY_ASSERT( thisIndex < numGenerators );
          
              BGring thisRing;
              if( orphanNeighbors.size() > 1 ){
                 // Walk the ordered nodes of the cell and build its boost.geometry ring
                 std::vector<IntPoint> cellBoundary;
                 for (std::vector<int>::const_iterator faceItr = submesh.cells[subIndex].begin();
                      faceItr != submesh.cells[subIndex].end(); ++faceItr){
                    const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
                    POLY_ASSERT(iface < submesh.faceCells.size());
                    POLY_ASSERT(submesh.faces[iface].size() == 2);
                    const unsigned inode = *faceItr < 0 ? submesh.faces[iface][1] : submesh.faces[iface][0];
                    cellBoundary.push_back(IntPoint(submesh.nodes[2*inode  ], 
                                                    submesh.nodes[2*inode+1],
                                                    mLow[0], mLow[1], mdx));
                 }
                 cellBoundary.push_back( cellBoundary[0] );  // Close the ring
                 boost::geometry::assign(thisRing, BGring(cellBoundary.begin(), 
                                                          cellBoundary.end()) );
            
            
                 // // Blago!
                 // cerr << endl << "Cell " << thisIndex << endl;
                 // cerr << endl << "SUBMESH CELL:" << endl;
                 // for (typename BGring::const_iterator itr = thisRing.begin();
                 //      itr != thisRing.end(); ++itr){
                 //    cerr << (*itr).realx(mLow[0],mdx) << " " 
                 //         << (*itr).realy(mLow[1],mdx) << endl;
                 // }
                 // for (typename BGring::const_iterator itr = thisRing.begin();
                 //      itr != thisRing.end(); ++itr){
                 //    cerr << *itr << endl;
                 // }
                 // // Blago!

            
                 // Simplify the resulting ring. Removes points that are within some minimum
                 // distance to their neighbors. Setting distance = 1 merges ring elements
                 // that are within one quantized mesh spacing. This essentially removes
                 // repeated cell nodes having length-zero cell faces.
                 BGring simplifiedRing;
                 boost::geometry::simplify(thisRing, simplifiedRing, 1);
                 thisRing = simplifiedRing;
              }
        
              // If the orphan has only a single neighbor, just compute its union with
              // that neighbor's cell ring from the full tessellation
              else{
                 thisRing = orphan;
              }
             
              // Union this new cell ring with the original cell ring from the full tessellation
              std::vector<BGring> unionRing;
              boost::geometry::union_( thisRing, cellRings[thisIndex], unionRing );
              POLY_ASSERT(unionRing.size() == 1);
              thisRing = unionRing[0];

              // Simplify the final ring. 
              BGring simplifiedRing;
              boost::geometry::simplify(thisRing, simplifiedRing, 1);
              thisRing = simplifiedRing;

        
              // // Blago!
              // cerr << endl << "Cell " << thisIndex << endl;
              // cerr << endl << "FINAL SUBMESH CELL:" << endl;
              // for (typename BGring::const_iterator itr = thisRing.begin();
              //      itr != thisRing.end(); ++itr){
              //    cerr << (*itr).realx(mLow[0],mdx) << " " 
              //         << (*itr).realy(mLow[1],mdx) << endl;
              // }
              // for (typename BGring::const_iterator itr = thisRing.begin();
              //      itr != thisRing.end(); ++itr){
              //    cerr << *itr << endl;
              // }
              // // Blago!


              cellRings[thisIndex] = thisRing;
           }
        }
     }
  }
  //*********************** End Adoption Algorithm ************************

  // Now build the unique mesh nodes and cell info.
  int iedge;
  for (i = 0; i != numGenerators; ++i) { 
     for (typename BGring::const_iterator itr = cellRings[i].begin();
          itr != cellRings[i].end() - 1;
          ++itr) {
        const IntPoint& pX1 = *itr;
        const IntPoint& pX2 = *(itr + 1);
        POLY_ASSERT(*itr != *(itr + 1));
        j = internal::addKeyToMap(pX1, point2node);
        k = internal::addKeyToMap(pX2, point2node);
        POLY_ASSERT(j != k);
        iedge = internal::addKeyToMap(internal::hashEdge(j, k), edgeHash2id);
        edgeCells[iedge].push_back(j < k ? i : ~i);
        mesh.cells[i].push_back(j < k ? iedge : ~iedge);
        // cerr << "Cell " << i << " adding edge " << iedge << " : " << pX1 << " " << pX2 << " : (" 
        //      << pX1.realx(mLow[0], mdx) << " " << pX1.realy(mLow[1], mdx) << ") ("
        //      << pX2.realx(mLow[0], mdx) << " " << pX2.realy(mLow[1], mdx) << ")" << endl;
     }
     POLY_ASSERT(mesh.cells[i].size() >= 3);
  }
  POLY_ASSERT(edgeCells.size() == edgeHash2id.size());


  // Fill in the mesh nodes.
  RealType node[2];
  mesh.nodes = vector<RealType>(2*point2node.size());
  for (typename map<IntPoint, int>::const_iterator itr = point2node.begin();
       itr != point2node.end();
       ++itr) {
     const IntPoint& p = itr->first;
     i = itr->second;
     POLY_ASSERT(i < mesh.nodes.size()/2);
     node[0] = p.realx(mLow[0],mdx);
     node[1] = p.realy(mLow[1],mdx);
    
     cerr << "Node (" << node[0] << "," << node[1] << ")" << endl;

     // Check if nodes are inside boundary (either bounding box or PLC, if defined)
     //bool inside = boost::geometry::within(realBGpoint(node[0],node[1]), realBoundary);
     bool inside = within(node, numPLCpoints, vertices, geometry);
     
     // if( !inside ){
     //   RealType result[2];
     //   RealType dist = nearestPoint( node, numPLCpoints, vertices, geometry, result );
     //   // Check the node has not moved more than 2.5 quantized mesh spacings. NOTE: this is
     //   // not a sharp estimate. Theoreticallly, the distance ought to be at most sqrt(2)*cdx, 
     //   // but nodes will fail this strict of a test.
     //   POLY_ASSERT2( dist < 2.5*mdx,
     //                 dist << " " << 2.5*mdx << " : (" << node[0] << " " << node[1] 
     //                 << ") (" << result[0] << " " << result[1] << ")\n" << geometry);
     //   node[0] = result[0];
     //   node[1] = result[1];
     // }
    
     //POLY_ASSERT( node[0] >= low[0] and node[0] <= high[0] );
     //POLY_ASSERT( node[1] >= low[1] and node[1] <= high[1] );
     mesh.nodes[2*i]   = node[0];
     mesh.nodes[2*i+1] = node[1];
  }
  
  // Fill in the mesh faces.
  mesh.faces = vector<vector<unsigned> >(edgeHash2id.size());
  for (typename map<EdgeHash, int>::const_iterator itr = edgeHash2id.begin();
       itr != edgeHash2id.end();
       ++itr) {
     const EdgeHash& ehash = itr->first;
     i = itr->second;
     POLY_ASSERT(i < mesh.faces.size());
     POLY_ASSERT(mesh.faces[i].size() == 0);
     mesh.faces[i].push_back(ehash.first);
     mesh.faces[i].push_back(ehash.second);
  }

  // Fill in the mesh faceCells.
  mesh.faceCells = vector<vector<int> >(mesh.faces.size());
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
  
  // // Clean up.
  // trifree((VOID*)delaunay.pointlist);
  // trifree((VOID*)delaunay.pointmarkerlist);
  // trifree((VOID*)delaunay.trianglelist);
  // trifree((VOID*)delaunay.edgelist);
  // trifree((VOID*)delaunay.edgemarkerlist);
  // trifree((VOID*) delaunay.segmentlist);
  // trifree((VOID*) delaunay.segmentmarkerlist);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TriangleTessellator<double>;




//------------------------------------------------------------------------------
//PRIVATE STUFF:
//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeDelaunay(const vector<RealType>& points,
                const vector<RealType>& PLCpoints,
                triangulateio& delaunay) const 
{
  triangulateio in;
   
  // Find the range of the generator points.
  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints = PLCpoints.size()/2;
  RealType vertices[2*numPLCpoints];
  for( unsigned i=0; i<PLCpoints.size(); ++i) vertices[i] = PLCpoints[i];
  RealType  low[2] = { numeric_limits<RealType>::max(),  numeric_limits<RealType>::max()};
  RealType high[2] = {-numeric_limits<RealType>::max(), -numeric_limits<RealType>::max()};
  int i, j;
  for (i = 0; i != numGenerators; ++i) {
    low[0] = min(low[0], points[2*i]);
    low[1] = min(low[1], points[2*i+1]);
    high[0] = max(high[0], points[2*i]);
    high[1] = max(high[1], points[2*i+1]);
  }
  for (i = 0; i != numPLCpoints; ++i) {
    low[0] = min(low[0], PLCpoints[2*i]);
    low[1] = min(low[1], PLCpoints[2*i+1]);
    high[0] = max(high[0], PLCpoints[2*i]);
    high[1] = max(high[1], PLCpoints[2*i+1]);
  }
  POLY_ASSERT(low[0] < high[0] and low[1] < high[1]);
  RealType box[2] = {high[0] - low[0], 
                     high[1] - low[1]};
  const double boxsize = 2.0*max(box[0], box[1]);

  // If no PLC boundary, simply add the generators
  if( numPLCpoints == 0 ){
    in.numberofpoints = numGenerators;
    in.pointlist = new RealType[2*in.numberofpoints];
    copy(points.begin(), points.end(), in.pointlist);
    in.numberofsegments = 0;
  }
  // Include false external generators if creating a bounded tessellation
  else{
    in.numberofpoints = numGenerators + 4;
    in.pointlist = new RealType[2*in.numberofpoints];
    copy(points.begin(), points.end(), in.pointlist);
    const double xmin = 0.5*(low[0] + high[0]) - boxsize;
    const double ymin = 0.5*(low[1] + high[1]) - boxsize;
    const double xmax = 0.5*(low[0] + high[0]) + boxsize;
    const double ymax = 0.5*(low[1] + high[1]) + boxsize;
    in.pointlist[2*numGenerators  ] = xmin; in.pointlist[2*numGenerators+1] = ymin;
    in.pointlist[2*numGenerators+2] = xmax; in.pointlist[2*numGenerators+3] = ymin;
    in.pointlist[2*numGenerators+4] = xmax; in.pointlist[2*numGenerators+5] = ymax;
    in.pointlist[2*numGenerators+6] = xmin; in.pointlist[2*numGenerators+7] = ymax;
    
    // Fill in Triangle's boundary info. Use an imposed outer box of fake generators.
    in.numberofsegments = 4;
    in.segmentlist = new int[2*in.numberofsegments];
    j = 0;
    for (i = 0; i != 4; ++i) {
      in.segmentlist[j++] = numGenerators + i;
      in.segmentlist[j++] = numGenerators + ((i + 1) % 4);
    }
  }
  
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
  if (numPLCpoints == 0 ) triangulate((char*)"Qzec", &in, &delaunay, 0);
  else                    triangulate((char*)"Qzep", &in, &delaunay, 0);

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
  if( numPLCpoints > 0 ) delete [] in.segmentlist;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
something(triangulateio& delaunay) const {
  // Create the Voronoi nodes from the list of triangles. Each triangle 
  // has 3 nodes p, q, r, and corresponds to a Voronoi node at (X,Y), say.
  // The Voronoi node is located at the center of the triangle, though things
  // get a little squirrely at boundaries.  On boundary edges we create
  // a vertex at the edge center.

  // Find the circumcenters of each triangle, and build the set of triangles
  // associated with each generator.
  RealType  clow[2] = { numeric_limits<RealType>::max(),  numeric_limits<RealType>::max()};
  RealType chigh[2] = {-numeric_limits<RealType>::max(), -numeric_limits<RealType>::max()};
  vector<RealPoint> circumcenters(delaunay.numberoftriangles);
  map<EdgeHash, vector<unsigned> > edge2tri;
  map<int, set<int> > gen2tri;
  int i, pindex, qindex, rindex;
  EdgeHash pq, pr, qr;
  for (i = 0; i != delaunay.numberoftriangles; ++i) {
    pindex = delaunay.trianglelist[3*i];
    qindex = delaunay.trianglelist[3*i + 1];
    rindex = delaunay.trianglelist[3*i + 2];
    pq = internal::hashEdge(pindex, qindex);
    pr = internal::hashEdge(pindex, rindex);
    qr = internal::hashEdge(qindex, rindex);
    geometry::computeCircumcenter2d(&delaunay.pointlist[2*pindex],
                                    &delaunay.pointlist[2*qindex],
                                    &delaunay.pointlist[2*rindex],
                                    &circumcenters[i].x);
    gen2tri[pindex].insert(i);
    gen2tri[qindex].insert(i);
    gen2tri[rindex].insert(i);
    edge2tri[pq].push_back(i);
    edge2tri[pr].push_back(i);
    edge2tri[qr].push_back(i);
    clow[0] = min(clow[0], circumcenters[i].x);
    clow[1] = min(clow[1], circumcenters[i].y);
    chigh[0] = max(chigh[0], circumcenters[i].x);
    chigh[1] = max(chigh[1], circumcenters[i].y);
  }
  
  POLY_ASSERT(circumcenters.size() == delaunay.numberoftriangles);
  POLY_ASSERT(clow[0] < chigh[0] and clow[1] < chigh[1]);
  RealType cbox[2] = {chigh[0] - clow[0], 
                      chigh[1] - clow[1]};
  const double cboxsize = 2.0*max(cbox[0], cbox[1]);
  const double cdx = max(degeneracy, cboxsize/coordMax);
}

}
