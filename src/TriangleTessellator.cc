//------------------------------------------------------------------------
// TriangleTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <set>
#include <list>
#include <map>
#include "float.h"
#include <limits>

#include "polytope.hh" // Pulls in ASSERT and TriangleTessellator.hh.
#include "convexHull_2d.hh"
#include "Polygon.hh"  // For polygon projections, etc.

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

//------------------------------------------------------------------------
// This function computes the circumcenter of a triangle with vertices
// A = (Ax, Ay), B = (Bx, By), and C = (Cx, Cy), and places the result 
// in X.
//------------------------------------------------------------------------
void 
computeCircumcenter(double* A, double* B, double* C, double* X)
{
  // This solution was taken from Wikipedia's entry:
  // http://en.wikipedia.org/wiki/Circumscribed_circle
  double D = 2.0*(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1]));
  X[0] = ((A[0]*A[0] + A[1]*A[1])*(B[1]-C[1]) + (B[0]*B[0] + B[1]*B[1])*(C[1]-A[1]) + 
          (C[0]*C[0] + C[1]*C[1])*(A[1]-B[1]))/D;
  X[1] = ((A[0]*A[0] + A[1]*A[1])*(C[0]-B[0]) + (B[0]*B[0] + B[1]*B[1])*(A[0]-C[0]) + 
          (C[0]*C[0] + C[1]*C[1])*(B[0]-A[0]))/D;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
// Bound a position (X) to be inside the triangle (A,B,C).  If X is 
// outside, it is projected to the midpoint of the closest edge.
//------------------------------------------------------------------------
int
_boundPointProjector(double*A, double* B, double* X) {
  const int checkAB = orient2d(A, B, X);
  if (checkAB <= 0) {
    cerr << "Projecting point " << X[0] << " " << X[1];
    X[0] = 0.5*(A[0] + B[0]);
    X[1] = 0.5*(A[1] + B[1]);
    cerr << " to " << X[0] << " " << X[1] << endl;
    return true;
  }
  return false;
}

void
boundPoint(double* A, double* B, double* C, double* X) {
  if (_boundPointProjector(A, B, X)) return;
  if (_boundPointProjector(B, C, X)) return;
  if (_boundPointProjector(C, A, X)) return;
}

//------------------------------------------------------------------------------
// An implementation of the map specialized to help constructing counters.
// This thing just overloads the index operator to start the count at zero
// for new key values.
//------------------------------------------------------------------------------
template<typename Key, 
         typename Comparator = std::less<Key> >
class CounterMap: public std::map<Key, unsigned> {
public:
  CounterMap(): std::map<Key, unsigned>() {}
  virtual ~CounterMap() {}
  unsigned& operator[](const Key& key) {
    typename std::map<Key, unsigned>::iterator itr = this->find(key);
    if (itr == this->end()) {
      std::map<Key, unsigned>::operator[](key) = 0U;
      itr = this->find(key);
    }
    ASSERT(itr != this->end());
    return itr->second;
  }
};

//------------------------------------------------------------------------------
// An implementation of the map specialized for our edge counting algorithm.
// This changes the following things:
//  this->flip(Key) : if Key is present, flip value.  If Key is new, set to true.
//  this->operator[Key]: if Key exists, return value.  Otherwise throw error!
//------------------------------------------------------------------------------
template<typename Key, 
         typename Comparator = std::less<Key> >
class BoolMap: public std::map<Key, bool> {
public:
  BoolMap(): std::map<Key, bool>() {}
  virtual ~BoolMap() {}
  bool operator[](const Key& key) const {
    typename std::map<Key, bool>::const_iterator itr = this->find(key);
    ASSERT(itr != this->end());
    return itr->second;
  }
  void flip(const Key& key) {
    typename std::map<Key, bool>::iterator itr = this->find(key);
    if (itr == this->end()) {
      this->insert(make_pair(key, true));
    } else {
      itr->second = !(itr->second);
    }
  }
};

//------------------------------------------------------------------------------
// Hash two node indices uniquely to represent an edge.
//------------------------------------------------------------------------------
std::pair<int, int>
hashEdge(const int i, const int j) {
  ASSERT(i != j);
  return i < j ? std::make_pair(i, j) : std::make_pair(j, i);
}

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
// Predicate to compare points based on distance from an origin.
//------------------------------------------------------------------------------
template<typename Point>
struct ComparePointDistances: public std::binary_function<Point, Point, bool> {
  const Point& origin;
  ComparePointDistances(const Point& x): origin(x) {}
  bool operator()(const Point& lhs, const Point& rhs) const {
    typedef typename Point::CoordType Coord;
    const Coord dx1 = lhs.x - origin.x;
    const Coord dy1 = lhs.y - origin.y;
    const Coord dx2 = rhs.x - origin.x;
    const Coord dy2 = rhs.y - origin.y;
    return ((dx1*dx1 + dy1*dy1) < (dx2*dx2 + dy2*dy2));
  }
};

//------------------------------------------------------------------------------
// Update the map of points to unique indices.
//------------------------------------------------------------------------------
template<typename Point>
int
addMeshNode(const Point& pX, map<Point, int>& point2node) {
  const typename map<Point, int>::const_iterator itr = point2node.find(pX);
  int result;
  if (itr == point2node.end()) {
    result = point2node.size();
    point2node[pX] = result;
  } else {
    result = itr->second;
  }
  return result;
}

//------------------------------------------------------------------------------
// Compute the intermediate point along the line segment (e1, e2) that is 
// equidistant from e1 and a third point p.
//------------------------------------------------------------------------------
double square(const double& a) { return a*a; }

void
computeEquidistantPoint(double* e1, double* e2, double* p, double* X) {
  double ehatx = e2[0] - e1[0];
  double ehaty = e2[1] - e1[1];
  const double e = sqrt(ehatx*ehatx + ehaty*ehaty);
  if (e < 1.0e-10) {
    X[0] = 0.5*(e1[0] + e2[0]);
    X[1] = 0.5*(e1[1] + e2[1]);
  } else {
    ehatx /= e;
    ehaty /= e;
    const double ax = p[0] - e1[0];
    const double ay = p[1] - e1[1];
    const double d = sqrt(square(ehatx*ax) + square(ehaty*ay));
    const double b = min(e, 0.5*(ax*ax + ay*ay)/d);
    X[0] = e1[0] + ehatx*b;
    X[1] = e1[1] + ehaty*b;
  }
}

} // end anonymous namespace

//------------------------------------------------------------------------------
template<typename RealType>
TriangleTessellator<RealType>::
TriangleTessellator():
  Tessellator<2, RealType>()
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
TriangleTessellator<RealType>::
~TriangleTessellator() 
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           Tessellation<2, RealType>& mesh) const 
{
  // Find the range of the input.
  RealType  low[2] = { numeric_limits<RealType>::max(),  numeric_limits<RealType>::max()};
  RealType high[2] = {-numeric_limits<RealType>::max(), -numeric_limits<RealType>::max()};
  unsigned i, j, n = points.size();
  for (i = 0; i != n; ++i)
  {
    j = i % 2;
    low[j] = min(low[j], points[i]);
    high[j] = max(high[j], points[i]);
  }
  ASSERT(low[0] < high[0] and low[1] < high[1]);
  RealType box[2] = {high[0] - low[0], 
                     high[1] - low[1]};
  const double dx = 1.0e-10*max(box[0], box[1]);

  // Build the convex hull of the input points.
  PLC<2, RealType> hull = convexHull_2d(points, low, dx);

  // Now tessellate in the hull.
  tessellate(points, hull, mesh);
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
  // First do a null boundary tessellation.
  this->tessellate(points, mesh);
  const unsigned numGens = points.size()/2;

  // Now extrapolate pseudo-generators from the convex hull to the bounding
  // box.
  unsigned i, j;
  for (i = 0; i != mesh.convexHull.facets.size(); ++i)
  {
    j = mesh.convexHull.facets[i][0];
  }

}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const PLC<2, RealType>& geometry,
           Tessellation<2, RealType>& mesh) const 
{
  ASSERT(!points.empty());

  // Make sure we're not modifying an existing tessellation.
  ASSERT(mesh.empty());

  typedef std::pair<int, int> EdgeHash;

  triangulateio in, delaunay;

  // Define input points.
  const unsigned numGenerators = points.size()/2;
  in.numberofpoints = numGenerators;
  in.pointlist = new RealType[points.size()];
  copy(points.begin(), points.end(), in.pointlist);

  // No point attributes or markers.
  in.numberofpointattributes = 0;
  in.pointattributelist = 0; 
  in.pointmarkerlist = 0; // No point markers.

  // Segments and/or holes.
  if (geometry.empty())
  {
    in.numberofsegments = 0;
    in.segmentlist = 0;
    in.segmentmarkerlist = 0;
    in.numberofholes = 0;
    in.holelist = 0;
  }
  else
  {
    in.numberofsegments = geometry.facets.size();
    in.segmentlist = new int[2*in.numberofsegments];
    int s = 0;
    for (int f = 0; f < geometry.facets.size(); ++f)
    {
      in.segmentlist[s++] = geometry.facets[f][0];
      in.segmentlist[s++] = geometry.facets[f][1];
    }
    in.segmentmarkerlist = 0;
    in.numberofholes = geometry.holes.size();
    in.holelist = new double[2*in.numberofholes];
    copy(geometry.holes.begin(), geometry.holes.end(), in.holelist);
  }

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
  if (geometry.empty())
  {
    delaunay.segmentlist = 0;
    delaunay.segmentmarkerlist = 0;
  }
  else
  {
    delaunay.numberofsegments = geometry.facets.size();
    delaunay.segmentlist = new int[2*delaunay.numberofsegments];
    delaunay.segmentmarkerlist = new int[delaunay.numberofsegments];
  }
  delaunay.edgelist = 0;
  delaunay.edgemarkerlist = 0;
  delaunay.holelist = 0;

  // Do the triangulation. Switches pass to triangle are:
  // -Q : Quiet (shaddap!), no output on the terminal except errors.
  // -z : Indices are all numbered from zero.
  // -e : Generates edges and places them in out.edgelist.
  // -c : Generates convex hull and places it in out.segmentlist.
  // -p : Uses the given PLC information.
  if (geometry.empty())
    triangulate((char*)"Qzec", &in, &delaunay, 0);
  else
    triangulate((char*)"Qzep", &in, &delaunay, 0);

  // Make sure we got something.
  if (delaunay.numberoftriangles == 0)
    error("TriangleTessellator: Delauney triangulation produced 0 triangles!");
  if (delaunay.numberofpoints != numGenerators)
  {
    char err[1024];
    snprintf(err, 1024, "TriangleTessellator: Delauney triangulation produced %d triangles\n(%d generating points given)", 
             delaunay.numberofpoints, (int)numGenerators);
    error(err);
  }

  // Transfer the convex hull data and build a Polygon representing 
  // the hull.
  mesh.convexHull.facets.resize(delaunay.numberofsegments);
  vector<double> hullVertices(2*delaunay.numberofsegments);
  RealType  low[2] = { numeric_limits<RealType>::max(),  numeric_limits<RealType>::max()};
  RealType high[2] = {-numeric_limits<RealType>::max(), -numeric_limits<RealType>::max()};
  for (int i = 0; i < delaunay.numberofsegments; ++i)
  {
    mesh.convexHull.facets[i].resize(2);
    mesh.convexHull.facets[i][0] = delaunay.segmentlist[2*i];
    mesh.convexHull.facets[i][1] = delaunay.segmentlist[2*i+1];
    hullVertices[2*i]   = points[2*delaunay.segmentlist[2*i]];
    hullVertices[2*i+1] = points[2*delaunay.segmentlist[2*i]+1];
    low[0] = min(low[0], hullVertices[2*i]);
    low[1] = min(low[1], hullVertices[2*i + 1]);
    high[0] = max(high[0], hullVertices[2*i]);
    high[1] = max(high[1], hullVertices[2*i + 1]);
  }
  Polygon<double> convexHull(hullVertices);

  //--------------------------------------------------------
  // Create the Voronoi tessellation from the triangulation.
  //--------------------------------------------------------

  // Create the Voronoi nodes from the list of triangles. Each triangle 
  // has 3 nodes p, q, r, and corresponds to a Voronoi node at (X,Y), say.
  // The Voronoi node is located at the center of the triangle, though things
  // get a little squirrely at boundaries.  On boundary edges we create
  // a vertex at the edge center.

  // Find the circumcenters of each triangle.  
  CounterMap<EdgeHash> edgeCounter;
  map<Point2<uint64_t>, set<int> > point2tri;
  vector<double> triCircumCenters(2*delaunay.numberoftriangles);
  double X[2];
  int i, j, k, pindex, qindex, rindex;
  for (i = 0; i != delaunay.numberoftriangles; ++i)
  {
    pindex = delaunay.trianglelist[3*i];
    qindex = delaunay.trianglelist[3*i + 1];
    rindex = delaunay.trianglelist[3*i + 2];
    ++edgeCounter[hashEdge(pindex, qindex)];
    ++edgeCounter[hashEdge(qindex, rindex)];
    ++edgeCounter[hashEdge(rindex, pindex)];
    computeCircumcenter(&delaunay.pointlist[2*pindex],
                        &delaunay.pointlist[2*qindex],
                        &delaunay.pointlist[2*rindex],
                        &triCircumCenters[2*i]);
    low[0] = min(low[0], triCircumCenters[2*i]);
    low[1] = min(low[1], triCircumCenters[2*i+1]);
    high[0] = max(high[0], triCircumCenters[2*i]);
    high[1] = max(high[1], triCircumCenters[2*i+1]);
    cerr << "circumcenter : " << triCircumCenters[2*i] << " " << triCircumCenters[2*i+1] << endl;
  }

  ASSERT(triCircumCenters.size() == 2*delaunay.numberoftriangles);
  ASSERT(low[0] < high[0] and low[1] < high[1]);
  RealType box[2] = {high[0] - low[0], 
                     high[1] - low[1]};
  const double dx = 1.0e-10*max(box[0], box[1]);

  // Flag any generators on the edge of the tessellation.
  vector<bool> exteriorGenerators(numGenerators, false);
  list<EdgeHash> exteriorEdges;
  for (typename CounterMap<EdgeHash>::const_iterator itr = edgeCounter.begin();
       itr != edgeCounter.end();
       ++itr) {
    ASSERT(itr->second == 1 or itr->second == 2);
    if (itr->second == 1) {
      i = itr->first.first;
      j = itr->first.second;
      exteriorGenerators[i] = true;
      exteriorGenerators[j] = true;
      exteriorEdges.push_back(itr->first);
    }
  }

  // Connect the exterior edges into an ordered Polygon representing the boundary.
  const unsigned numBoundaryPoints = exteriorEdges.size();
  vector<RealType> boundaryPoints;
  vector<EdgeHash> boundaryEdgeOrder;
  boundaryEdgeOrder.push_back(exteriorEdges.front());
  exteriorEdges.pop_front();
  i = boundaryEdgeOrder.front().first;
  boundaryPoints.reserve(2*numBoundaryPoints);
  copy(&points[2*i], &points[2*i+2], back_inserter(boundaryPoints));
  i = boundaryEdgeOrder.front().second;
  copy(&points[2*i], &points[2*i+2], back_inserter(boundaryPoints));
  while (exteriorEdges.size() > 0) {
    list<EdgeHash>::iterator itr = find_if(exteriorEdges.begin(), exteriorEdges.end(),
                                           MatchEitherPairValue<int>(i));
    ASSERT(itr != exteriorEdges.end());
    i = (itr->first == i ? itr->second : itr->first);
    boundaryEdgeOrder.push_back(*itr);
    exteriorEdges.erase(itr);
    if (exteriorEdges.size() > 0) copy(&points[2*i], &points[2*i+2], back_inserter(boundaryPoints));
  }
  ASSERT(boundaryEdgeOrder.size() == numBoundaryPoints);
  ASSERT(MatchEitherPairValue<int>(boundaryEdgeOrder.front().first)(boundaryEdgeOrder.back()));
  ASSERT(boundaryPoints.size() == 2*numBoundaryPoints);
  Polygon<RealType> boundary(boundaryPoints);

  // Walk the triangles again, and take care of any circumcenters that are outside
  // the allowed surface.
  // We also compute the integer representations of the mesh node coordinates, and
  // build up the mesh nodes as these integer versions.
  int e1, e2, iedge;
  Point2<uint64_t> pX;
  double X1[2], X2[2];
  map<Point2<uint64_t>, int> point2node;
  map<int, set<int> > gen2nodes;
  vector<vector<Point2<uint64_t> > > boundaryEdgePoints(numBoundaryPoints);
  for (i = 0; i != delaunay.numberoftriangles; ++i)
  {
    pindex = delaunay.trianglelist[3*i];
    qindex = delaunay.trianglelist[3*i + 1];
    rindex = delaunay.trianglelist[3*i + 2];

    // if (boundary.queryPoint(&triCircumCenters[2*i]) == Polygon<double>::POINT_OUTSIDE)
    // {
    //   cerr << "Outside circumcenter : " << triCircumCenters[2*i] << " " << triCircumCenters[2*i+1] << endl;
    //   cerr << "Triangle generator edges : (" << delaunay.pointlist[2*pindex] << " " << delaunay.pointlist[2*pindex+1] << ") ("
    //        << delaunay.pointlist[2*qindex] << " " << delaunay.pointlist[2*qindex+1] << ") ("
    //        << delaunay.pointlist[2*rindex] << " " << delaunay.pointlist[2*rindex+1] << ")" << endl;

    //   // This circumcenter is outside the tessellation boundary, and needs to be truncated.
    //   // It will turn into two new nodes along the mesh boundary.
    //   const double inside_orient = orient2d(&delaunay.pointlist[2*pindex], 
    //                                         &delaunay.pointlist[2*qindex], 
    //                                         &delaunay.pointlist[2*rindex]);
    //   if (orient2d(&delaunay.pointlist[2*pindex], 
    //                &delaunay.pointlist[2*qindex], 
    //                &triCircumCenters[2*i])*inside_orient < 0.0) {
    //     e1 = pindex; e2 = qindex; j = rindex;
    //     cerr << "Option pq : " << orient2d(&delaunay.pointlist[2*pindex], 
    //                                        &delaunay.pointlist[2*qindex], 
    //                                        &triCircumCenters[2*i]) << " " << inside_orient << endl;
    //   } else if (orient2d(&delaunay.pointlist[2*qindex], 
    //                       &delaunay.pointlist[2*rindex], 
    //                       &triCircumCenters[2*i])*inside_orient < 0.0) {
    //     e1 = qindex; e2 = rindex; j = pindex;
    //     cerr << "Option qr : " << orient2d(&delaunay.pointlist[2*qindex], 
    //                                        &delaunay.pointlist[2*rindex], 
    //                                        &triCircumCenters[2*i]) << " " << inside_orient << endl;
    //   } else {
    //     ASSERT(orient2d(&delaunay.pointlist[2*rindex], 
    //                     &delaunay.pointlist[2*pindex], 
    //                     &triCircumCenters[2*i])*inside_orient < 0.0);
    //     e1 = rindex; e2 = pindex; j = qindex;
    //     cerr << "Option rp : " << orient2d(&delaunay.pointlist[2*rindex], 
    //                                        &delaunay.pointlist[2*pindex], 
    //                                        &triCircumCenters[2*i]) << " " << inside_orient << endl;
    //   }

    //   X1[0] = 0.5*(delaunay.pointlist[2*e1]     + delaunay.pointlist[2*j]);
    //   X1[1] = 0.5*(delaunay.pointlist[2*e1 + 1] + delaunay.pointlist[2*j + 1]);
    //   X2[0] = 0.5*(delaunay.pointlist[2*e2]     + delaunay.pointlist[2*j]);
    //   X2[1] = 0.5*(delaunay.pointlist[2*e2 + 1] + delaunay.pointlist[2*j + 1]);

    //   // Bump the initial points a bit to get them off of the legs of the initial triangle.
    //   X1[0] += 1.0e-10*(triCircumCenters[2*i]   - X1[0]);
    //   X1[1] += 1.0e-10*(triCircumCenters[2*i+1] - X1[1]);
    //   X2[0] += 1.0e-10*(triCircumCenters[2*i]   - X2[0]);
    //   X2[1] += 1.0e-10*(triCircumCenters[2*i+1] - X2[1]);

    //   // First intersection.
    //   iedge = boundary.closestIntersection(X1, &triCircumCenters[2*i], X);
    //   cerr << "Intersect edge " << iedge << endl;
    //   pX = Point2<uint64_t>(X[0] - low[0], X[1] - low[1], dx);
    //   k = addMeshNode(pX, point2node);
    //   gen2nodes[e1].insert(k);
    //   gen2nodes[j].insert(k);
    //   ASSERT(iedge >= 0 and iedge < numBoundaryPoints);
    //   boundaryEdgePoints[iedge].push_back(pX);
    //   cerr << "         Projecting to " << X[0] << " " << X[1] << endl;

    //   // Second intersection.
    //   iedge = boundary.closestIntersection(X2, &triCircumCenters[2*i], X);
    //   cerr << "Intersect edge " << iedge << endl;
    //   pX = Point2<uint64_t>(X[0] - low[0], X[1] - low[1], dx);
    //   k = addMeshNode(pX, point2node);
    //   gen2nodes[e2].insert(k);
    //   gen2nodes[j].insert(k);
    //   ASSERT(iedge >= 0 and iedge < numBoundaryPoints);
    //   boundaryEdgePoints[iedge].push_back(pX);
    //   cerr << "         Projecting to " << X[0] << " " << X[1] << endl;

    // } else {
      // This circumcenter is on the interior.
      pX = Point2<uint64_t>(triCircumCenters[2*i] - low[0],
                            triCircumCenters[2*i+1] - low[1], dx);
      k = addMeshNode(pX, point2node);
      gen2nodes[pindex].insert(k);
      gen2nodes[qindex].insert(k);
      gen2nodes[rindex].insert(k);
    // }
  }

  // // Look for any remaining edges that need to be split.
  // for (iedge = 0; iedge != numBoundaryPoints; ++iedge) {
  //   i = boundaryEdgeOrder[iedge].first;
  //   j = boundaryEdgeOrder[iedge].second;
  //   cerr << "Bounding edge (" << points[2*i] << " " << points[2*i+1] << ") (" << points[2*j] << " " << points[2*j+1] << ") : " << boundaryEdgePoints[iedge].size() << endl;

  //   if (boundaryEdgePoints[iedge].size() == 0) {
  //     // This edge does not have any points inserted, so split it in half between the
  //     // two generator endpoints.
  //     pX = Point2<uint64_t>(0.5*(points[2*i]   + points[2*j])   - low[0],
  //                           0.5*(points[2*i+1] + points[2*j+1]) - low[1],
  //                           dx);
  //     k = addMeshNode(pX, point2node);
  //     cerr << "Adding mesh node " << k << " @ " << pX << " to generators " << i << " " << j << endl;
  //     gen2nodes[i].insert(k);
  //     gen2nodes[j].insert(k);

  //   } else {
  //     // This edge has already had boundary points added, so figure out which ones 
  //     // should be associated with each generator endpoint.
  //     pX = Point2<uint64_t>(points[2*i] - low[0], points[2*i+1] - low[1], dx);
  //     sort(boundaryEdgePoints[iedge].begin(), boundaryEdgePoints[iedge].end(),
  //          ComparePointDistances<Point2<uint64_t> >(pX));
  //     k = addMeshNode(boundaryEdgePoints[iedge].front(), point2node);
  //     gen2nodes[i].insert(k);
  //     k = addMeshNode(boundaryEdgePoints[iedge].back(), point2node);
  //     gen2nodes[j].insert(k);
  //   }
  // }

  // Create mesh nodes for each exterior generator.
  for (i = 0; i != numGenerators; ++i) {
    if (exteriorGenerators[i]) {
      pX = Point2<uint64_t>(points[2*i] - low[0], points[2*i+1] - low[1], dx);
      k = addMeshNode(pX, point2node);
      gen2nodes[i].insert(k);
    }
  }

  // Build the unique mesh nodes.
  mesh.nodes = vector<double>(2*point2node.size());
  for (typename map<Point2<uint64_t>, int>::const_iterator itr = point2node.begin();
       itr != point2node.end();
       ++itr) {
    i = itr->second;
    ASSERT(2*i < mesh.nodes.size());
    mesh.nodes[2*i]   = itr->first.realx(low[0], dx);
    mesh.nodes[2*i+1] = itr->first.realy(low[1], dx);
  }

  // Sort the vertices of each generator counter-clockwise.
  unsigned nv, i0, i1;
  vector<uint64_t> vertices;
  vector<unsigned> vindices;
  PLC<2, uint64_t> hull;
  mesh.cells.resize(numGenerators);
  EdgeHash ehash;
  map<EdgeHash, int> ehash2face;
  uint64_t ilow[2] = {0U, 0U}, imax = numeric_limits<uint64_t>::max()/4;
  double norm, nhat[2];
  for (i = 0; i != numGenerators; ++i)
  {
    cerr << "Generator @ " << points[2*i] << " " << points[2*i + 1] << endl;
    nv = gen2nodes[i].size();
    vindices = vector<unsigned>(gen2nodes[i].begin(), gen2nodes[i].end());
    vertices = vector<uint64_t>();
    vertices.reserve(2*nv);
    low[0] = numeric_limits<RealType>::max();
    low[1] = numeric_limits<RealType>::max();
    high[0] = -numeric_limits<RealType>::max();
    high[1] = -numeric_limits<RealType>::max();
    X[0] = 0.0;
    X[1] = 0.0;
    for (j = 0; j != nv; ++j)
    {
      k = vindices[j];
      ASSERT(k < mesh.nodes.size()/2);
      low[0] = min(low[0], mesh.nodes[2*k]);
      low[1] = min(low[1], mesh.nodes[2*k + 1]);
      high[0] = max(high[0], mesh.nodes[2*k]);
      high[1] = max(high[1], mesh.nodes[2*k + 1]);
      X[0] += mesh.nodes[2*k];
      X[1] += mesh.nodes[2*k + 1];
      cerr << "   --> mesh node : " << mesh.nodes[2*k] << " " << mesh.nodes[2*k+1] << endl;
    }
    cerr << "  low/high : (" << low[0] << " " << low[1] << ") (" << high[0] << " " << high[1] << ")" << endl;
    ASSERT(low[0] < high[0] and low[1] < high[1]);
    X[0] /= nv;
    X[1] /= nv;
    for (j = 0; j != nv; ++j)
    {
      k = vindices[j];
      ASSERT(k < mesh.nodes.size()/2);
      nhat[0] = mesh.nodes[2*k    ] - X[0];
      nhat[1] = mesh.nodes[2*k + 1] - X[1];
      norm = sqrt(nhat[0]*nhat[0] + nhat[1]*nhat[1]);
      ASSERT(norm > 0.0);
      nhat[0] /= norm;
      nhat[1] /= norm;
      vertices.push_back(uint64_t(nhat[0]*imax + 2.1*imax));
      vertices.push_back(uint64_t(nhat[1]*imax + 2.1*imax));
      // vertices.push_back(uint64_t((mesh.nodes[2*k    ] - low[0])/box[0]*numeric_limits<uint64_t>::max()/2));
      // vertices.push_back(uint64_t((mesh.nodes[2*k + 1] - low[1])/box[1]*numeric_limits<uint64_t>::max()/2));
      // cerr << "Adding vertex : " << k << " (" << mesh.nodes[2*k] << " " << mesh.nodes[2*k + 1] << ") (" << vertices[2*k] << " " << vertices[2*k + 1] << ")" << endl;
    }
    ASSERT(vertices.size() == 2*nv);
    hull = convexHull_2d(vertices, ilow, uint64_t(1));
    if (!(hull.facets.size() == nv)) {
      cerr << "Blago : " << hull.facets.size() << " " << nv << endl;
      for (j = 0; j != nv; ++j) {
        k = vindices[j];
        cerr << "  (" << vertices[2*j] << " " << vertices[2*j+1] << ")    ("
             << mesh.nodes[2*k] << " " << mesh.nodes[2*k+1] << ")" << endl;
      }
    }
    cerr << "Building cell " << i << " :";
    for (j = 0; j != nv; ++j) cerr << " (" << vindices[hull.facets[j][0]] << " " << vindices[hull.facets[j][1]] << ")";
    cerr << endl;
    ASSERT(hull.facets.size() == nv);
    mesh.cells[i].reserve(nv);
    for (j = 0; j != nv; ++j) {
      ASSERT(hull.facets[j].size() == 2);
      ASSERT(hull.facets[j][0] < nv and hull.facets[j][1] < nv);
      i0 = vindices[hull.facets[j][0]];
      i1 = vindices[hull.facets[j][1]];
      ehash = hashEdge(i0, i1);
      // cerr << "Checking for edge hash : " << ehash.first << " " << ehash.second << endl;
      if (ehash2face.find(ehash) == ehash2face.end()) {
        k = mesh.faces.size();
        vector<unsigned> thpt(2);
        thpt[0] = i0;
        thpt[1] = i1;
        mesh.faces.push_back(thpt);
        mesh.faceCells.push_back(vector<unsigned>());
        ehash2face[ehash] = k;
      } else {
        k = ehash2face[ehash];
        ASSERT((i0 == mesh.faces[k][1] and i1 == mesh.faces[k][0]) or
               (i0 == mesh.faces[k][0] and i1 == mesh.faces[k][1]));
        if (i0 == mesh.faces[k][1] and i1 == mesh.faces[k][0]) k = ~k;
      }
      ASSERT(mesh.faces.size() == mesh.faceCells.size());
      // cerr << "Adding face " << k << " to cell " << i << endl;
      mesh.cells[i].push_back(k);
      mesh.faceCells[k < 0 ? ~k : k].push_back(i);
      ASSERT(mesh.faceCells[k < 0 ? ~k : k].size() <= 2);
    }
  }

//   // The Voronoi node is located at the center of the circle that contains 
//   // p, q, and r. 
//   vector<vector<int> > cellNodes(delaunay.numberofpoints);

//   // A table of Voronoi nodes falling on Delaunay edges. This prevents
//   // us from double-counting nodes.
//   map<pair<int, int>, int> nodesOnEdges; 
//   for (int i = 0; i < delaunay.numberoftriangles; ++i)
//   {
//     // Coordinates of the triangle's vertices p, q, r.
//     int pindex = delaunay.trianglelist[3*i];
//     double p[2];
//     p[0] = delaunay.pointlist[2*pindex];
//     p[1] = delaunay.pointlist[2*pindex+1];

//     int qindex = delaunay.trianglelist[3*i+1];
//     double q[2];
//     q[0] = delaunay.pointlist[2*qindex],
//     q[1] = delaunay.pointlist[2*qindex+1];

//     int rindex = delaunay.trianglelist[3*i+2];
//     double r[2];
//     r[0] = delaunay.pointlist[2*rindex],
//     r[1] = delaunay.pointlist[2*rindex+1];

//     // // Get the circumcenter.
//     // double X[2];
//     // computeCircumcenter(p, q, r, X);

//     // Compute the center of the triangle.
//     double X[2];
//     X[0] = (p[0] + q[0] + r[0])/3.0;
//     X[1] = (p[1] + q[1] + r[1])/3.0;

//     // If the point is outside the convex hull, project it to 
//     // an edge.
//     if (convexHull.queryPoint(X) == Polygon<double>::POINT_OUTSIDE)
//     {
//       double Y[2];
//       convexHull.project(X, Y);
// //cout << "Projected (" << X[0] << ", " << X[1] << ", r = " << sqrt(X[0]*X[0]+X[1]*X[1]) << " -> (" << Y[0] << ", " << Y[1] << ", r = " << sqrt(Y[0]*Y[0]+Y[1]*Y[1]) << ")\n";
//       X[0] = Y[0];
//       X[1] = Y[1];
//     }

//     // If this node lies on an edge, add it to the table.
//     // This uses Jonathan Shewchuck's fast geometry predicate orient2d().
//     if (orient2d(p, q, X) == 0.0)
//     {
//       pair<int, int> key(min(pindex, qindex), max(pindex, qindex));
//       map<pair<int, int>, int>::const_iterator iter = nodesOnEdges.find(key);
//       if (iter != nodesOnEdges.end())
//       {
//         // Add this existing node to the missing cell and proceed.
//         cellNodes[rindex].push_back(iter->second);
//         continue;
//       }
//       else
//         nodesOnEdges[key] = i;
//     }
//     else if (orient2d(p, r, X) == 0.0) 
//     {
//       pair<int, int> key(min(pindex, rindex), max(pindex, rindex));
//       map<pair<int, int>, int>::const_iterator iter = nodesOnEdges.find(key);
//       if (iter != nodesOnEdges.end())
//       {
//         // Add this existing node to the missing cell and proceed.
//         cellNodes[qindex].push_back(iter->second);
//         continue;
//       }
//       else
//         nodesOnEdges[key] = i;
//     }
//     else if (orient2d(q, r, X) == 0.0)
//     {
//       pair<int, int> key(min(qindex, rindex), max(qindex, rindex));
//       map<pair<int, int>, int>::const_iterator iter = nodesOnEdges.find(key);
//       if (iter != nodesOnEdges.end())
//       {
//         // Add this existing node to the missing cell and proceed.
//         cellNodes[pindex].push_back(iter->second);
//         continue;
//       }
//       else
//         nodesOnEdges[key] = i;
//     }

//     // Assign the node coordinate.
// //cout << "Computed node " << mesh.nodes.size()/2 << " at (" << X[0] << ", " << X[1] << ")" << endl;
// //cout << " for cells at (" << p[0] << ", " << p[1] << "), (" << q[0] << ", " << q[1] << "), (" << r[0] << ", " << r[1] << ")\n";
// //cout << " for cells " << pindex << ", " << qindex << ", " << rindex << endl;
//     mesh.nodes.push_back(RealType(X[0]));
//     mesh.nodes.push_back(RealType(X[1]));

//     // Associate this node with its Voronoi cells.
//     cellNodes[pindex].push_back(i);
//     cellNodes[qindex].push_back(i);
//     cellNodes[rindex].push_back(i);
//   }

//   // Compute the Voronoi faces between cells. These have a one-to-one
//   // correspondence with Delaunay edges of triangles unless we have placed 
//   // nodes along edges.
//   mesh.cells.resize(delaunay.numberofpoints);
//   for (int i = 0; i < delaunay.numberofedges; ++i)
//   {
//     int cell1 = delaunay.edgelist[2*i],
//         cell2 = delaunay.edgelist[2*i+1];

//     // Have we placed any nodes along this edge? If so, there's no 
//     // face here.
//     pair<int, int> key(min(cell1, cell2), max(cell1, cell2));
//     map<pair<int, int>, int>::const_iterator iter = nodesOnEdges.find(key);
//     if (iter != nodesOnEdges.end())
// //{
// //cout << "cells " << cell1 << ", " << cell2 << " have a node on their edge." << endl;
//       continue;
// //}

//     // Hook the cells up to the faces.
//     mesh.faces.push_back(vector<unsigned>());
//     mesh.faceCells.push_back(vector<unsigned>());
//     mesh.faceCells.back().reserve(2);
//     mesh.faceCells.back().push_back(cell1);
//     mesh.faceCells.back().push_back(cell2);

//     // Hook the faces up to the cells.
//     mesh.cells[cell1].push_back(mesh.faces.size()-1);
//     mesh.cells[cell2].push_back(mesh.faces.size()-1);

//     // Nodes that are attached to both cell 1 and cell2 are nodes 
//     // of this face.
//     set<int> cell1Nodes(cellNodes[cell1].begin(), cellNodes[cell1].end()),
//              cell2Nodes(cellNodes[cell2].begin(), cellNodes[cell2].end()),
//              faceNodes;
//     set_intersection(cell1Nodes.begin(), cell1Nodes.end(),
//                      cell2Nodes.begin(), cell2Nodes.end(),
//                      inserter(faceNodes, faceNodes.begin()));
//     mesh.faces.back().resize(faceNodes.size());
// //cout << "Face " << mesh.faces.size()-1 << " (between " << cell1 << " and " << cell2 << ") has " << faceNodes.size() << " nodes\n";
//     copy(faceNodes.begin(), faceNodes.end(), mesh.faces.back().begin());
//   }

//   // If no boundary was specified, leave.
//   if (geometry.empty()) 
//   {
//     // Clean up.
//     delete [] in.pointlist;
//     trifree((VOID*)delaunay.pointlist);
//     trifree((VOID*)delaunay.pointmarkerlist);
//     trifree((VOID*)delaunay.trianglelist);
//     trifree((VOID*)delaunay.edgelist);
//     trifree((VOID*)delaunay.edgemarkerlist);
//     if (geometry.empty())
//     {
//       trifree((VOID*)delaunay.segmentlist);
//       trifree((VOID*)delaunay.segmentmarkerlist);
//     }
//     else
//     {
//       delete [] in.segmentlist;
//       delete [] in.holelist;
//       delete [] delaunay.segmentlist;
//       delete [] delaunay.segmentmarkerlist;
//     }
//     return;
//   }

//   // At this point, all of the interior cells are set up properly, and all 
//   // nodes exist. However, we still have to finish constructing the 
//   // Voronoi cells that sit at the boundary, since they don't have "back" 
//   // faces. This Tessellator constructs faces that coincide with the 
//   // convex hull of the triangulation. 
  
//   set<pair<int, int> > patchedFaces;
//   for (int i = 0; i < delaunay.numberofsegments; ++i)
//   {
//     // A segment of the convex hull consists of the two cells/generators
//     // connected by a face.
//     int cell1 = delaunay.segmentlist[2*i];
//     int cell2 = delaunay.segmentlist[2*i+1];
// //cout << "Inspecting segment for " << cell1 << ", " << cell2 << endl;

//     // Add nodes to all faces attached to these cell with only 
//     // one node.
//     for (size_t f = 0; f < mesh.cells[cell1].size(); ++f)
//     {
//       int face = mesh.cells[cell1][f];
// //cout << "face " << face << " has " << mesh.faces[face].size() << " nodes\n";
//       if (mesh.faces[face].size() == 1)
//       {
//         // We've found a face with a single node. Find out where the
//         // other node should go and create it.
//         for (size_t f1 = 0; f1 < mesh.cells[cell2].size(); ++f1)
//         {
//           int face1 = mesh.cells[cell2][f1];
//           if (face1 == face)
//           {
//             // The node should bisect the segment on the convex hull 
//             // between cell1 and cell2.
//             RealType nodex = 0.5*(points[2*cell1]+points[2*cell2]),
//                  nodey = 0.5*(points[2*cell1+1]+points[2*cell2+1]);
// //cout << "Adding node (" << nodex << ", " << nodey << ") for cells " << cell1 << ", " << cell2 << endl;
//             mesh.nodes.push_back(nodex);
//             mesh.nodes.push_back(nodey);
//             mesh.faces[face].push_back(mesh.nodes.size()/2 - 1);
//             break;
//           }
//         }
//       }
//     }
//   }

//   // Create extra nodes for the boundary faces. These nodes coincide with 
//   // the Voronoi generater points.
//   size_t oldNumNodes = mesh.nodes.size() / 2;
// //cout << "Mesh had " << oldNumNodes << ", now adding " << mesh.convexHull.facets.size() << endl;
//   mesh.nodes.resize(2*(oldNumNodes + mesh.convexHull.facets.size()));
//   for (size_t i = 0; i < mesh.convexHull.facets.size(); ++i)
//   {
//     unsigned pindex = mesh.convexHull.facets[i][0];
//     mesh.nodes[2*(oldNumNodes+i)]   = points[2*pindex];
//     mesh.nodes[2*(oldNumNodes+i)+1] = points[2*pindex+1];
// //cout << "Adding node (" << points[2*pindex] << ", " << points[2*pindex+1] << ")\n";
//   }

//   // Now we construct the remaining faces for the boundary cells.
//   set<int> addedBoundaryFaces;
//   for (int i = 0; i < mesh.convexHull.facets.size(); ++i)
//   {
//     int cell = mesh.convexHull.facets[i][0]; // The boundary cell.
//     if (addedBoundaryFaces.find(cell) != addedBoundaryFaces.end())
//       continue;
//     else
//       addedBoundaryFaces.insert(cell);

//     // Add two new faces for this cell.
//     mesh.faces.resize(mesh.faces.size()+2);
//     mesh.faceCells.resize(mesh.faceCells.size()+2);

//     // Each face has one node that it shares with an existing face and 
//     // one that is still to be hooked up to the new node coinciding 
//     // with the generator point (and still has value -1).
//     map<unsigned, int> numFacesForNode;
//     for (size_t f = 0; f < mesh.cells[cell].size(); ++f)
//     {
//       int face = mesh.cells[cell][f];
//       for (size_t n = 0; n < mesh.faces[face].size(); ++n)
//         ++numFacesForNode[mesh.faces[face][n]];
//     }
//     int node1 = -1, node2 = -1;
//     for (map<unsigned, int>::const_iterator iter = numFacesForNode.begin();
//          iter != numFacesForNode.end(); ++iter)
//     {
//       if ((iter->second < 2) and ((node1 == -1) or (node2 == -1)))
//       {
//         if (node1 == -1)
//         {
//           node1 = iter->first;
//         }
//         else
//         {
//           ASSERT(node2 == -1);
//           node2 = iter->first;
//         }
//       }
//     }
//     ASSERT(node1 >= 0);
//     ASSERT(node2 >= 0);
//     mesh.faces[mesh.faces.size()-2].push_back(node1);
//     mesh.faces[mesh.faces.size()-1].push_back(node2);

//     // The remaining node on each face is the one coinciding with 
//     // the generator point.
//     int genPointNode = oldNumNodes+i;
//     mesh.faces[mesh.faces.size()-2].push_back(genPointNode);
//     mesh.faces[mesh.faces.size()-1].push_back(genPointNode);

//     // Each face has only the 1 cell.
//     mesh.faceCells[mesh.faceCells.size()-2].push_back(cell);
//     mesh.faceCells[mesh.faceCells.size()-1].push_back(cell);
//     mesh.cells[cell].push_back(mesh.faces.size()-2);
//     mesh.cells[cell].push_back(mesh.faces.size()-1);
//   }

  // Clean up.
  delete [] in.pointlist;
  trifree((VOID*)delaunay.pointlist);
  trifree((VOID*)delaunay.pointmarkerlist);
  trifree((VOID*)delaunay.trianglelist);
  trifree((VOID*)delaunay.edgelist);
  trifree((VOID*)delaunay.edgemarkerlist);
  if (geometry.empty())
  {
    trifree((VOID*)delaunay.segmentlist);
    trifree((VOID*)delaunay.segmentmarkerlist);
  }
  else
  {
    delete [] in.segmentlist;
    delete [] in.holelist;
    delete [] delaunay.segmentlist;
    delete [] delaunay.segmentmarkerlist;
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TriangleTessellator<double>;

}
