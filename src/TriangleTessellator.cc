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

// We use the Boost.Geometry library to handle polygon intersections and such.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/multi/geometries/multi_point.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>


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


//------------------------------------------------------------------------------
// Teach Boost.Geometry how to handle our Point2 class with appropriate traits.
//------------------------------------------------------------------------------
typedef int64_t CoordHash;
BOOST_GEOMETRY_REGISTER_POINT_2D(polytope::Point2<CoordHash>, CoordHash, boost::geometry::cs::cartesian, x, y)
BOOST_GEOMETRY_REGISTER_POINT_2D(polytope::Point2<double>, double, boost::geometry::cs::cartesian, x, y)


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
  // Use the PLC method with an empty geometry.
  PLC<2, RealType> geometry;
  tessellate(points, points, geometry, mesh);
  
  /*
  // Use an empty PLCpoints vector to compute Delaunay
  std::vector<RealType> PLCpoints;
  triangulateio delaunay;
  computeDelaunay(points, PLCpoints, delaunay);
  */
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

  typedef std::pair<int, int> EdgeHash;
  typedef Point2<CoordHash> IntPoint;
  typedef Point2<double> RealPoint;
  typedef boost::geometry::model::point<RealType, 2, boost::geometry::cs::cartesian> 
    realBGpoint;
  typedef boost::geometry::model::polygon<IntPoint,    // point type
                                          false>       // clockwise
    BGpolygon;
  typedef boost::geometry::model::ring<IntPoint,       // point type
                                       false>          // clockwise
    BGring;
  typedef boost::geometry::model::polygon<realBGpoint, // point type
                                          false>       // clockwise
    realBGpolygon;
  
  
  typedef boost::geometry::model::multi_point<IntPoint> BGmulti_point;
  typedef boost::geometry::model::multi_polygon<BGpolygon> BGmulti_polygon;

  const CoordHash coordMax = (1LL << 32); // numeric_limits<CoordHash>::max() >> 32U;
  const double degeneracy = 1.0e-12;
  
  triangulateio in, delaunay;

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

  // Define input points, including our false external generators.
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
  // cerr << "Chose bounding range : (" << xmin << " " << ymin << ") (" << xmax << " " << ymax << ")" << endl;

  // No point attributes or markers.
  in.numberofpointattributes = 0;
  in.pointattributelist = 0; 
  in.pointmarkerlist = 0; // No point markers.

  // Fill in Triangle's boundary info.  We use our imposed outer box of fake
  // generators.
  in.numberofsegments = 4;
  in.segmentlist = new int[2*in.numberofsegments];
  j = 0;
  for (i = 0; i != 4; ++i) {
    in.segmentlist[j++] = numGenerators + i;
    in.segmentlist[j++] = numGenerators + ((i + 1) % 4);
  }
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
  // if (geometry.empty())
  //   triangulate((char*)"Qzec", &in, &delaunay, 0);
  // else
  //   triangulate((char*)"Qzep", &in, &delaunay, 0);
  triangulate((char*)"Qzep", &in, &delaunay, 0);

  // Make sure we got something.
  if (delaunay.numberoftriangles == 0)
    error("TriangleTessellator: Delauney triangulation produced 0 triangles!");
  if (delaunay.numberofpoints != numGenerators + 4) {
    char err[1024];
    snprintf(err, 1024, "TriangleTessellator: Delauney triangulation produced %d triangles\n(%d generating points given)", 
             delaunay.numberofpoints, (int)numGenerators);
    error(err);
  }

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
  RealType  clow[2] = { numeric_limits<RealType>::max(),  numeric_limits<RealType>::max()};
  RealType chigh[2] = {-numeric_limits<RealType>::max(), -numeric_limits<RealType>::max()};
  internal::CounterMap<EdgeHash> edgeCounter;
  vector<RealPoint> circumcenters(delaunay.numberoftriangles);
  map<int, set<int> > gen2tri;
  int k, pindex, qindex, rindex, iedge;
  for (i = 0; i != delaunay.numberoftriangles; ++i) {
    pindex = delaunay.trianglelist[3*i];
    qindex = delaunay.trianglelist[3*i + 1];
    rindex = delaunay.trianglelist[3*i + 2];
    if (pindex < numGenerators and qindex < numGenerators and rindex < numGenerators) {
      ++edgeCounter[internal::hashEdge(pindex, qindex)];
      ++edgeCounter[internal::hashEdge(qindex, rindex)];
      ++edgeCounter[internal::hashEdge(rindex, pindex)];
    }
    geometry::computeCircumcenter2d(&delaunay.pointlist[2*pindex],
                                    &delaunay.pointlist[2*qindex],
                                    &delaunay.pointlist[2*rindex],
                                    &circumcenters[i].x);
    gen2tri[pindex].insert(i);
    gen2tri[qindex].insert(i);
    gen2tri[rindex].insert(i);
    // cerr << "circumcenter : " << circumcenters[i] << endl;
    // cerr << delaunay.pointlist[2*pindex]   << " " 
    //      << delaunay.pointlist[2*pindex+1] << " " 
    //      << delaunay.pointlist[2*qindex]   << " " 
    //      << delaunay.pointlist[2*qindex+1] << " " 
    //      << delaunay.pointlist[2*rindex]   << " " 
    //      << delaunay.pointlist[2*rindex+1] << " "
    //      << circumcenters[i].x             << " " 
    //      << circumcenters[i].y             << endl; 
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

  if( recursionDepth == 1 ){
     mLow.resize(2);     mHigh.resize(2);
     mLow[0] = clow[0];  mHigh[0] = chigh[0];
     mLow[1] = clow[1];  mHigh[1] = chigh[1];
     mdx = cdx;
  }
  
  // Flag any generators on the edge of the tessellation.  Here we mean the actual
  // generators, not our added boundary ones.
  //vector<bool> exteriorGenerators(numGenerators, false);
  //BoolMap<EdgeHash> exteriorEdgeTest;
  list<EdgeHash> exteriorEdges;
  //vector<vector<EdgeHash> > exteriorEdgesOfGen(numGenerators);
  for (typename internal::CounterMap<EdgeHash>::const_iterator itr = edgeCounter.begin();
       itr != edgeCounter.end();
       ++itr) {
    POLY_ASSERT(itr->second == 1 or itr->second == 2);
    if (itr->second == 1) {
      //i = itr->first.first;
      //j = itr->first.second;
      //exteriorGenerators[i] = true;
      //exteriorGenerators[j] = true;
      exteriorEdges.push_back(itr->first);
      //exteriorEdgeTest.insert(make_pair(itr->first, true));
      //exteriorEdgesOfGen[i].push_back(itr->first);
      //exteriorEdgesOfGen[j].push_back(itr->first);
    }
  }



  // Build the polygon representing our boundaries.
  BGpolygon boundary;
  realBGpolygon realBoundary;
  if (geometry.empty()) {
    // The user did not provide a boundary, so use our local edge use count to 
    // find the bounding edges.  Note in this case there will be no holes.
    const unsigned numBoundaryPoints = exteriorEdges.size();
    vector<IntPoint> boundaryPoints;
    vector<EdgeHash> boundaryEdgeOrder;
    boundaryEdgeOrder.push_back(exteriorEdges.front());
    exteriorEdges.pop_front();
    
    i = boundaryEdgeOrder.front().first;
    boundaryPoints.reserve(numBoundaryPoints + 1);
    boundaryPoints.push_back(IntPoint(points[2*i], points[2*i+1],
                                      mLow[0], mLow[1], mdx));
    boost::geometry::append( realBoundary, realBGpoint(points[2*i], points[2*i+1]));
    
    i = boundaryEdgeOrder.front().second;
    boundaryPoints.push_back(IntPoint(points[2*i], points[2*i+1],
                                      mLow[0], mLow[1], mdx));
    boost::geometry::append( realBoundary, realBGpoint(points[2*i], points[2*i+1]));
    while (exteriorEdges.size() > 0) {
      list<EdgeHash>::iterator itr = find_if(exteriorEdges.begin(), exteriorEdges.end(),
                                             MatchEitherPairValue<int>(i));
      POLY_ASSERT(itr != exteriorEdges.end());
      i = (itr->first == i ? itr->second : itr->first);
      boundaryEdgeOrder.push_back(*itr);
      exteriorEdges.erase(itr);
      boundaryPoints.push_back(IntPoint(points[2*i], points[2*i+1],
                                        mLow[0], mLow[1], mdx));
      boost::geometry::append( realBoundary, realBGpoint(points[2*i], points[2*i+1]));
    }
    POLY_ASSERT(boundaryEdgeOrder.size() == numBoundaryPoints);
    POLY_ASSERT(MatchEitherPairValue<int>(boundaryEdgeOrder.front().first)(boundaryEdgeOrder.back()));
    POLY_ASSERT(boundaryPoints.size() == numBoundaryPoints + 1);
    POLY_ASSERT(boundaryPoints.front() == boundaryPoints.back());
    boost::geometry::assign(boundary, BGring(boundaryPoints.begin(), boundaryPoints.end()));

  } else {
    // Copy the PLC provided boundary information into a Boost.Geometry polygon.
    vector<IntPoint> boundaryPoints;
    boundaryPoints.reserve(geometry.facets.size() + 1);
    i = geometry.facets[0][0];
    boundaryPoints.push_back(IntPoint(PLCpoints[2*i], PLCpoints[2*i+1],
                                      mLow[0], mLow[1], mdx));
    boost::geometry::append( realBoundary, realBGpoint(PLCpoints[2*i], PLCpoints[2*i+1]));

    for (j = 0; j != geometry.facets.size(); ++j) {
      POLY_ASSERT(geometry.facets[j].size() == 2);
      i =  geometry.facets[j][1];
      boundaryPoints.push_back(IntPoint(PLCpoints[2*i], PLCpoints[2*i+1],
                                        mLow[0], mLow[1], mdx));
      boost::geometry::append( realBoundary, realBGpoint(PLCpoints[2*i], PLCpoints[2*i+1]));
    }
    POLY_ASSERT(boundaryPoints.size() == geometry.facets.size() + 1);
    POLY_ASSERT(boundaryPoints.front() == boundaryPoints.back());
    boost::geometry::assign(boundary, BGring(boundaryPoints.begin(), boundaryPoints.end()));

    // Add any interior holes.
    const unsigned numHoles = geometry.holes.size();
    if (numHoles > 0) {
      typename BGpolygon::inner_container_type& holes = boundary.inners();
      holes.resize(numHoles);
      typename realBGpolygon::inner_container_type& realHoles = realBoundary.inners();
      realHoles.resize(numHoles);
      for (k = 0; k != numHoles; ++k) {
        boundaryPoints = vector<IntPoint>();
        boundaryPoints.reserve(geometry.holes[k].size() + 1);
        i = geometry.holes[k][0][0];
        boundaryPoints.push_back(IntPoint(PLCpoints[2*i], PLCpoints[2*i+1],
                                          mLow[0], mLow[1], mdx));
        boost::geometry::append( realHoles[k], realBGpoint(PLCpoints[2*i], PLCpoints[2*i+1]));
        for (j = 0; j != geometry.holes[k].size(); ++j) {
          POLY_ASSERT(geometry.holes[k][j].size() == 2);
          i =  geometry.holes[k][j][1];
          boundaryPoints.push_back(IntPoint(PLCpoints[2*i], PLCpoints[2*i+1],
                                            mLow[0], mLow[1], mdx));
          boost::geometry::append( realHoles[k], realBGpoint(PLCpoints[2*i], PLCpoints[2*i+1]));
        }
        POLY_ASSERT(boundaryPoints.size() == geometry.holes[k].size() + 1);
        POLY_ASSERT(boundaryPoints.front() == boundaryPoints.back());
        boost::geometry::assign(holes[k], BGring(boundaryPoints.begin(), boundaryPoints.end()));
      }
    }
  }

  // Walk each generator and build up it's unique nodes and faces.
  mesh.cells.resize(numGenerators);
  //CoordHash minR, thpt;
  IntPoint X;
  bool inside;
  map<IntPoint, int> point2node;
  map<IntPoint, set<int> > point2neighbors;
  map<EdgeHash, int> edgeHash2id;
  map<int, vector<int> > edgeCells;
  map<int, BGring> cellRings;
  map<int, vector<BGring> > orphanage;
  for (i = 0; i != numGenerators; ++i) {

    // Add the circumcenters as points for the cell.
    set<IntPoint> cellPointSet;
    for (set<int>::const_iterator triItr = gen2tri[i].begin();
         triItr != gen2tri[i].end();
         ++triItr) {
      cellPointSet.insert(IntPoint(circumcenters[*triItr].x, circumcenters[*triItr].y,
                                   mLow[0], mLow[1], mdx));
      // cerr << "Cell " << i << " adding circumcenter (" 
      //      << circumcenters[*triItr].x << " " << circumcenters[*triItr].y << ") "
      //      << IntPoint(circumcenters[*triItr].x, circumcenters[*triItr].y,
      //                  mLow[0], mLow[1], mdx) << endl;
    }
    POLY_ASSERT2(cellPointSet.size() >= 3, cellPointSet.size());

    // // Build the convex hull of the cell points.
    // vector<double> cellPointCoords;
    // for (j = 0; j != cellPoints.size(); ++j) {
    //   cellPointCoords.push_back(cellPoints[j].x);
    //   cellPointCoords.push_back(cellPoints[j].y);
    // }
    // POLY_ASSERT(cellPointCoords.size() == 2*cellPoints.size());
    // PLC<2, double> hull = convexHull_2d<double>(cellPointCoords, low, dx);
    // POLY_ASSERT(hull.facets.size() >= 3);
    // POLY_ASSERT(hull.facets[0][0] < cellPoints.size());
    // vector<RealPoint> ringPoints;
    // ringPoints.push_back(cellPoints[hull.facets[0][0]]);
    // for (j = 0; j != hull.facets.size(); ++j) {
    //   POLY_ASSERT(hull.facets[j].size() == 2);
    //   POLY_ASSERT(hull.facets[j][1] < cellPoints.size());
    //   ringPoints.push_back(cellPoints[hull.facets[j][1]]);
    // }
    // cellRings[i] = BGring(ringPoints.begin(), ringPoints.end());

    // Build the convex hull of the cell points.
    BGmulti_point mpoints(cellPointSet.begin(), cellPointSet.end());
    // for (typename set<IntPoint>::const_iterator itr = cellPointSet.begin();
    //      itr != cellPointSet.end();
    //      ++itr) mpoints.push_back(RealPoint(itr->realx(mLow[0], mdx),
    //                                         itr->realy(mLow[1], mdx)));
    boost::geometry::convex_hull(mpoints, cellRings[i]);

    // Intersect with the boundary to get the bounded cell.
    // Since for complex boundaries this may return more than one polygon, we find
    // the one that contains the generator.
    vector<BGring> cellIntersections;
    boost::geometry::intersection(boundary, cellRings[i], cellIntersections);
    if (cellIntersections.size() == 0) {
      cerr << points[2*i] << " " << points[2*i+1] << endl 
           << boost::geometry::dsv(cellRings[i]) << endl
           << boost::geometry::dsv(mpoints) << endl
           << boost::geometry::dsv(boundary) << endl;
    }
    POLY_ASSERT(cellIntersections.size() > 0);
    if (cellIntersections.size() == 1) {
      cellRings[i] = cellIntersections[0];
    } else {
      X = IntPoint(points[2*i], points[2*i+1],
                   mLow[0], mLow[1], mdx);
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
    for (typename BGring::const_iterator itr = cellRings[i].begin();
         itr != cellRings[i].end() - 1;
         ++itr) {
      // const IntPoint& pX1 = *itr;
      // const IntPoint& pX2 = *(itr + 1);
      // POLY_ASSERT(*itr != *(itr + 1));
      // j = internal::addKeyToMap(pX1, point2node);
      // k = internal::addKeyToMap(pX2, point2node);
      // POLY_ASSERT(j != k);
      // iedge = internal::addKeyToMap(internal::hashEdge(j, k), edgeHash2id);
      // edgeCells[iedge].push_back(j < k ? i : ~i);
      // mesh.cells[i].push_back(j < k ? iedge : ~iedge);

      std::map<IntPoint, set<int> >::iterator p2nItr = point2neighbors.find(*itr);
      point2neighbors[*itr].insert(i);
    }
    // POLY_ASSERT(mesh.cells[i].size() >= 3);
  }
  // POLY_ASSERT(edgeCells.size() == edgeHash2id.size());

  // Build the map from cell to set of neighboring cells  
  map<int, set<int> > neighbors;
  for (i = 0; i != numGenerators; ++i){
    for (typename BGring::const_iterator itr = cellRings[i].begin();
         itr != cellRings[i].end() - 1; ++itr) {
      neighbors[i].insert( point2neighbors[*itr].begin(), point2neighbors[*itr].end() );
    }
    neighbors[i].erase(i);
  }
  
  
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
  if( orphanage.size() > 0 && recursionDepth == 1 ){
    cerr << orphanage.size() << " orphaned cells" << endl;
    for (std::map<int, std::vector<BGring> >::const_iterator orphanItr = orphanage.begin();
         orphanItr != orphanage.end(); ++orphanItr){
      int parent = orphanItr->first;
      cerr << "Parent cell: " << parent << ":" << endl;
      for (unsigned iorphan = 0; iorphan != orphanItr->second.size(); ++iorphan){
        BGring orphan = orphanItr->second[iorphan];

        // Build the neighboring cells of the orphaned chunk
        std::set<int> orphanNeighbors;
        std::set<int> orphanNeighborhood;
        for (typename BGring::const_iterator pointItr = orphan.begin();
             pointItr != orphan.end() - 1; ++pointItr) {
          std::map<IntPoint, std::set<int> >::iterator it = point2neighbors.find( *pointItr );
          if (it != point2neighbors.end()){
            std::set<int> neighborSet = it->second;
            for (std::set<int>::const_iterator setItr = neighborSet.begin();
                 setItr != neighborSet.end(); ++setItr){
               orphanNeighbors.insert(*setItr);
               orphanNeighborhood.insert(*setItr);
               orphanNeighborhood.insert(neighbors[*setItr].begin(), neighbors[*setItr].end());
            }
          }
        }
        POLY_ASSERT( orphanNeighbors.size() > 0 );
        POLY_ASSERT( orphanNeighborhood.size() > 0 );
        
        
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
                        << "(" << (*itr).realx(clow[0],mdx) 
                        << "," << (*itr).realy(clow[1],mdx) << ")" << endl;
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
             subPLCpoints.push_back( (*itr).realx(clow[0],mdx) );
             subPLCpoints.push_back( (*itr).realy(clow[1],mdx) );
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
    
    // Check if nodes are inside boundary (either bounding box or PLC, if defined)
    bool inside = boost::geometry::within(realBGpoint(node[0],node[1]), realBoundary);
    if( !inside ){
       if (geometry.empty()) {
         node[0] = max(low[0], min(high[0], node[0]));
         node[1] = max(low[1], min(high[1], node[1]));
       } else {
         RealType result[2];
         RealType dist = nearestPoint( node, numPLCpoints, vertices, geometry, result );
         // Check the node has not moved more than 2.5 quantized mesh spacings. NOTE: this is
         // not a sharp estimate. Theoreticallly, the distance ought to be at most sqrt(2)*cdx, 
         // but nodes will fail this strict of a test.
         POLY_ASSERT2( dist < 2.5*mdx,
                       dist << " " << 2.5*mdx << " : (" << node[0] << " " << node[1] << ") (" << result[0] << " " << result[1] << ")\n" << geometry);
         node[0] = result[0];
         node[1] = result[1];
       }
    }
    
    POLY_ASSERT( node[0] >= low[0] && node[0] <= high[0] );
    POLY_ASSERT( node[1] >= low[1] && node[1] <= high[1] );
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
    trifree((VOID*) delaunay.segmentlist);
    trifree((VOID*) delaunay.segmentmarkerlist);
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TriangleTessellator<double>;





/*
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

  // Define input points, including our false external generators.
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
  // cerr << "Chose bounding range : (" << xmin << " " << ymin << ") (" << xmax << " " << ymax << ")" << endl;

  // No point attributes or markers.
  in.numberofpointattributes = 0;
  in.pointattributelist = 0; 
  in.pointmarkerlist = 0; // No point markers.

  // Fill in Triangle's boundary info.  We use our imposed outer box of fake
  // generators.
  in.numberofsegments = 4;
  in.segmentlist = new int[2*in.numberofsegments];
  j = 0;
  for (i = 0; i != 4; ++i) {
    in.segmentlist[j++] = numGenerators + i;
    in.segmentlist[j++] = numGenerators + ((i + 1) % 4);
  }
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
  // if (geometry.empty())
  //   triangulate((char*)"Qzec", &in, &delaunay, 0);
  // else
  //   triangulate((char*)"Qzep", &in, &delaunay, 0);
  triangulate((char*)"Qzep", &in, &delaunay, 0);

  // Make sure we got something.
  if (delaunay.numberoftriangles == 0)
    error("TriangleTessellator: Delauney triangulation produced 0 triangles!");
  if (delaunay.numberofpoints != numGenerators + 4) {
    char err[1024];
    snprintf(err, 1024, "TriangleTessellator: Delauney triangulation produced %d triangles\n(%d generating points given)", 
             delaunay.numberofpoints, (int)numGenerators);
    error(err);
  }



template<typename RealType>
void
TriangleTessellator<RealType>::
computeDelaunay(const vector<RealType>& points,
                const vector<RealType>& PLCpoints,
                triangulateio& delaunay) const 
{
  // Create the Voronoi nodes from the list of triangles. Each triangle 
  // has 3 nodes p, q, r, and corresponds to a Voronoi node at (X,Y), say.
  // The Voronoi node is located at the center of the triangle, though things
  // get a little squirrely at boundaries.  On boundary edges we create
  // a vertex at the edge center.

  // Find the circumcenters of each triangle, and build the set of triangles
  // associated with each generator.
  RealType  clow[2] = { numeric_limits<RealType>::max(),  numeric_limits<RealType>::max()};
  RealType chigh[2] = {-numeric_limits<RealType>::max(), -numeric_limits<RealType>::max()};
  internal::CounterMap<EdgeHash> edgeCounter;
  vector<RealPoint> circumcenters(delaunay.numberoftriangles);
  map<int, set<int> > gen2tri;
  int k, pindex, qindex, rindex, iedge;
  for (i = 0; i != delaunay.numberoftriangles; ++i) {
    pindex = delaunay.trianglelist[3*i];
    qindex = delaunay.trianglelist[3*i + 1];
    rindex = delaunay.trianglelist[3*i + 2];
    if (pindex < numGenerators and qindex < numGenerators and rindex < numGenerators) {
      ++edgeCounter[internal::hashEdge(pindex, qindex)];
      ++edgeCounter[internal::hashEdge(qindex, rindex)];
      ++edgeCounter[internal::hashEdge(rindex, pindex)];
    }
    geometry::computeCircumcenter2d(&delaunay.pointlist[2*pindex],
                                    &delaunay.pointlist[2*qindex],
                                    &delaunay.pointlist[2*rindex],
                                    &circumcenters[i].x);
    gen2tri[pindex].insert(i);
    gen2tri[qindex].insert(i);
    gen2tri[rindex].insert(i);
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

  if( recursionDepth == 1 ){
     mLow.resize(2);     mHigh.resize(2);
     mLow[0] = clow[0];  mHigh[0] = chigh[0];
     mLow[1] = clow[1];  mHigh[1] = chigh[1];
     mdx = cdx;
  }
  
  // Flag any generators on the edge of the tessellation.  Here we mean the actual
  list<EdgeHash> exteriorEdges;
  for (typename internal::CounterMap<EdgeHash>::const_iterator itr = edgeCounter.begin();
       itr != edgeCounter.end();
       ++itr) {
    POLY_ASSERT(itr->second == 1 or itr->second == 2);
    if (itr->second == 1) {
      exteriorEdges.push_back(itr->first);
    }
  }

*/

}
