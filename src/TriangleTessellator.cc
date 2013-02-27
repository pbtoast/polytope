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
int
findOtherTriIndex(const int* indices,
                  const int a,
                  const int b) {
  POLY_ASSERT(a == indices[0] or a == indices[1] or a == indices[2]);
  POLY_ASSERT(b == indices[0] or b == indices[1] or b == indices[2]);
  POLY_ASSERT(indices[0] != indices[1] and 
              indices[0] != indices[2] and 
              indices[1] != indices[2]);
  int c;
  if (a != indices[0] and b != indices[0]) {
    c = indices[0];
  }else {
    c = ((a == indices[1] or b == indices[1]) ? indices[2] : indices[1]);
  }
  return c;
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
// Find the circumcenters of each triangle, and build the set of triangles
// associated with each generator.
//------------------------------------------------------------------------------
template<typename RealType>
void
computeTriangleMaps(const RealType* pointList,
                    const RealType* triangleList,
                    const unsigned numberOfTriangles,
                    map<pair<int,int>, vector<unsigned> >& edge2tri,
                    map<int, set<unsigned> >& gen2tri,
                    vector<Point2<RealType> >& circumcenters) {
  int pindex, qindex, rindex, i, j;
  pair<int,int> pq, pr, qr;
  for (i = 0; i != numberOfTriangles; ++i) {
    pindex = triangleList[3*i  ];
    qindex = triangleList[3*i+1];
    rindex = triangleList[3*i+2];
    pq = internal::hashEdge(pindex, qindex);
    pr = internal::hashEdge(pindex, rindex);
    qr = internal::hashEdge(qindex, rindex);
    geometry::computeCircumcenter2d(pointList[2*pindex],
                                    pointList[2*qindex],
                                    pointList[2*rindex],
                                    &circumcenters[i].x);
    gen2tri[pindex].insert(i);
    gen2tri[qindex].insert(i);
    gen2tri[rindex].insert(i);
    edge2tri[pq].push_back(i);
    edge2tri[pr].push_back(i);
    edge2tri[qr].push_back(i);
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
// Compute bounding box which contains the "infinite" sphere for 
// unbounded tessellations
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
computeUnboundedBox(const RealType* points,
                    const unsigned numPoints,
                    const bool globalReduce,
                    RealType low[Dimension],
                    RealType high[Dimension]) {
  RealType xmin[Dimension], xmax[Dimension];
  geometry::computeBoundingBox(points, numPoints, globalReduce, xmin, xmax);
  const RealType rinf = -std::numeric_limits<RealType>::max();
  const RealType cent[Dimension];
  for(unsigned j=0; j != Dimension; ++j) {
    rinf = max( rinf, xmax[j]-xmin[j] );
    cent[j] = 0.5*( xmin[j]+xmax[j] );
  }
  for(unsigned j=0; j != Dimension; ++j) {
    low[j]  = cent[j] - rinf;
    high[j] = cent[j] + rinf;
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
  const unsigned numGenerators = points.size()/2;
  POLY_ASSERT(numGenerators > 1 );
  
  // Reset the bounding box and quantized mesh spacing
  mLow.assign(  2,  numeric_limits<RealType>::max() );
  mHigh.assign( 2, -numeric_limits<RealType>::max() );
  mdx = 0;
  
  if( numGenerators == 2 ){
    this->computeVoronoiFromCollinearPoints(points,mesh);
  }else{
    // Check that the points are not all collinear. If they are, the mesh is 1D
    // degenerate, and we cannot compute the Delaunay.
    bool collinear = true;
    int i = 2;
    while (collinear and i != numGenerators){
      collinear *= geometry::collinear<2,RealType>(&points[0], &points[2], &points[2*i], 1.0e-10);
      ++i;
    }

    if (collinear){
      this->computeVoronoiFromCollinearPoints(points,mesh);
    }else{
      this->computeVoronoi(points,mesh);
    }
  }
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
  // Reset the bounding box and quantized mesh spacing
  mLow.assign(  2,  numeric_limits<RealType>::max() );
  mHigh.assign( 2, -numeric_limits<RealType>::max() );
  mdx = 0;
  
  POLY_ASSERT(!points.empty());
  
  // Make sure we're not modifying an existing tessellation.
  POLY_ASSERT(mesh.empty());

  // Find the range of the generator points.
  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints = PLCpoints.size()/2;
  int i, j, k;

  // compute bounded cell rings from an unbounded tessellation
  vector<BGring> cellRings;
  map<int, vector<BGring> > orphanage;
  this->computeCellRings(points, PLCpoints, geometry, 
                         cellRings, orphanage);
  POLY_ASSERT( cellRings.size() == numGenerators );
  
  //*********************** Begin Adoption Algorithm ************************
  // Intersecting cells with the boundary has created orphaned cell pieces. ("Won't 
  // somebody please think of the children!?") Make sub-tessellations using the 
  // generators neighboring the orphaned pieces. Construct a PLC boundary for the 
  // sub-tessellation by using the geometry obtained by union-ing the orphan 
  // with its neighboring cells. The way we compute cell neighbors should ensure
  // that the union gives a contiguous geometry with no holes
  if( orphanage.size() > 0 ){
    
    // Construct map from node points to neighboring cells
    map<IntPoint, set<int> > point2neighbors;
    for (i = 0; i != cellRings.size(); ++i){
      for (typename BGring::const_iterator itr = cellRings[i].begin();
	   itr != cellRings[i].end() - 1; ++itr){
	point2neighbors[*itr].insert(i);
      }
    }
    
    // // First aglomerate all orphaned pieces that neighbor one another by searching orphan pairs
    // // TODO: Figure out a more efficient way to do this operation
    // for (i = 0; i < orphanage.size()-1; ++i){
    //   for (j = i+1; j < orphanage.size(); ++j){
    //     BGring unionRing, orphan1 = orphanage[i], orphan2 = orphanage[j];
    //     for (typename BGring::const_iterator pointItr = orphan1.begin();
    //          pointItr != orphan1.end()-1; ++pointItr) {
    //       typename BGring::const_iterator it = orphan2.find(*pointItr);
    //       if (it != orphan2.end()){
    //         boost::geometry::union_( orphan1, orphan2, unionRing );
    //         // TODO: store the unionRing into the resized orphanage
    //       }
    //     }
    //   }
    // }

    for (map<int,vector<BGring> >::const_iterator itr = orphanage.begin();
         itr != orphanage.end(); ++itr){
      int parent = itr->first;
      // cerr << "***Cell " << parent << " has " << itr->second.size() << " orphans:" << endl;
      for (i = 0; i != itr->second.size(); ++i){	
	BGring orphan = itr->second[i];
	std::set<int> orphanNeighbors;
	for (typename BGring::const_iterator pointItr = orphan.begin();
	     pointItr != orphan.end()-1; ++pointItr) {
	  std::map<IntPoint, std::set<int> >::iterator it = point2neighbors.find(*pointItr);
	  if (it != point2neighbors.end()){
	    std::set<int> neighborSet = it->second;
	    for (std::set<int>::const_iterator setItr = neighborSet.begin();
		 setItr != neighborSet.end(); ++setItr){
	      orphanNeighbors.insert(*setItr);
	    }
	    orphanNeighbors.erase(parent);
	  }
	}
	POLY_ASSERT( orphanNeighbors.size() > 0 );
	
	
	// // Blago!
	// cerr << "Orphaned piece has neighbors";
	// for( std::set<int>::const_iterator iii = orphanNeighbors.begin();
	//      iii != orphanNeighbors.end(); ++iii){
	//   cerr << " " << *iii;
	// }
	// cerr << endl;
	// // Blago!

	
	// If the orphan only has a single neighbor, we can skip a lot of work.
	// No need to tessellate - simply union the orphan with its neighbor cell.
	vector<BGring> subCellRings;
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
	  for (unsigned ii = 0; ii < nSides; ++ii) {
	    subPLC.facets[ii][0] = ii;
	    subPLC.facets[ii][1] = (ii+1) % nSides;
	  }
	  
	  map<int, vector<BGring> > subOrphanage;
	  this->computeCellRings(subpoints, subPLCpoints, subPLC,
                                 subCellRings, subOrphanage);
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
	  
	  if (orphanNeighbors.size() > 1){
	    thisRing = subCellRings[subIndex];
	    
	    // // Blago!
	    // cerr << endl << "Cell " << thisIndex << endl;
	    // cerr << endl << "SUBMESH CELL:" << endl;
	    // for (typename BGring::const_iterator itr = thisRing.begin();
	    //      itr != thisRing.end(); ++itr){
	    //   cerr << (*itr).realx(mLow[0],mdx) << " " 
	    //        << (*itr).realy(mLow[1],mdx) << endl;
	    // }
	    // // Blago!
	    
	    
	    // Simplify the resulting ring. Removes points that are within some minimum
	    // distance to their neighbors. Setting distance = 1 merges ring elements
	    // that are within one quantized mesh spacing. This essentially removes
	    // repeated cell nodes having length-zero cell faces. An unfortunate
	    // consequence of using a third-party lib to do all our unions/intersections
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
          if(unionRing.size() > 1){
             cerr << "Blago!" << endl << "Cell " << thisIndex
                  << " has more than one cell ring:" << endl;
             for( i=0; i<unionRing.size(); ++i){
                cerr << endl << "Ring " << i << ":" << endl;
                for (typename BGring::const_iterator itr = thisRing.begin();
                     itr != thisRing.end(); ++itr){
                   cerr << (*itr).realx(mLow[0],mdx) << " " 
                        << (*itr).realy(mLow[1],mdx) << " "
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
	  
	  
	  // // Blago!
	  // cerr << endl << "Cell " << thisIndex << endl;
	  // cerr << endl << "FINAL SUBMESH CELL:" << endl;
	  // for (typename BGring::const_iterator itr = thisRing.begin();
	  //      itr != thisRing.end(); ++itr){
	  //   cerr << (*itr).realx(mLow[0],mdx) << " " 
	  //        << (*itr).realy(mLow[1],mdx) << endl;
	  // }
	  // for (typename BGring::const_iterator itr = thisRing.begin();
	  //      itr != thisRing.end(); ++itr){
	  //   cerr << *itr << endl;
	  // }
	  // // Blago!
	  
	  
	  cellRings[thisIndex] = thisRing;
	}
      }
    }
  }
  //*********************** End Adoption Algorithm ************************




  // Now build the unique mesh nodes and cell info.
  IntPoint X, IntNode;
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
    node[0] = p.realx(mLow[0],mdx);
    node[1] = p.realy(mLow[1],mdx);
    
    // Check if nodes are inside boundary (either bounding box or PLC, if defined)
    //bool inside = boost::geometry::within(realBGpoint(node[0],node[1]), realBoundary);
    bool inside = within(node, numPLCpoints, &PLCpoints[0], geometry);
    
    if( !inside ){
      RealType result[2];
      RealType dist = nearestPoint( node, numPLCpoints, &PLCpoints[0], geometry, result );
      // Check the node has not moved more than 2.5 quantized mesh spacings. NOTE: this is
      // not a sharp estimate. Theoreticallly, the distance ought to be at most sqrt(2)*cdx, 
      // but nodes will fail this strict of a test.
      POLY_ASSERT2( dist < 2.5*mdx,
                    dist << " " << 2.5*mdx << " : "
                    << "(" << node[0]   << " " << node[1]   << "( "
                    << "(" << result[0] << " " << result[1] << ")\n" << geometry );
      node[0] = result[0];
      node[1] = result[1];
    }
     
    // POLY_ASSERT( node[0] >= low[0] and node[0] <= high[0] );
    // POLY_ASSERT( node[1] >= low[1] and node[1] <= high[1] );
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
//PRIVATE STUFF:
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeVoronoiFromCollinearPoints(const vector<RealType>& points,
				  Tessellation<2, RealType>& mesh) const {
  POLY_ASSERT(mesh.empty());
  
  const unsigned numGenerators = points.size()/2;
  
  const CoordHash coordMax = (1LL << 31); // numeric_limits<CoordHash>::max() >> 32U;
  const double degeneracy = 4.0e-10;
  int i;
  
  vector<pair<RealPoint,int> > pointIndexPairs;
  for (i = 0; i != numGenerators; ++i){
    pointIndexPairs.push_back(make_pair(RealPoint(points[2*i], points[2*i+1]), i));
  }
  sort( pointIndexPairs.begin(), 
	pointIndexPairs.end(), 
	internal::pairCompareFirst<RealPoint,int> );
  
  // Bounding box for the points
  RealType low[2], high[2];
  geometry::computeBoundingBox<2,RealType>(points, true, low, high);
  POLY_ASSERT(low[0] <= high[0] and low[1] <= high[1]);

  // The bounding box which contains PLC, and all circumcenters and generators
  RealType cbox[2] = {high[0] - low[0], high[1] - low[1]};
  
  // The bounding circle onto which we project the "infinite" rays of the 
  // unbounded faces of the tessellation.
  const RealType rtmp    = 2.0*max(cbox[0], cbox[1]);
  const RealType ctmp[2] = {0.5*(low[0]+high[0]), 0.5*(low[1]+high[1])};

  // We resize mLow and boxsize so that the bounding box
  // contains the "infinite" sphere. mHigh is not really needed.
  mLow [0] = min(mLow [0], ctmp[0]-rtmp);
  mLow [1] = min(mLow [1], ctmp[1]-rtmp);
  mHigh[0] = max(mHigh[0], ctmp[0]+rtmp);  
  mHigh[1] = max(mHigh[1], ctmp[1]+rtmp);
  
  const RealType rinf     = 0.5*(mHigh[0]-mLow[0]);
  const RealType cboxc[2] = {0.5*(mLow[0]+mHigh[0]), 0.5*(mLow[1]+mHigh[1])};

  const double cboxsize = 2.0*rinf;
  mdx = max(degeneracy, cboxsize/coordMax);

  // Sizes for the node, face, cell lists
  mesh.cells.resize(numGenerators);
  mesh.nodes.resize(4*numGenerators);
  mesh.faces.resize(numGenerators+3);
  mesh.faceCells.resize(numGenerators+3);
  
  bool test;
  unsigned inode, iface, icell1, icell2;
  RealPoint p1, p2, r1, r2, node, midpt;
  IntPoint IntNode;
  vector<int> cellFaces(2);
  for (i = 0; i != numGenerators-1; ++i){
    inode  = 2*i;
    iface  = i;
    icell1 = pointIndexPairs[i  ].second;
    icell2 = pointIndexPairs[i+1].second;

    p1    = pointIndexPairs[i  ].first;
    p2    = pointIndexPairs[i+1].first;
    midpt = RealPoint( 0.5*(p1.x + p2.x),
		       0.5*(p1.y + p2.y) );
    r1.x = p2.x - p1.x;
    r1.y = p2.y - p1.y;
    geometry::unitVector<2,RealType>(&r1.x);
    r2.x = -r1.y;
    r2.y =  r1.x;

    // Node 2*i: endpt of interior face
    test = geometry::rayCircleIntersection(&midpt.x, &r2.x, cboxc, rinf, 1.0e-10, &node.x);
    POLY_ASSERT(test);
    IntNode = IntPoint(node.x, node.y, mLow[0], mLow[1], mdx);
    mesh.nodes[2*inode  ] = IntNode.realx(mLow[0],mdx);
    mesh.nodes[2*inode+1] = IntNode.realy(mLow[1],mdx);
    mesh.infNodes.push_back(1);
    
    // Node 2*i+1: other endpt of interior face
    r2 *= -1.0;
    test = geometry::rayCircleIntersection(&midpt.x, &r2.x, cboxc, rinf, 1.0e-10, &node.x);
    POLY_ASSERT(test);
    IntNode = IntPoint(node.x,node.y,mLow[0],mLow[1],mdx);
    mesh.nodes[2*(inode+1)  ] = IntNode.realx(mLow[0],mdx);
    mesh.nodes[2*(inode+1)+1] = IntNode.realy(mLow[1],mdx);
    mesh.infNodes.push_back(1);
    
    // Mesh cells: the interior face between cells i and i+1. Construct mesh
    // so that this face is oriented positively for cell i
    mesh.cells[icell1].push_back( iface);
    mesh.cells[icell2].push_back(~iface);
    
    // Nodes around the interior mesh faces
    mesh.faces[iface].push_back(inode  );
    mesh.faces[iface].push_back(inode+1);
    mesh.infFaces.push_back(0);
    
    // Cells adjacent to each interior face. Positively oriented for cell i.
    mesh.faceCells[iface].push_back( icell1);
    mesh.faceCells[iface].push_back(~icell2);
  }

  // ------------------ Extra node at index 0 ----------------- //

  inode  = 2*(numGenerators-1);
  iface  = numGenerators-1;
  icell1 = pointIndexPairs[0].second;
  icell2 = pointIndexPairs[1].second;

  // Node position
  p1   = pointIndexPairs[icell1].first;
  p2   = pointIndexPairs[icell2].first;
  r1.x = p1.x - p2.x;
  r1.y = p1.y - p2.y;
  geometry::unitVector<2,RealType>(&r1.x);
  
  test = geometry::rayCircleIntersection(&p1.x, &r1.x, cboxc, rinf, 1.0e-10, &node.x);
  POLY_ASSERT(test);
  IntNode = IntPoint(node.x, node.y, mLow[0], mLow[1], mdx);
  mesh.nodes[2*inode  ] = IntNode.realx(mLow[0],mdx);
  mesh.nodes[2*inode+1] = IntNode.realy(mLow[1],mdx);
  mesh.infNodes.push_back(1);
  
  // Extra faces for cell 0
  mesh.cells[icell1].push_back(iface  );
  mesh.cells[icell1].push_back(iface+1);

  // Nodes for those two faces
  mesh.faces[iface].push_back(1    );
  mesh.faces[iface].push_back(inode);
  mesh.faces[iface+1].push_back(inode);
  mesh.faces[iface+1].push_back(0    );
  
  // Both faces are inf faces
  mesh.infFaces.push_back(1);
  mesh.infFaces.push_back(1);

  // Only cell 0 is around those faces
  mesh.faceCells[iface  ].push_back(icell1);
  mesh.faceCells[iface+1].push_back(icell1);


  // ------------------ Extra node at index N ----------------- //
  
  inode  = 2*numGenerators-1;
  iface  = numGenerators+1;
  icell1 = pointIndexPairs[numGenerators-1].second;
  icell2 = pointIndexPairs[numGenerators-2].second;

  // Node position
  p1   = pointIndexPairs[numGenerators-1].first;
  p2   = pointIndexPairs[numGenerators-2].first;
  r1.x = p1.x - p2.x;
  r1.y = p1.y - p2.y;
  geometry::unitVector<2,RealType>(&r1.x);
  
  test = geometry::rayCircleIntersection(&p2.x, &r1.x, cboxc, rinf, 1.0e-10, &node.x);
  POLY_ASSERT(test);
  IntNode = IntPoint(node.x, node.y, mLow[0], mLow[1], mdx);
  mesh.nodes[2*inode  ] = IntNode.realx(mLow[0],mdx);
  mesh.nodes[2*inode+1] = IntNode.realy(mLow[1],mdx);
  mesh.infNodes.push_back(1);
  
  // Extra faces for cell N
  mesh.cells[icell1].push_back(iface  );
  mesh.cells[icell1].push_back(iface+1);

  // Nodes for those two faces
  mesh.faces[iface  ].push_back(2*numGenerators-4 );
  mesh.faces[iface  ].push_back(inode             );
  mesh.faces[iface+1].push_back(inode             );
  mesh.faces[iface+1].push_back(2*numGenerators-3 );
  
  // Both faces are inf faces
  mesh.infFaces.push_back(1);
  mesh.infFaces.push_back(1);

  // Only cell 0 is around those faces
  mesh.faceCells[iface  ].push_back(icell1);
  mesh.faceCells[iface+1].push_back(icell1);

  // Blago!
  cerr << mesh << endl;
  // Blago!

}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
computeVoronoi(const vector<RealType>& points,
               Tessellation<2, RealType>& mesh) const {

  POLY_ASSERT(!points.empty());
  POLY_ASSERT(points.size() != 2);

  // Make sure we're not modifying an existing tessellation.
  POLY_ASSERT(mesh.empty());

  const CoordHash coordMax = (1LL << 31); // numeric_limits<CoordHash>::max() >> 32U;
  const double degeneracy = 4.0e-10;
  
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
  RealType low [2] = { numeric_limits<RealType>::max(), 
		       numeric_limits<RealType>::max()};
  RealType high[2] = {-numeric_limits<RealType>::max(), 
		      -numeric_limits<RealType>::max()};  
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
    low[0]  = min(low[0], circumcenters[i].x);
    low[1]  = min(low[1], circumcenters[i].y);
    high[0] = max(high[0], circumcenters[i].x);
    high[1] = max(high[1], circumcenters[i].y);
  }
  POLY_ASSERT(circumcenters.size() == delaunay.numberoftriangles);

  // The circumcenters may all lie inside the convex hull of the
  // generators for an unbounded tessellation. Include the generator
  // locations in the high/low search
  geometry::expandBoundingBox<2,RealType>(&delaunay.pointlist[0],
					  2*delaunay.numberofpoints,
					  true, low, high);
  POLY_ASSERT(low[0] <= high[0] and low[1] <= high[1]);

  // The bounding box which contains all circumcenters and generators
  RealType cbox[2] = {high[0] - low[0], high[1] - low[1]};
  
  // The bounding circle onto which we project the "infinite" rays of the 
  // unbounded faces of the tessellation.
  const RealType rtmp    = 2.0*max(cbox[0], cbox[1]);
  const RealType ctmp[2] = {0.5*(low[0]+high[0]), 0.5*(low[1]+high[1])};
  
  // We resize mLow and boxsize so that the bounding box
  // contains the "infinite" sphere. mHigh is not really needed.
  // cerr << "mLow  = (" << mLow[0]      << "," << mLow[1]      << ")" << endl;
  // cerr << "mHigh = (" << mHigh[0]     << "," << mHigh[1]     << ")" << endl;
  // cerr << "c-r   = (" << ctmp[0]-rtmp << "," << ctmp[1]-rtmp << ")" << endl;
  // cerr << "c+r   = (" << ctmp[0]+rtmp << "," << ctmp[1]+rtmp << ")" << endl;
  mLow [0] = min(mLow [0], ctmp[0]-rtmp);
  mLow [1] = min(mLow [1], ctmp[1]-rtmp);
  mHigh[0] = max(mHigh[0], ctmp[0]+rtmp);  
  mHigh[1] = max(mHigh[1], ctmp[1]+rtmp);
  
  const RealType rinf     = 0.5*(mHigh[0]-mLow[0]);
  const RealType cboxc[2] = {0.5*(mLow[0]+mHigh[0]), 0.5*(mLow[1]+mHigh[1])};

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
  
  // // Blago!
  // cerr << "Voronoi Bounding box:"
  //      << "(" << mLow[0] << "," << mHigh[0] << ")x"
  //      << "(" << mLow[1] << "," << mHigh[1] << ")" << endl;
  // cerr << "Box size    = " << cboxsize << endl;
  // cerr << "spacing     = " << mdx << endl;
  // cerr << "radius      = " << rinf << endl;
  // cerr << "center      = (" << cboxc[0] << "," << cboxc[1] << ")" << endl;
  // // Blago!

  // The exterior edges of the triangularization have "unbounded" rays, originating
  // at the circumcenter of the corresponding triangle and passing perpendicular to
  // the edge
  bool test;
  RealPoint ehat, pinf;
  map<EdgeHash, unsigned> edge2id;
  int i1, i2, ivert, k;
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
      ivert = findOtherTriIndex(&delaunay.trianglelist[3*i], i1, i2);
      
      computeEdgeUnitVector(&delaunay.pointlist[2*i1],
			    &delaunay.pointlist[2*i2],
			    &delaunay.pointlist[2*ivert],
			    &ehat.x);
       
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
  vector<unsigned> faceVec(2);
  internal::CounterMap<unsigned> faceCounter;
  mesh.cells = vector<vector<int> >(numGenerators);
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
    
    // Get the face sorted nodes
    vector<unsigned> faceNodes = 
       computeSortedFaceNodes(vector<EdgeHash>(meshEdges.begin(), meshEdges.end()));
    
    // Check the orientation of the nodes around the generator. Reverse the order
    // if they're ordered clockwise. It is enough to check the order of the first
    // two nodes in the faceNodes list, since its node list is sequential
    POLY_ASSERT(faceNodes.size() > 1);
    if( orient2d(&mesh.nodes[2*faceNodes[0]], 
		 &mesh.nodes[2*faceNodes[1]], 
		 &delaunay.pointlist[2*pindex]) < 0 ){
      reverse(faceNodes.begin(), faceNodes.end());
    }
    
    // The ordered mesh nodes around a given generator
    POLY_ASSERT(faceNodes.size() > 2);
    for (i = 0; i != faceNodes.size(); ++i){
      i1 = faceNodes[i];
      i2 = faceNodes[(i+1) % faceNodes.size()];
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
        if ( mesh.infNodes[i1] == 1 and mesh.infNodes[i2] == 1 ){
          mesh.infFaces.push_back(1);
        }else{
          mesh.infFaces.push_back(0);
        }
      }
      
      // Store the cell-face info based on face orientation
      if( orient2d(&mesh.nodes[2*face.first], 
                   &mesh.nodes[2*face.second], 
                   &delaunay.pointlist[2*pindex]) > 0 ){
        mesh.cells[pindex].push_back(iface);
        mesh.faceCells[iface].push_back(pindex);
      }else{
        mesh.cells[pindex].push_back(~iface);
        mesh.faceCells[iface].push_back(~int(pindex));
      }
    }
  }
  
  // // Blago!
  // cerr << mesh;
  // // Blago!
  
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
		 const vector<RealType>& PLCpoints,
		 const PLC<2, RealType>& geometry,
		 vector<BGring>& cellRings,
		 map<int, vector<BGring> >& orphanage) const 
{
  POLY_ASSERT(!points.empty());
  POLY_ASSERT2(!PLCpoints.empty(), "Error: attempting to create a bounded "
               << "tessellation with no bounding points");
  POLY_ASSERT2(!geometry.empty(),  "Error: attempting to create a bounded "
               << "tessellation with no bounding PLC");

  // Find the range of the generator points.
  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints = PLCpoints.size()/2;
  int i, j, k;
  
  RealType low[2], high[2];
  geometry::computeBoundingBox<2,RealType>(PLCpoints, true, low, high);
  POLY_ASSERT(low[0] < high[0] and low[1] < high[1]);  

  mLow[0] = min(mLow[0], low[0]);  mHigh[0] = max(mHigh[0], high[0]);
  mLow[1] = min(mLow[1], low[1]);  mHigh[1] = max(mHigh[1], high[1]);

  // Start by creating an unbounded tessellation
  Tessellation<2,RealType> mesh;
  if( numGenerators == 2 ){
     this->computeVoronoiFromCollinearPoints(points,mesh);
  }else{
     this->computeVoronoi(points,mesh);
  }
  
  // Bounding circle radius and center
  const RealType rinf = 0.5*(mHigh[0]-mLow[0]);
  const RealType cboxc[2] = {0.5*(mLow[0]+mHigh[0]), 0.5*(mLow[1]+mHigh[1])};
     
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
  IntPoint X, IntNode;
  RealPoint ehat, pinf;
  bool inside, test;
  cellRings.resize(numGenerators);
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
      cellBoundary.push_back(IntNode);
      // If both face nodes are infNodes, the inf face that connects them crosses the
      // interior of the "infinite" bounding circle. Test if the face intersects
      // the PLC boundary
      
      if (mesh.infFaces[iface] == 1) {
	vector<RealType> result;
	unsigned nints = intersect(&mesh.nodes[2*inode1], &mesh.nodes[2*inode2], 
				   numPLCpoints, &PLCpoints[0], geometry, result);
	// If it intersects, create a new inf node between them by projecting
	// the generator position out to the bounding circle along the ray
	// pointing normal to the original inf face.
	if (nints > 0){
          // // Blago!
          // cerr << "Segment connecting inf nodes has crossed the PLC boundary!" << endl;
          // cerr << "Cell " << i << endl;
          // cerr << "  inf node 1 : (" << mesh.nodes[2*inode1] << "," 
          //      << mesh.nodes[2*inode1+1] << ")" << endl;
          // cerr << "  inf node 2 : (" << mesh.nodes[2*inode2] << "," 
          //      << mesh.nodes[2*inode2+1] << ")" << endl;
          // cerr << "  " << nints << " intersections:" << endl;
          // for (j=0; j<result.size()/2; ++j){
          //    cerr << "    (" << result[2*j] << "," << result[2*j+1] << ")" << endl;
          // }
          // // Blago!
	  computeEdgeUnitVector(&mesh.nodes[2*inode1],
				&mesh.nodes[2*inode2],
				(RealType*)&points[2*i],
				&ehat.x);
	  test = geometry::rayCircleIntersection(&points[2*i],
						 &ehat.x,
						 cboxc,
						 rinf,
						 1.0e-10,
						 &pinf.x);
	  POLY_ASSERT(test);
          // cerr << "  new node   : (" << pinf.x << "," << pinf.y << ")" << endl;
	  
          // Add the fictitious inf node
	  IntNode = IntPoint(pinf.x, pinf.y, mLow[0], mLow[1], mdx);
	  cellBoundary.push_back(IntNode);
	}
      }
    }
    POLY_ASSERT(cellBoundary.size() > 0);
    cellBoundary.push_back( cellBoundary[0] );  // Close the ring
    boost::geometry::assign(cellRings[i], BGring(cellBoundary.begin(), cellBoundary.end()));
    POLY_ASSERT(cellRings[i].size() > 0);

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
