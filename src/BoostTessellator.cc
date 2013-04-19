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
#include "polytope_tessellator_utilities.hh"
#include "within.hh"
#include "nearestPoint.hh"
#include "intersect.hh"
#include "Clipper2d.hh"

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
// The vertex of an edge is infinite if its pointer is null (Boost
// criteria) OR its coordinate is out of the bounding box (our criteria)
//------------------------------------------------------------------------
inline
bool isFinite(const VD::vertex_type* v, const unsigned coordMax) {
  if(v) return ((v->x() >= 0) and (v->x() <= coordMax) and 
                (v->y() >= 0) and (v->y() <= coordMax));
  else return false;
}


//------------------------------------------------------------------------
// Comparison operator for voronoi vertices
//------------------------------------------------------------------------
struct vertexCompare {
   bool operator()(const VD::vertex_type v1, const VD::vertex_type v2) {
      return (v1.x() < v2.x()                      ? true :
              v1.x() == v2.x() and v1.y() < v2.y() ? true :
              false);
   }
};


//------------------------------------------------------------------------
// compute the unit vector direction of an infinite edge in the Boost
// Voronoi diagram
//------------------------------------------------------------------------
template<typename RealType>
void
computeInfiniteEdgeDirection(const VD::edge_type* edge,
                             const vector<IntPoint>& points,
			     const RealType* low,
			     const RealType delta,
                             RealType* direction) {
  POLY_ASSERT(edge->is_infinite());
  const VD::cell_type* cell1 = edge->cell();
  const VD::cell_type* cell2 = edge->twin()->cell();
  // Assume only point-generators for the time being
  POLY_ASSERT(cell1->contains_point() and cell2->contains_point());
  size_t index1 = cell1->source_index();
  size_t index2 = cell2->source_index();
  POLY_ASSERT(index1 < points.size() and index2 < points.size());
  const IntPoint p1 = points[index1];
  const IntPoint p2 = points[index2];
  RealType r[2] = {p2.realx(low[0],delta) - p1.realx(low[0],delta),
                   p2.realy(low[1],delta) - p1.realy(low[1],delta)};
  if (edge->vertex0()) {direction[0] = -r[1];  direction[1] =  r[0];}
  else                 {direction[0] =  r[1];  direction[1] = -r[0];}
  geometry::unitVector<2, RealType>(direction);
}


} //end anonymous namespace



//------------------------------------------------------------------------------
template<typename RealType>
BoostTessellator<RealType>::
BoostTessellator():
  Tessellator<2, RealType>() {
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
BoostTessellator<RealType>::
~BoostTessellator() {
}
//------------------------------------------------------------------------------

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

  const CoordHash coordMax = (1LL << 26);
  const double degeneracy  = 1.5e-8;

  // Compute the bounding box for this problem
  RealType low[2], high[2];
  geometry::computeBoundingBox<2,RealType>(points, true, low, high);
  POLY_ASSERT(low[0] <= high[0] and low[1] <= high[1]);

  // Infinite-radius circle for inf nodes
  const RealType box[2] = {high[0] - low[0],
			   high[1] - low[1]};
  const RealType rinf   = 4.0*max(box[0], box[1]);
  const RealType cen[2] = {0.5*(low[0] + high[0]),
			   0.5*(low[1] + high[1])};
  
  // Final bounding box dimensions for quantization
  mLow.resize(2);  mHigh.resize(2);
  mLow [0] = min(low [0], cen[0] - rinf);
  mHigh[0] = max(high[0], cen[0] + rinf);
  mLow [1] = min(low [1], cen[1] - rinf);
  mHigh[1] = max(high[1], cen[1] + rinf);
  mDelta = max(degeneracy, 2.0*rinf/coordMax);
  POLY_ASSERT(mLow[0] <= mHigh[0] and mLow[1] <= mHigh[1]);

  mRinf = rinf;
  mCenter.push_back(0.5*(mLow[0] + mHigh[0]));
  mCenter.push_back(0.5*(mLow[1] + mHigh[1]));
  mCoordMax = coordMax;

  this->computeVoronoiUnbounded(points, mesh);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
tessellate(const vector<RealType>& points,
           RealType* low,
           RealType* high,
           Tessellation<2, RealType>& mesh) const {
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(low != 0 and high != 0);
  POLY_ASSERT(points.size() > 0);
  POLY_ASSERT(points.size() % 2);
  
  // Build a PLC with the bounding box, and then use the PLC method.
  ReducedPLC<2, RealType> box = this->boundingBox(low, high);
  this->tessellate(points, box.points, box, mesh);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const vector<RealType>& PLCpoints,
           const PLC<2, RealType>& geometry,
           Tessellation<2, RealType>& mesh) const {
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(!geometry.empty());
  POLY_ASSERT(points.size() > 0 and PLCpoints.size() > 0);
  POLY_ASSERT(points.size() % 2 == 0 and PLCpoints.size() % 2 == 0);

  const CoordHash coordMax = (1LL << 26);
  const double degeneracy  = 1.5e-8;

  // Compute the bounding box for this problem
  RealType low[2], high[2];
  geometry::computeBoundingBox<2,RealType>(PLCpoints, true, low, high);
  POLY_ASSERT(low[0] <= high[0] and low[1] <= high[1]);

  // Infinite-radius circle for inf nodes
  const RealType box[2] = {high[0] - low[0],
			   high[1] - low[1]};
  const RealType rinf   = 4.0*max(box[0], box[1]);
  const RealType cen[2] = {0.5*(low[0] + high[0]),
			   0.5*(low[1] + high[1])};
  
  // Final bounding box dimensions for quantization
  mLow.resize(2);  mHigh.resize(2);
  mLow [0] = min(low [0], cen[0]-rinf);
  mHigh[0] = max(high[0], cen[0]+rinf);
  mLow [1] = min(low [1], cen[1]-rinf);
  mHigh[1] = max(high[1], cen[1]+rinf);
  mDelta = max(degeneracy, 2.0*rinf/coordMax);
  POLY_ASSERT(mLow[0] <= mHigh[0] and mLow[1] <= mHigh[1]);

  mRinf = rinf;
  mCenter.push_back(0.5*(mLow[0] + mHigh[0]));
  mCenter.push_back(0.5*(mLow[1] + mHigh[1]));
  mCoordMax = coordMax;

  this->computeVoronoiBounded(points, PLCpoints, geometry, mesh);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Private routines called by tessellate:
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
computeCellNodesCollinear(const vector<RealType>& points,
                          map<IntPoint, pair<int,int> >& nodeMap,
                          map<int, vector<unsigned> >& cellNodes) const{
  const unsigned numGenerators = points.size()/2;
  int i;

  // Sort the genertaors but keep their original indices
  vector<pair<RealPoint,int> > pointIndexPairs;
  for (i = 0; i != numGenerators; ++i){
    pointIndexPairs.push_back(make_pair(RealPoint(points[2*i], points[2*i+1]), i));
  }
  sort( pointIndexPairs.begin(), pointIndexPairs.end(),
	internal::pairCompareFirst<RealPoint,int> );

  // Number of nodes
  const int nnodes = 2*numGenerators;

  bool test;
  unsigned inode, icell1, icell2;
  RealPoint p1, p2, r1, r2, node, midpt;
  vector<RealPoint> nodeList(nnodes);

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
  test = geometry::rayCircleIntersection(&p1.x, &r1.x, &mCenter[0],
                                         mRinf, 1.0e-10, &node.x);
  POLY_ASSERT(test);
  nodeList[inode] = node;

  // Node 1: endpt of first interior face
  test = geometry::rayCircleIntersection(&midpt.x, &r2.x, &mCenter[0],
                                         mRinf, 1.0e-10, &node.x);
  POLY_ASSERT(test);
  nodeList[inode+1] = node;
  
  // Node 2: other endpt of first interior face
  r2 *= -1.0;
  test = geometry::rayCircleIntersection(&midpt.x, &r2.x, &mCenter[0],
                                         mRinf, 1.0e-10, &node.x);
  POLY_ASSERT(test);
  nodeList[inode+2] = node;

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
    test = geometry::rayCircleIntersection(&midpt.x, &r2.x, &mCenter[0], 
                                           mRinf, 1.0e-10, &node.x);
    POLY_ASSERT(test);
    nodeList[inode] = node;
    
    // Node 1: other endpt of interior face
    r2 *= -1.0;
    test = geometry::rayCircleIntersection(&midpt.x, &r2.x, &mCenter[0],
                                           mRinf, 1.0e-10, &node.x);
    POLY_ASSERT(test);
    nodeList[inode+1] = node;

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
  
  test = geometry::rayCircleIntersection(&p2.x, &r1.x, &mCenter[0], 
                                         mRinf, 1.0e-10, &node.x);
  POLY_ASSERT(test);
  nodeList[inode] = node;
    
  // Last node for final cell
  cellNodes[icell1].push_back(inode);

  POLY_ASSERT(nodeList.size() == nnodes);
  IntPoint IntNode;
  for (i = 0; i != nnodes; ++i) {
    IntNode = IntPoint(nodeList[i].x, nodeList[i].y, mLow[0], mLow[1], mDelta);
    nodeMap[IntNode] = make_pair(i,0);
  }

  // Post-conditions
  POLY_ASSERT(nodeMap.size()   == nnodes       );
  POLY_ASSERT(cellNodes.size() == numGenerators);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
computeCellNodes(const vector<RealType>& points,
                 map<IntPoint, pair<int,int> >& nodeMap,
                 map<int, vector<unsigned> >& cellNodes) const{
  const int numGenerators = points.size()/2;
  int i, j;

  // Convert point generators to Polytope integer points
  vector<pair<IntPoint, int> > generatorToIndex(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
    generatorToIndex[i] = make_pair(IntPoint(points[2*i], points[2*i+1],
                                             mLow[0], mLow[1], mDelta), i);
  }
  
  // Sort the input points by the first element of the generator-index pair.
  // The second element provides the pre-sort generator index. Shame on you,
  // Boost, for making us do this.
  sort(generatorToIndex.begin(), generatorToIndex.end(),
       internal::pairCompareFirst<IntPoint, int> );

  // Some Boost voodoo to iterate over the first element of each generator-index pair
  typedef vector<pair<IntPoint, int> >::value_type value_type;
  boost::function<IntPoint(value_type&)> f = boost::bind(&value_type::first, _1);
  
  // The Boost.Polygon Voronoi Diagram object
  VD voronoi;
  construct_voronoi(boost::make_transform_iterator(generatorToIndex.begin(), f),
                    boost::make_transform_iterator(generatorToIndex.end(),   f),
                    &voronoi);
  POLY_ASSERT(voronoi.num_cells() == numGenerators);

  // Set the "color" of each edge and use it as a local index
  bool test;
  int sortedIndex=0, cellIndex;
  IntPoint node;
  RealPoint direction, pinf, endpt;
  map<IntPoint, int> node2id;
  for (VD::const_cell_iterator cellItr = voronoi.cells().begin(); 
       cellItr != voronoi.cells().end(); ++cellItr, ++sortedIndex) {
    const VD::edge_type* edge = cellItr->incident_edge();
    vector<unsigned> nodeChain;
    do {
      // Some pre-conditions
      POLY_ASSERT(sortedIndex == cellItr->source_index());
      POLY_ASSERT(sortedIndex <  numGenerators);
      cellIndex = generatorToIndex[sortedIndex].second;
      POLY_ASSERT(cellIndex   <  numGenerators);
      
      // // Blago!
      // cout << endl << "Cell " << cellIndex << endl;
      // cout << "   Infinite edge? " << ((edge->is_infinite()) ? "yes" : "no") << endl;
      // cout << "   Vertex 0:" << endl;
      // if (edge->vertex0()) {
      //    IntPoint vert = IntPoint(edge->vertex0()->x(),edge->vertex0()->y());
      //    cout << "      position = (" << vert.realx(mLow[0],mDelta) << ","
      //         << vert.realy(mLow[1],mDelta) << ")" << vert << endl;
      // } else {
      //    cout << "      Inf node" << endl;
      // }
      // cout << "   Vertex 1:" << endl;
      // if (edge->vertex1()) {
      //    IntPoint vert = IntPoint(edge->vertex1()->x(),edge->vertex1()->y());
      //    cout << "      position = (" << vert.realx(mLow[0],mDelta) << ","
      //         << vert.realy(mLow[1],mDelta) << ")" << vert << endl;
      // } else {
      //    cout << "      Inf node" << endl;
      // }
      // // Blago!

      // The two vertex pointers for this edge
      // NOTE: If edge is infinite, one of these pointers is null
      const VD::vertex_type* v0 = edge->vertex0();
      const VD::vertex_type* v1 = edge->vertex1();

      // Finite edge: just add vertex 0 to the cell nodes
      if (isFinite(v0,mCoordMax) and isFinite(v1,mCoordMax)) {
        node = IntPoint(v0->x(), v0->y());
        j = internal::addKeyToMap(node, node2id);
        nodeChain.push_back(j);
        if (j == node2id.size() - 1) nodeMap[node] = make_pair(j,1);
      }
      
      // Infinite edge: Determine the direction of the ray pointing to infinity.
      // Add the origin vertex of the ray and the projected point
      else if (isFinite(v0,mCoordMax) or isFinite(v1,mCoordMax)) {
        const VD::vertex_type* vfin = (isFinite(v0,mCoordMax)) ? v0 : v1;
        node = IntPoint(vfin->x(), vfin->y());
        endpt = RealPoint(node.realx(mLow[0], mDelta), node.realy(mLow[1], mDelta));
        
        // Determine the edge direction pointing to infinity
        const VD::cell_type *cell1 = edge->cell(), *cell2 = edge->twin()->cell();
        // Assume only point-generators for the time being
        POLY_ASSERT(cell1->contains_point() and cell2->contains_point());
        size_t index1 = cell1->source_index(), index2 = cell2->source_index();
        POLY_ASSERT(index1 < numGenerators and index2 < numGenerators);
        const IntPoint p1 = generatorToIndex[index1].first;
        const IntPoint p2 = generatorToIndex[index2].first;
        RealType r[2] = {p2.realx(mLow[0],mDelta) - p1.realx(mLow[0],mDelta),
                         p2.realy(mLow[1],mDelta) - p1.realy(mLow[1],mDelta)};
        if (isFinite(v0,mCoordMax)) {direction.x = -r[1];  direction.y =  r[0];}
        else                        {direction.x =  r[1];  direction.y = -r[0];}
	geometry::unitVector<2, RealType>(&direction.x);

        // // Blago!
        // cerr << "Projected Inf Node:" << endl
        //      << "   endpoint    : " << endpt << endl
        //      << "   direction   : " << direction << endl
        //      << "   cell Index 1: " << index1 << endl
        //      << "   cell Index 2: " << index2 << endl
        //      << "   point 1     : " << p1.realx(mLow[0], mDelta)
        //      << " " << p1.realy(mLow[1], mDelta) << endl
        //      << "   point 2     : " << p2.realx(mLow[0], mDelta)
        //      << " " << p2.realy(mLow[1], mDelta) << endl;
        // // Blago!

        // Project the finite vertex to the infinite shell
        test = geometry::rayCircleIntersection(&endpt.x, &direction.x,
                                               &mCenter[0], mRinf, 1.0e-8, &pinf.x);
        POLY_ASSERT(test);

        // Vertex 0 is finite, vertex 1 is the projected infNode. Add them in order
        if (isFinite(v0, mCoordMax)) {
          // Vertex 0
          j = internal::addKeyToMap(node, node2id);
          nodeChain.push_back(j);
          if (j == node2id.size() - 1) nodeMap[node] = make_pair(j,1);
          // Vertex 1
          node = IntPoint(pinf.x, pinf.y, mLow[0], mLow[1], mDelta);
          j = internal::addKeyToMap(node, node2id);
          nodeChain.push_back(j);
          if (j == node2id.size() - 1) nodeMap[node] = make_pair(j,0);
        }
        
        // Vertex 0 is the projected infNode. Only add vertex 0.
        else {
          node = IntPoint(pinf.x, pinf.y, mLow[0], mLow[1], mDelta);
          j = internal::addKeyToMap(node, node2id);
          nodeChain.push_back(j);
          if (j == node2id.size() - 1) nodeMap[node] = make_pair(j,0);
        }
      }
      edge = edge->next();
    } while (edge != cellItr->incident_edge());

    // Remove repeated node indices in the chain
    vector<unsigned>::iterator it = std::unique(nodeChain.begin(), nodeChain.end());
    nodeChain.resize(std::distance(nodeChain.begin(), it));
    if (nodeChain.front() == nodeChain.back()) nodeChain.resize(nodeChain.size()-1);

    // Some post-conditions
    POLY_ASSERT(!nodeChain.empty());
    
    // // Blago!
    // cerr << "Cell " << cellIndex << endl;
    // for (j = 0; j != nodeChain.size(); ++j) cerr << nodeChain[j] << " ";
    // cerr << endl;
    // // Blago!

    cellNodes[cellIndex] = nodeChain;
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);
  POLY_ASSERT(!nodeMap.empty());

  // // Blago!
  // for (i = 0; i != numGenerators; ++i) { 
  //   cerr << "Cell " << i << endl;
  //   int len = cellNodes[i].size();
  //   for (j = 0; j != cellNodes[i].size(); ++j) {
  //      cerr << "(" << cellNodes[i][j] << "," << cellNodes[i][(j+1)%len] << ") ";
  //   }
  //   cerr << endl << endl;
  // }
  // // Blago!
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
computeCellRings(const vector<RealType>& points,
		 const vector<RealType>& PLCpoints,
		 const PLC<2, RealType>& geometry,
		 vector<BGring>& cellRings,
		 bool performCellAdoption) const {
  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints  = PLCpoints.size()/2;
  map<IntPoint, pair<int, int> > nodeMap;
  map<int, vector<unsigned> > cellNodes;
  int i, j;

  // Check for collinearity and use the appropriate routine
  if (numGenerators == 2)
    this->computeCellNodesCollinear(points, nodeMap, cellNodes);
  else {
    bool collinear = true;
    i = 2;
    while (collinear and i != numGenerators) {
      collinear *= geometry::collinear<2,RealType>(&points[0], &points[2], &points[2*i], 1.0e-10);
      ++i;
    }

    if (collinear)
      this->computeCellNodesCollinear(points, nodeMap, cellNodes);
    else
      this->computeCellNodes(points, nodeMap, cellNodes);
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);
  POLY_ASSERT(!nodeMap.empty());

  // Quantize the PLC points
  vector<IntPoint> IntPLCPoints(numPLCpoints);
  for (i = 0; i != numPLCpoints; ++i) {
    IntPLCPoints[i] = IntPoint(PLCpoints[2*i], PLCpoints[2*i+1],
			       mLow[0], mLow[1], mDelta);
  }

  // Generate the quantized boundary to handle Boost.Geometry intersections
  BGpolygon boundary;
  constructBoostBoundary(IntPLCPoints, geometry, boundary);
  
  // Initialize the object to handle cell intersections
  Clipper2d<CoordHash> clipper(boundary);

  // Create a reverse look-up map of IDs to circumcenters
  POLY_ASSERT(!nodeMap.empty());
  map<int, IntPoint> id2nodes;
  for (map<IntPoint, pair<int,int> >::const_iterator itr = nodeMap.begin();
       itr != nodeMap.end(); ++itr) {
    i = itr->second.first;
    id2nodes[i] = itr->first;
  }
  POLY_ASSERT(id2nodes.size() == nodeMap.size());  

  // Construct the cell rings
  unsigned nodeIndex;
  IntPoint node;
  vector<BGring> orphans;
  cellRings.resize(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
    vector<IntPoint> cellBoundary;
    for (j = 0; j != cellNodes[i].size(); ++j) {
      nodeIndex = cellNodes[i][j];
      POLY_ASSERT(nodeIndex < id2nodes.size());
      node = id2nodes[nodeIndex];
      cellBoundary.push_back(node);
    }
    POLY_ASSERT(!cellBoundary.empty());
    cellBoundary.push_back(cellBoundary[0]);
    boost::geometry::assign(cellRings[i], BGring(cellBoundary.begin(), 
                                                 cellBoundary.end()));
    boost::geometry::correct(cellRings[i]);
    POLY_ASSERT(!cellRings[i].empty());
    POLY_ASSERT(cellRings[i].front() == cellRings[i].back());
    
    // Compute the boundary intersections
    clipper.clipCell(IntPoint(points[2*i], points[2*i+1], mLow[0], mLow[1], mDelta),
                     cellRings[i],
                     orphans);
    
    // Remove any repeated points
    boost::geometry::unique(cellRings[i]);
  }

  // If any orphaned cells exist, run the adoption algorithm
  // and modify the neighboring cell rings
  if (!orphans.empty()) {
    BoostOrphanage<RealType> orphanage(this);
    orphanage.adoptOrphans(points, &mLow[0], &mHigh[0], mDelta, cellRings, orphans);
  }
  
  // Post-conditions
  POLY_ASSERT(cellRings.size() == numGenerators);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
computeVoronoiUnbounded(const vector<RealType>& points,
			Tessellation<2, RealType>& mesh) const{
  const int numGenerators = points.size()/2;
  map<IntPoint, pair<int, int> > nodeMap;
  map<int, vector<unsigned> > cellNodes;
  int i;

  // Check for collinearity and use the appropriate routine
  if (numGenerators == 2)
    this->computeCellNodesCollinear(points, nodeMap, cellNodes);
  else {
    bool collinear = true;
    i = 2;
    while (collinear and i != numGenerators) {
      collinear *= geometry::collinear<2,RealType>(&points[0], &points[2], &points[2*i], 1.0e-10);
      ++i;
    }

    if (collinear)
      this->computeCellNodesCollinear(points, nodeMap, cellNodes);
    else
      this->computeCellNodes(points, nodeMap, cellNodes);
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);
  POLY_ASSERT(!nodeMap.empty());

  // Copy the quantized nodes to the final tessellation.
  int inside;
  const unsigned numNodes = nodeMap.size();
  mesh.nodes.resize(2*numNodes);
  mesh.infNodes.resize(numNodes);
  for (map<IntPoint, pair<int,int> >::const_iterator itr = nodeMap.begin();
       itr != nodeMap.end(); ++itr) {
    i = itr->second.first;
    inside = itr->second.second;
    POLY_ASSERT(i >= 0 and i < numNodes);
    POLY_ASSERT(inside == 0 or inside == 1);
    mesh.nodes[2*i  ] = itr->first.realx(mLow[0], mDelta);
    mesh.nodes[2*i+1] = itr->first.realy(mLow[1], mDelta);
    mesh.infNodes[i] = (inside == 1) ? 0 : 1;
  }
  POLY_ASSERT(mesh.infNodes.size() == mesh.nodes.size()/2);
  
  // Finish constructing the cell-face-node topology
  constructUnboundedMeshTopology(cellNodes, points, mesh);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
computeVoronoiBounded(const vector<RealType>& points,
		      const vector<RealType>& PLCpoints,
		      const PLC<2, RealType>& geometry,
		      Tessellation<2, RealType>& mesh) const {
  vector<BGring> cellRings;
  this->computeCellRings(points, PLCpoints, geometry, cellRings, false);

  constructBoundedMeshTopology(cellRings, points, PLCpoints, geometry, &mLow[0], mDelta, mesh);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Private tessellate routines
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
tessellate(const std::vector<RealType>& points,
           const std::vector<RealType>& PLCpoints,
           const PLC<2, RealType>& geometry,
           const RealType* low,
           const RealType* high,
           const RealType dx,
           Tessellation<2, RealType>& mesh) const {
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(!geometry.empty());
  POLY_ASSERT(points.size() > 0 and PLCpoints.size() > 0);
  POLY_ASSERT(points.size() % 2 == 0 and PLCpoints.size() % 2 == 0);
  POLY_ASSERT(low != 0 and high != 0);
  POLY_ASSERT(low[0] <= high[0] and low[1] <= high[1]);

  const CoordHash coordMax = (1LL << 26);
  const double degeneracy  = 1.5e-8;

  // Infinite-radius circle for inf nodes
  const RealType box[2] = {high[0] - low[0],
			   high[1] - low[1]};
  const RealType rinf   = 0.5*max(box[0], box[1]);
  const RealType cen[2] = {0.5*(low[0] + high[0]),
			   0.5*(low[1] + high[1])};
  
  // Bounding box dimensions for quantization
  mLow.resize(2);  mHigh.resize(2);
  mLow [0] = min(low [0], cen[0]-rinf);
  mHigh[0] = max(high[0], cen[0]+rinf);
  mLow [1] = min(low [1], cen[1]-rinf);
  mHigh[1] = max(high[1], cen[1]+rinf);
  mDelta = dx;

  // Store the bounding circle data
  mRinf = rinf;
  mCenter.push_back(0.5*(mLow[0] + mHigh[0]));
  mCenter.push_back(0.5*(mLow[1] + mHigh[1]));
  mCoordMax = coordMax;

  // Post-conditions
  POLY_ASSERT(dx == max(degeneracy, 2.0*rinf/coordMax));
  POLY_ASSERT(mLow[0] <= mHigh[0] and mLow[1] <= mHigh[1]);

  this->computeVoronoiBounded(points, PLCpoints, geometry, mesh);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template class BoostTessellator<double>;


} //end polytope namespace
