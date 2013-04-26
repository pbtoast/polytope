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

#include "within.hh"
#include "nearestPoint.hh"
#include "intersect.hh"

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
                             const vector<pair<IntPoint, int> >& generatorToIndex,
			     const RealType* low,
			     const RealType delta,
                             RealType* direction) {
  const VD::cell_type* cell1 = edge->cell();
  const VD::cell_type* cell2 = edge->twin()->cell();
  // Assume only point-generators for the time being
  POLY_ASSERT(cell1->contains_point() and cell2->contains_point());
  size_t index1 = cell1->source_index();
  size_t index2 = cell2->source_index();
  POLY_ASSERT(index1 < generatorToIndex.size() and 
	      index2 < generatorToIndex.size());
  const IntPoint p1 = generatorToIndex[index1].first;
  const IntPoint p2 = generatorToIndex[index2].first;
  RealType r[2] = {p2.realx(low[0],delta) - p1.realx(low[0],delta),
                   p2.realy(low[1],delta) - p1.realy(low[1],delta)};
  if (edge->vertex0()) { direction[0] = -r[1];  direction[1] =  r[0]; }
  else                 { direction[0] =  r[1];  direction[1] = -r[0]; }
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
  mLow [0] = cen[0] - rinf;
  mHigh[0] = cen[0] + rinf;
  mLow [1] = cen[1] - rinf;
  mHigh[1] = cen[1] + rinf;
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
                          map<PointType, pair<int,int> >& nodeMap,
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
  PointType IntNode;
  for (i = 0; i != nnodes; ++i) {
    // IntNode = PointType(nodeList[i].x, nodeList[i].y, mLow[0], mLow[1], mDelta);
    IntNode = PointType(nodeList[i].x, nodeList[i].y);
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
                 map<PointType, pair<int,int> >& nodeMap,
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


  // Blago!
  RealType low [2] = { numeric_limits<RealType>::max(),  numeric_limits<RealType>::max()};
  RealType high[2] = {-numeric_limits<RealType>::max(), -numeric_limits<RealType>::max()};
  for (i = 0; i != numGenerators; ++i) {
    low [0] = min(low [0], points[2*i  ]);
    low [1] = min(low [1], points[2*i+1]);
    high[0] = max(high[0], points[2*i  ]);
    high[1] = max(high[1], points[2*i+1]);
  }

  // cerr << "\nVertices:" << endl;
  for (VD::const_vertex_iterator itr = voronoi.vertices().begin();
       itr != voronoi.vertices().end(); ++itr) {
    low [0] = min(low [0], (mLow[0] + mDelta*itr->x()));
    low [1] = min(low [1], (mLow[1] + mDelta*itr->y()));
    high[0] = max(high[0], (mLow[0] + mDelta*itr->x()));
    high[1] = max(high[1], (mLow[1] + mDelta*itr->y()));
    // cerr << itr->x() << " " << itr->y() << endl;
  }

  RealType center[2]  = {0.5*(low[0] + high[0]), 0.5*(low[1] + high[1])};
  RealType box[2]     = {high[0] - low[0], high[1] - low[1]};
  const RealType rinf = 2.0*max(box[0], box[1]);

  mCenter[0] = center[0];
  mCenter[1] = center[1];
  mRinf = rinf;
  // Blago!



  // Iterate over the edges. Boost has organized them CCW around each generator.
  bool test;
  int sortedIndex=0, cellIndex;
  PointType vert, node;
  RealPoint direction, pinf, endpt, midpt;
  map<PointType, int> node2id;
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
      // if (cellIndex == 2 or cellIndex == 2350 or cellIndex == 2497) {
      // cerr << endl << "Cell " << cellIndex << endl;
      // cerr << "   Infinite edge? " << ((edge->is_infinite()) ? "yes" : "no") << endl;
      // cerr << "   Vertex 0:" << endl;
      // if (edge->vertex0()) {
      //    IntPoint vert = IntPoint(edge->vertex0()->x(),edge->vertex0()->y());
      //    cerr << "      position = (" << vert.realx(mLow[0],mDelta) << ","
      //         << vert.realy(mLow[1],mDelta) << ")" << vert << endl;
      // } else {
      //    cerr << "      Inf node" << endl;
      // }
      // cerr << "   Vertex 1:" << endl;
      // if (edge->vertex1()) {
      //    IntPoint vert = IntPoint(edge->vertex1()->x(),edge->vertex1()->y());
      //    cerr << "      position = (" << vert.realx(mLow[0],mDelta) << ","
      //         << vert.realy(mLow[1],mDelta) << ")" << vert << endl;
      // } else {
      //    cerr << "      Inf node" << endl;
      // }
      // }
      // // Blago!

      // The two vertex pointers for this edge
      // NOTE: If edge is infinite, one of these pointers is null
      const VD::vertex_type* v0 = edge->vertex0();
      const VD::vertex_type* v1 = edge->vertex1();

      // Finite edge: just add vertex 0 to the cell nodes
      // if (isFinite(v0,mCoordMax) and isFinite(v1,mCoordMax)) {
      if (v0 and v1) {
        vert = PointType(v0->x(), v0->y());
        node = PointType(vert.realx(mLow[0], mDelta), vert.realy(mLow[1], mDelta));
        j = internal::addKeyToMap(node, node2id);
        nodeChain.push_back(j);
        if (j == node2id.size() - 1) nodeMap[node] = make_pair(j,1);
      }
      
      // Infinite edge: Determine the direction of the ray pointing to infinity.
      // Add the origin vertex of the ray and the projected point
      // else if (isFinite(v0,mCoordMax) or isFinite(v1,mCoordMax)) {
      else {
        POLY_ASSERT(v0 or v1);
        const VD::vertex_type* vfin = (v0) ? v0 : v1;
        vert  = PointType(vfin->x(), vfin->y());
        node  = PointType(vert.realx(mLow[0], mDelta), vert.realy(mLow[1], mDelta));
        endpt = RealPoint(vert.realx(mLow[0], mDelta), vert.realy(mLow[1], mDelta));
        
        // Determine the edge direction pointing to infinity
        const VD::cell_type *cell1 = edge->cell(), *cell2 = edge->twin()->cell();
        // Assume only point-generators for the time being
        POLY_ASSERT(cell1->contains_point() and cell2->contains_point());
        size_t index1 = cell1->source_index(), index2 = cell2->source_index();
        POLY_ASSERT(index1 < numGenerators and index2 < numGenerators);
        const unsigned cellIndex1 = generatorToIndex[index1].second;
        const unsigned cellIndex2 = generatorToIndex[index2].second;
        RealType r[2] = {points[2*cellIndex2  ] - points[2*cellIndex1  ],
                         points[2*cellIndex2+1] - points[2*cellIndex1+1]};
        if (v0) {direction.x = -r[1];  direction.y =  r[0];}
        else    {direction.x =  r[1];  direction.y = -r[0];}
	geometry::unitVector<2, RealType>(&direction.x);
        
        // Project the finite vertex to the infinite shell
        // test = geometry::rayCircleIntersection(&endpt.x, &direction.x,
        //                                        &mCenter[0], mRinf, 1.0e-8, &pinf.x);
        test = geometry::rayCircleIntersection(&endpt.x, &direction.x,
                                               center, rinf, 1.0e-10, &pinf.x);
        POLY_ASSERT(test);

        // // Blago!
        // if (cellIndex == 2 or cellIndex == 2350 or cellIndex == 2497) {
        // cerr << "Projected Inf Node:" << endl
        //      << "   endpoint    : " << endpt << endl
        //      << "   direction   : " << direction << endl
        //      << "   cell Index 1: " << cellIndex1 << endl
        //      << "   cell Index 2: " << cellIndex2 << endl
        //      << "   point 1     : " << points[2*cellIndex1] << " " << points[2*cellIndex1+1] << endl
        //      << "   point 2     : " << points[2*cellIndex2] << " " << points[2*cellIndex2+1] << endl
        //      << "   projected pt: " << pinf << endl;
        // }
        // // Blago!
        
        
        // Vertex 0 is finite, vertex 1 is the projected infNode. Add them in order
        // if (isFinite(v0, mCoordMax)) {
        if (v0) {
          // Vertex 0
          j = internal::addKeyToMap(node, node2id);
          nodeChain.push_back(j);
          if (j == node2id.size() - 1) nodeMap[node] = make_pair(j,1);
          // Vertex 1
          node = PointType(pinf.x, pinf.y);
          j = internal::addKeyToMap(node, node2id);
          nodeChain.push_back(j);
          if (j == node2id.size() - 1) nodeMap[node] = make_pair(j,0);
        }
        
        // Vertex 0 is the projected infNode. Only add vertex 0.
        else {
          node = PointType(pinf.x, pinf.y);
          j = internal::addKeyToMap(node, node2id);
          nodeChain.push_back(j);
          if (j == node2id.size() - 1) nodeMap[node] = make_pair(j,0);
        }
      }

      edge = edge->next();
    } while (edge != cellItr->incident_edge());
    POLY_ASSERT(!nodeChain.empty());
    
    // Remove repeated node indices in the chain
    vector<unsigned>::iterator it = std::unique(nodeChain.begin(), nodeChain.end());
    nodeChain.resize(std::distance(nodeChain.begin(), it));
    if (nodeChain.front() == nodeChain.back()) nodeChain.resize(nodeChain.size()-1);
    POLY_ASSERT(!nodeChain.empty());
    
    // // Blago!
    // cerr << "\nCell " << cellIndex << endl;
    // for (j = 0; j != nodeChain.size(); ++j) cerr << nodeChain[j] << " ";
    // // Blago!

    cellNodes[cellIndex] = nodeChain;
  }

  // // Blago!
  // cerr << endl << "nodeMap:" << endl;
  // for (map<PointType, pair<int,int> >::const_iterator itr = nodeMap.begin();
  //      itr != nodeMap.end(); ++itr) {
  //   cerr << itr->first << "\t<-->\t" << itr->second.first << endl;
  // }
  // // Blago!

  POLY_ASSERT(cellNodes.size() == numGenerators);
  POLY_ASSERT(!nodeMap.empty());
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
computeCellRings(const vector<RealType>& points,
		 const vector<RealType>& PLCpoints,
		 const PLC<2, RealType>& geometry,
		 vector<IntRing>& cellRings,
		 bool performCellAdoption) const {
  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints  = PLCpoints.size()/2;
  map<PointType, pair<int, int> > nodeMap;
  map<int, vector<unsigned> > cellNodes;
  int i, j, k;

  // Check for collinearity and use the appropriate routine
  bool collinear = true;
  if (numGenerators == 2) {
    this->computeCellNodesCollinear(points, nodeMap, cellNodes);
  } else {
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
  vector<PointType> IntPLCPoints(numPLCpoints);
  for (i = 0; i != numPLCpoints; ++i) {
    // IntPLCPoints[i] = PointType(PLCpoints[2*i], PLCpoints[2*i+1],
    //                             mLow[0], mLow[1], mDelta);
    IntPLCPoints[i] = PointType(PLCpoints[2*i], PLCpoints[2*i+1]);
  }
  
  // Generate the quantized boundary to handle Boost.Geometry intersections
  BGpolygon boundary;
  constructBoostBoundary(IntPLCPoints, geometry, boundary);
  
  // Initialize the object to handle cell intersections
  // Clipper2d<CoordHash> clipper(boundary);
  Clipper2d<RealType> clipper(boundary);
  
  // Create a reverse look-up map of IDs to nodes
  POLY_ASSERT(!nodeMap.empty());
  map<int, PointType> id2nodes;
  vector<int> projectedNodes(nodeMap.size());
  for (map<PointType, pair<int,int> >::const_iterator itr = nodeMap.begin();
       itr != nodeMap.end(); ++itr) {
    i = itr->second.first;
    POLY_ASSERT(i < nodeMap.size());
    id2nodes[i] = itr->first;
    projectedNodes[i] = (itr->second.second == 0) ? 1 : 0;
  }
  POLY_ASSERT(id2nodes.size() == nodeMap.size());  

  // Construct the cell rings
  unsigned i1, i2, numIntersections;
  PointType node;
  RealPoint n1, n2, ninf;
  vector<BGring> orphans;
  cellRings.resize(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
    vector<PointType> cellBoundary;
    BGring ring;
    POLY_ASSERT(cellNodes[i].size() > 1);
    
    // // Blago!
    // cerr << "\nCell " << i << endl;
    // for (j = 0; j != cellNodes[i].size(); ++j) {
    //    i1 = cellNodes[i][j];
    //    i2 = cellNodes[i][(j+1) % cellNodes[i].size()];
    //    n1 = id2nodes[i1];
    //    n2 = id2nodes[i2];
    //    cerr << orient2d(&n1.x, &n2.x, (double*)&points[2*i]) << " ";
    // }
    // // Blago!
    
    for (j = 0; j != cellNodes[i].size(); ++j) {
      i1 = cellNodes[i][j];
      i2 = cellNodes[i][(j+1) % cellNodes[i].size()];
      POLY_ASSERT(i1 < id2nodes.size() and i2 < id2nodes.size());
      node = id2nodes[i1];
      // Add the first node of this edge
      cellBoundary.push_back(node);

      // If both nodes of this edge are outside the boundary,
      // check that the edge does not intersect the PLC
      if (!collinear and projectedNodes[i1] == 1 and projectedNodes[i2] == 1) {
        n1 = id2nodes[i1];
        n2 = id2nodes[i2];
        vector<RealType> result;
        numIntersections = intersect(&n1.x, &n2.x, numPLCpoints,
                                     &PLCpoints[0], geometry, result);
        
        // If they intersect, project a new node outward perpendicular to this
        // edge so that the clipped cell is contained in the PLC boundary
        if (numIntersections > 0) {
          POLY_ASSERT(result.size()/2 == numIntersections);
          RealType r[2]     = {n2.x - n1.x, n2.y - n1.y};
          RealType rperp[2] = {r[1]       , -r[0]      };
          geometry::unitVector<2,RealType>(rperp);
          bool test = geometry::rayCircleIntersection(&points[2*i],
                                                      rperp,
                                                      &mCenter[0],
                                                      mRinf,
                                                      1.0e-8,
                                                      &ninf.x);
          POLY_ASSERT(test);
          // node = PointType(ninf.x, ninf.y, mLow[0], mLow[1], mDelta);
          node = PointType(ninf.x, ninf.y);
          cellBoundary.push_back(node);
        }
      }
    }
    POLY_ASSERT(!cellBoundary.empty());
    cellBoundary.push_back(cellBoundary[0]);
    boost::geometry::assign(ring, BGring(cellBoundary.begin(), 
                                         cellBoundary.end()));
    // boost::geometry::correct(ring);
    // boost::geometry::unique(ring);
    // BGring simpRing;
    // boost::geometry::simplify(ring, simpRing, 2*mDelta);
    // ring = simpRing;
    POLY_ASSERT(!ring.empty());
    POLY_ASSERT(ring.front() == ring.back());
    
    // // Blago!
    // if (i == 3 or i == 2350 or i == 2497) {
    // cerr << "Pre-clipped ring " << i << endl;
    // for (typename BGring::iterator itr = ring.begin(); itr != ring.end(); ++itr) 
    //    cerr << (*itr).x << " " << (*itr).y << endl;
    // }
    // // Blago!

    // Compute the boundary intersections
    clipper.clipCell(PointType(points[2*i], points[2*i+1]), ring, orphans);
    
    // Remove any repeated points
    // boost::geometry::unique(ring);
    
    // // Blago!
    // if (i == 3 or i == 2350 or i == 2497) {
    // cerr << "Post-clipped ring " << i << endl;
    // for (typename BGring::iterator itr = ring.begin(); itr != ring.end(); ++itr) 
    //    cerr << (*itr).x << " " << (*itr).y << endl;
    // }
    // // Blago!
    
    // Quantize the clipped ring 
    vector<IntPoint> IntCellBoundary;
    for (typename BGring::iterator itr = ring.begin(); itr!= ring.end(); ++itr)
      IntCellBoundary.push_back(IntPoint((*itr).x, (*itr).y, mLow[0], mLow[1], mDelta));
    boost::geometry::assign(cellRings[i], IntRing(IntCellBoundary.begin(),
                                                  IntCellBoundary.end()));
    // boost::geometry::correct(cellRings[i]);
    boost::geometry::unique(cellRings[i]);
    POLY_ASSERT(!cellRings[i].empty());
    POLY_ASSERT(cellRings[i].front() == cellRings[i].back());

    // // Blago!
    // if (i == 3 or i == 2350 or i == 2497) {
    // cerr << "Post-quantized ring " << i << endl;
    // for (typename IntRing::iterator itr = cellRings[i].begin(); 
    //      itr != cellRings[i].end(); ++itr) 
    //    cerr << (*itr).realx(mLow[0], mDelta) << " " 
    //         << (*itr).realy(mLow[1], mDelta) << endl;
    // }
    // // Blago!
  }

  // If any orphaned cells exist, run the adoption algorithm
  // and modify the neighboring cell rings
  if (!orphans.empty() and performCellAdoption) {      
    // Quantize the orphans
    vector<IntRing> IntOrphans(orphans.size());
    for (k = 0; k != orphans.size(); ++k) {
      vector<IntPoint> IntOrphanBoundary;
      for (typename BGring::iterator itr = orphans[k].begin();
           itr!= orphans[k].end(); ++itr) {
        IntOrphanBoundary.push_back(IntPoint((*itr).x, (*itr).y, 
                                             mLow[0], mLow[1], mDelta));
      } 
      boost::geometry::assign(IntOrphans[k], 
                              IntRing(IntOrphanBoundary.begin(),
                                      IntOrphanBoundary.end()));
      boost::geometry::correct(IntOrphans[k]);
      boost::geometry::unique(IntOrphans[k]);
      POLY_ASSERT(!IntOrphans[k].empty());
      POLY_ASSERT(IntOrphans[k].front() == IntOrphans[k].back());
    }
    
    BoostOrphanage<RealType> orphanage(this);
    orphanage.adoptOrphans(points, &mLow[0], &mHigh[0], mDelta, cellRings, IntOrphans);
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
  map<PointType, pair<int, int> > nodeMap;
  map<int, vector<unsigned> > cellNodes;
  int i;

  // Check for collinearity and use the appropriate routine
  bool collinear = true;
  if (numGenerators == 2)
    this->computeCellNodesCollinear(points, nodeMap, cellNodes);
  else {
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
  for (map<PointType, pair<int,int> >::const_iterator itr = nodeMap.begin();
       itr != nodeMap.end(); ++itr) {
    i = itr->second.first;
    inside = itr->second.second;
    POLY_ASSERT(i >= 0 and i < numNodes);
    POLY_ASSERT(inside == 0 or inside == 1);
    mesh.nodes[2*i  ] = itr->first.x;//.realx(mLow[0], mDelta);
    mesh.nodes[2*i+1] = itr->first.y;//.realy(mLow[1], mDelta);
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
  vector<IntRing> cellRings;
  this->computeCellRings(points, PLCpoints, geometry, cellRings, true);

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

  // // Blago!
  // cerr << "Stored radius   = " << mRinf << endl
  //      << "Computed radius = " << rinf << endl
  //      << "Stored center   = (" << mCenter[0] << "," << mCenter[1] << ")" << endl
  //      << "Computed center = (" << cen[0] << "," << cen[1] << ")" << endl;
  // // Blago!

  // Store the bounding circle data
  mRinf = rinf;
  mCenter.resize(2);
  mCenter[0] = 0.5*(mLow[0] + mHigh[0]);
  mCenter[1] = 0.5*(mLow[1] + mHigh[1]);
  mCoordMax = coordMax;

  // // Blago!
  // cerr << "Input Bounding Box = "
  //      << "(" << low[0] << "," << high[0] << ")X"
  //      << "(" << low[1] << "," << high[1] << ")" << endl;
  // cerr << "Inner Bounding Box = "
  //      << "(" << mLow[0] << "," << mHigh[0] << ")X"
  //      << "(" << mLow[1] << "," << mHigh[1] << ")" << endl;
  // cerr << "Inner Mesh Spacing = " << mDelta << endl;
  // cerr << "Bounding radius    = " << mRinf << endl;
  // cerr << "Circle center      = "
  //      << "(" << mCenter[0] << "," << mCenter[1] << ")" << endl;
  // // Blago!  

  // Post-conditions
  POLY_ASSERT2(dx == max(degeneracy, 2.0*rinf/coordMax),
               "\nInput spacing    : " << dx << 
               "\nComputed spacing : " << (2.0*rinf/coordMax));
  POLY_ASSERT(mLow[0] <= mHigh[0] and mLow[1] <= mHigh[1]);

  vector<IntRing> cellRings;
  this->computeCellRings(points, PLCpoints, geometry, cellRings, false);

  constructBoundedMeshTopology(cellRings, points, PLCpoints, 
                               geometry, &mLow[0], mDelta, mesh);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template class BoostTessellator<double>;


} //end polytope namespace
