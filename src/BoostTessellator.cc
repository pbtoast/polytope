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
                             const QuantizedCoordinates<2, RealType>& coords,
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
  const RealPoint r = coords.dequantize(p2) - coords.dequantize(p1);
  // RealType r[2] = {p2.realx(low[0],delta) - p1.realx(low[0],delta),
  //                  p2.realy(low[1],delta) - p1.realy(low[1],delta)};
  // if (edge->vertex0()) { direction[0] = -r[1];  direction[1] =  r[0]; }
  // else                 { direction[0] =  r[1];  direction[1] = -r[0]; }
  if (edge->vertex0()) { direction[0] = -r.y;  direction[1] =  r.x; }
  else                 { direction[0] =  r.y;  direction[1] = -r.x; }
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

  // Initialize quantized coordinate system
  mCoords.initialize(points);
  
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
  POLY_ASSERT(points.size() % 2 == 0);
  
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

  // Initialize quantized coordinate system
  mCoords.initialize(PLCpoints);

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
computeCellNodes(const vector<RealType>& points,
                 map<PointType, pair<int,int> >& nodeMap,
                 vector<vector<unsigned> >& cellNodes) const{
  const int numGenerators = points.size()/2;
  int i, j;

  // Convert point generators to Polytope integer points
  vector<pair<IntPoint, int> > generatorToIndex(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
     generatorToIndex[i] = make_pair(mCoords.quantize(&points[2*i]), i);
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

  // Compute a bounding box for the floating point vertex positions
  RealType vlow[2], vhigh[2];
  for (VD::const_vertex_iterator itr = voronoi.vertices().begin();
       itr != voronoi.vertices().end(); ++itr) {
    vlow [0] = min(mCoords.low [0], (mCoords.low[0] + mCoords.delta*itr->x()));
    vlow [1] = min(mCoords.low [1], (mCoords.low[1] + mCoords.delta*itr->y()));
    vhigh[0] = max(mCoords.high[0], (mCoords.low[0] + mCoords.delta*itr->x()));
    vhigh[1] = max(mCoords.high[1], (mCoords.low[1] + mCoords.delta*itr->y()));
  }

  //mCoords.expand(vlow, vhigh);
  
  RealType vcenter[2]  = {0.5*(vlow[0] + vhigh[0]), 0.5*(vlow[1] + vhigh[1])};
  RealType vbox[2]     = {vhigh[0] - vlow[0], vhigh[1] - vlow[1]};
  const RealType vrinf = 2.0*max(vbox[0], vbox[1]);

  // mCenter[0] = center[0];
  // mCenter[1] = center[1];
  // mRinf = rinf;
  // Blago!

  // Iterate over the edges. Boost has organized them CCW around each generator.
  int sortedIndex=0, cellIndex;
  PointType vert, node;
  RealPoint direction, pinf, endpt, midpt;
  map<PointType, int> node2id;
  cellNodes.resize(numGenerators);
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
      // cerr << endl << "Cell " << cellIndex << endl;
      // cerr << "   Infinite edge? " << ((edge->is_infinite()) ? "yes" : "no") << endl;
      // cerr << "   Vertex 0:" << endl;
      // if (edge->vertex0()) {
      //    IntPoint vert = IntPoint(edge->vertex0()->x(),edge->vertex0()->y());
      //    cerr << "      position = (" << mCoords.dequantize(vert) << endl;
      // } else {
      //    cerr << "      Inf node" << endl;
      // }
      // cerr << "   Vertex 1:" << endl;
      // if (edge->vertex1()) {
      //    IntPoint vert = IntPoint(edge->vertex1()->x(),edge->vertex1()->y());
      //    cerr << "      position = (" << mCoords.dequantize(vert) << endl;
      // } else {
      //    cerr << "      Inf node" << endl;
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
        node  = PointType(vert.realx(mCoords.low[0], mCoords.delta), 
                          vert.realy(mCoords.low[1], mCoords.delta));
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
        node  = PointType(vert.realx(mCoords.low[0], mCoords.delta), 
                          vert.realy(mCoords.low[1], mCoords.delta));
        endpt = PointType(vert.realx(mCoords.low[0], mCoords.delta), 
                          vert.realy(mCoords.low[1], mCoords.delta));
        
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
        //pinf = mCoords.projectPoint(&endpt.x, &direction.x);
        bool test = geometry::rayCircleIntersection(&endpt.x, &direction.x,
                                                    vcenter, vrinf, 1.0e-10, &pinf.x);
        POLY_ASSERT(test);

        // // Blago!
        // cerr << "Projected Inf Node:" << endl
        //      << "   endpoint    : " << endpt << endl
        //      << "   direction   : " << direction << endl
        //      << "   cell Index 1: " << cellIndex1 << endl
        //      << "   cell Index 2: " << cellIndex2 << endl
        //      << "   point 1     : " << points[2*cellIndex1] << " " << points[2*cellIndex1+1] << endl
        //      << "   point 2     : " << points[2*cellIndex2] << " " << points[2*cellIndex2+1] << endl
        //      << "   projected pt: " << pinf << endl;
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
    if (nodeChain.front() == nodeChain.back())// and nodeChain.size() != 1) 
       nodeChain.resize(nodeChain.size()-1);
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
		 const map<PointType, pair<int, int> >& nodeMap,
                 vector<vector<unsigned> >& cellNodes,
		 Clipper2d<RealType>& clipper,
                 vector<IntRing>& cellRings,
                 const bool collinear,
		 const bool performCellAdoption) const {
  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints = PLCpoints.size()/2;
  int i, j, k;
  
  // // Quantize the PLC points
  // const unsigned numPLCpoints = PLCpoints.size()/2;
  // vector<PointType> IntPLCPoints(numPLCpoints);
  // for (i = 0; i != numPLCpoints; ++i) {
  //   IntPLCPoints[i] = PointType(PLCpoints[2*i], PLCpoints[2*i+1]);
  // }
  
  // // Generate the quantized boundary to handle Boost.Geometry intersections
  // BGpolygon boundary;
  // constructBoostBoundary(IntPLCPoints, geometry, boundary);
  
  // // Initialize the object to handle cell intersections
  // Clipper2d<RealType> clipper(boundary);
  
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
          ninf = mCoords.projectPoint(&points[2*i], rperp);
          node = PointType(ninf.x, ninf.y);
          cellBoundary.push_back(node);
        }
      }
    }
    POLY_ASSERT(!cellBoundary.empty());
    cellBoundary.push_back(cellBoundary[0]);
    boost::geometry::assign(ring, BGring(cellBoundary.begin(), 
                                         cellBoundary.end()));
    POLY_ASSERT(!ring.empty());
    POLY_ASSERT(ring.front() == ring.back());
    
    // Compute the boundary intersections
    clipper.clipCell(PointType(points[2*i], points[2*i+1]), ring, orphans);
    
    // Quantize the clipped ring 
    vector<IntPoint> IntCellBoundary;
    for (typename BGring::iterator itr = ring.begin(); itr!= ring.end(); ++itr)
      IntCellBoundary.push_back(mCoords.quantize(&(*itr).x));
  
    boost::geometry::assign(cellRings[i], IntRing(IntCellBoundary.begin(),
                                                  IntCellBoundary.end()));
    boost::geometry::unique(cellRings[i]);
    POLY_ASSERT(!cellRings[i].empty());
    POLY_ASSERT(cellRings[i].front() == cellRings[i].back());
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
        IntOrphanBoundary.push_back(mCoords.quantize(&(*itr).x));
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
    orphanage.adoptOrphans(points, mCoords, cellRings, IntOrphans);
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
    vector<PointType> nodeList;
    computeCellNodesCollinear(points, mCoords, nodeList, cellNodes);
    for (i = 0; i != nodeList.size(); ++i) {
      nodeMap[nodeList[i]] = make_pair(i,0);
    }
  }
  else {
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
    mesh.nodes[2*i  ] = itr->first.x;
    mesh.nodes[2*i+1] = itr->first.y;
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
  POLY_ASSERT(!points.empty());
  POLY_ASSERT(mesh.empty());

  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints  = PLCpoints.size()/2;
  map<PointType, pair<int, int> > nodeMap;
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
    vector<PointType> nodeList;
    computeCellNodesCollinear(points, mCoords, nodeList, cellNodes);
    for (i = 0; i != nodeList.size(); ++i) {
      nodeMap[nodeList[i]] = make_pair(i,1);
    }
  }
  else {
    this->computeCellNodes(points, nodeMap, cellNodes);
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);
  POLY_ASSERT(!nodeMap.empty());

  // Quantize the PLC points
  vector<PointType> IntPLCPoints(numPLCpoints);
  for (i = 0; i != numPLCpoints; ++i) {
    IntPLCPoints[i] = PointType(PLCpoints[2*i], PLCpoints[2*i+1]);
  }
  
  // Generate the quantized boundary to handle Boost.Geometry intersections
  BGpolygon boundary;
  constructBoostBoundary(IntPLCPoints, geometry, boundary);
  
  // Initialize the object to handle cell intersections
  Clipper2d<RealType> clipper(boundary);

  // Compute bounded cell rings
  vector<IntRing> cellRings;
  this->computeCellRings(points, PLCpoints, geometry, nodeMap, 
                         cellNodes, clipper, cellRings, collinear, true);

  // Input nodes and construct the final mesh topology
  constructBoundedMeshTopology(cellRings, points, mCoords, mesh);
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
           const std::vector<CoordHash>& IntPLCpoints,
           const PLC<2, RealType>& geometry,
           const QuantizedCoordinates<2, RealType>& coords,
           vector<vector<vector<CoordHash> > >& IntCells) const {
  // Pre-conditions
  POLY_ASSERT(!geometry.empty());
  POLY_ASSERT(points.size() > 0 and IntPLCpoints.size() > 0);
  POLY_ASSERT(points.size() % 2 == 0 and IntPLCpoints.size() % 2 == 0);
  POLY_ASSERT(!coords.empty());

  // The Quantized coordinates
  mCoords = coords;
  
  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints = IntPLCpoints.size()/2;
  map<PointType, pair<int, int> > nodeMap;
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
    vector<PointType> nodeList;
    computeCellNodesCollinear(points, mCoords, nodeList, cellNodes);
    for (i = 0; i != nodeList.size(); ++i) {
      nodeMap[nodeList[i]] = make_pair(i,1);
    }
  }
  else {
    this->computeCellNodes(points, nodeMap, cellNodes);
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);
  POLY_ASSERT(!nodeMap.empty());
  
  // Store the input boundary as a Boost.Geometry polygon
  BGpolygon boundary;
  vector<PointType> boundaryPoints(numPLCpoints);
  vector<RealType> RealPLCpoints;
  for (i = 0; i != numPLCpoints; ++i) {
     //boundaryPoints[i] = PointType(IntPLCpoints[2*i], IntPLCpoints[2*i+1]);
    RealPoint pt = mCoords.dequantize(IntPoint(IntPLCpoints[2*i], IntPLCpoints[2*i+1]));
    boundaryPoints[i] = pt;
    RealPLCpoints.push_back(pt.x);
    RealPLCpoints.push_back(pt.y);
  }
  constructBoostBoundary(boundaryPoints, geometry, boundary);
  
  // Initialize the object to handle cell intersections
  Clipper2d<RealType> clipper(boundary);

  // Compute bounded cell rings
  vector<IntRing> cellRings;
  this->computeCellRings(points, RealPLCpoints, geometry, nodeMap, 
                         cellNodes, clipper, cellRings, collinear, false);

  // Store the rings in a non-Boost way
  IntCells.resize(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
    vector<CoordHash> node(2);
    int index = 0;
    IntCells[i].resize(cellRings[i].size());
    for (typename IntRing::const_iterator itr = cellRings[i].begin();
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
template class BoostTessellator<double>;


} //end polytope namespace
