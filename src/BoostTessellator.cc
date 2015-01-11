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
#include "polytope_plc_canned_geometries.hh"
#include "PLC_Boost_2d.hh"

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
// Comparison operator for two points. They're equal if one lives inside the
// 3x3 region around the other:    _ _ _
//                                |_|_|_|
//                                |_|_|_|
//                                |_|_|_|
//------------------------------------------------------------------------------
template<typename RealType>
class ThreeByThreeCompare {
public:
  bool operator()(const Point2<RealType> pt1, const Point2<RealType> pt2) const {
    return (pt1.x < pt2.x-1 ? true :
            (pt1.x == pt2.x-1 or pt1.x == pt2.x or pt1.x == pt2.x+1) and 
            pt1.y < pt2.y-1 ? true : false);
  }
};

//------------------------------------------------------------------------------
// Comparison operator for two points. They're equal if one lives inside the
// 3x3 region around the other:    _ _ _
//                                |_|_|_|
//                                |_|_|_|
//                                |_|_|_|
//------------------------------------------------------------------------------
template<typename RealType>
class ThreeByThreeTolCompare {
public:
  ThreeByThreeTolCompare(RealType tol_) : tol(tol_) {};
  bool operator()(const Point2<RealType> pt1, const Point2<RealType> pt2) const {
    return (pt1.x < pt2.x-tol ? true :
            (pt1.x >= pt2.x-tol and pt1.x <= pt2.x+tol) and 
            pt1.y < pt2.y-tol ? true : false);
  }
private:
  RealType tol;
};


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
  mCoords.initialize(points, mDegeneracy);
  
  // Check for collinear generators
  const bool collinear = geometry::collinear<2, RealType>(points, 1.0e-10);
  
  // Use the appropriate cell node routine
  vector<vector<unsigned> > cellNodes;
  map<int, PointType> id2node;
  vector<unsigned> infNodes;
  if (collinear) 
  {
    this->computeCellNodesCollinear(points, cellNodes, id2node, infNodes);
  }
  else 
  {
    this->computeCellNodes(points, cellNodes, id2node, infNodes);
  }
  POLY_ASSERT(cellNodes.size() == points.size()/2);

  // Copy the quantized nodes to the final tessellation.
  const unsigned numNodes = id2node.size();
  RealPoint node;
  mesh.nodes.resize(2*numNodes);
  mesh.infNodes = infNodes;
  for (typename map<int, PointType>::const_iterator itr = id2node.begin();
       itr != id2node.end();
       ++itr) {
    const unsigned inode = itr->first;
    POLY_ASSERT(inode >= 0 and inode < numNodes);
    node = BTT::dequantize(mCoords, itr->second);
    mesh.nodes[2*inode  ] = node.x;
    mesh.nodes[2*inode+1] = node.y;
  }

  // Finish constructing the cell-face-node topology
  constructUnboundedMeshTopology(cellNodes, points, mesh);
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
  POLY_ASSERT(low != 0 and high != 0);
  POLY_ASSERT(low[0] <= high[0] and low[1] <= high[1]);
  
  // Build a PLC with the bounding box, and then use the PLC method.
  ReducedPLC<2, RealType> box = plc_box<2, RealType>(low, high);
  this->tessellate(points, box, mesh);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
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
BoostTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const ReducedPLC<2, RealType>& geometry,
           Tessellation<2, RealType>& mesh) const {
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(not geometry.points.empty());

  // Initialize quantized coordinate system
  mCoords.initialize(geometry.points, mDegeneracy);

  // Check for collinear generators
  const bool collinear = geometry::collinear<2, RealType>(points, 1.0e-10);
  const unsigned numGenerators = points.size()/2;

  // Use the appropriate cell node routine
  vector<vector<unsigned> > cellNodes;
  map<int, PointType> id2node;
  vector<unsigned> infNodes;
  if (collinear) 
  {
    this->computeCellNodesCollinear(points, cellNodes, id2node, infNodes);
  }
  else 
  {
    this->computeCellNodes(points, cellNodes, id2node, infNodes);
  }
  POLY_ASSERT(cellNodes.size() == numGenerators);

  vector<ReducedPLC<2, CoordType> > orphans;
  vector<ReducedPLC<2, CoordType> > cells(numGenerators);
  ReducedPLC<2, CoordType> boundary;
  boundary.facets = geometry.facets;
  boundary.holes  = geometry.holes;
  boundary.points.resize(geometry.points.size());
  for (int j = 0; j < geometry.points.size()/2; ++j) {
    const PointType pp = BTT::quantize(mCoords, &geometry.points[2*j]);
    boundary.points[2*j  ] = pp.x;
    boundary.points[2*j+1] = pp.y;
  }
  for (int i = 0; i < numGenerators; ++i) {
    const PointType generator = BTT::quantize(mCoords, &points[2*i]);
    ReducedPLC<2, CoordType> cell = plcOfCell<CoordType>(cellNodes[i], id2node);


    // // Blago!
    // RealType pta[2] = {points[2*i], points[2*i+1]};
    // RealType orientation = 0.0;
    // for (int ifacet = 0; ifacet < cell.facets.size(); ++ifacet) {
    //    RealType ptb[2] = {cell.points[2*cell.facets[ifacet][0]  ],
    //                       cell.points[2*cell.facets[ifacet][0]+1]};
    //    RealType ptc[2] = {cell.points[2*cell.facets[ifacet][1]  ],
    //                       cell.points[2*cell.facets[ifacet][1]+1]};
    //    orientation += orient2d(pta, ptb, ptc);
    // }
    // if (orientation <= 0.0) {
    //    cerr << i << " " << orientation << "\n" << cell << endl;
    //    std::reverse(cell.points.begin(), cell.points.end());
    //    cerr << i << " " << orientation << "\n" << cell << endl;
    // }
    // // POLY_ASSERT2(orientation > 0.0,
    // //              i << " " << orientation << "\n" << cell);
    // // Blago!

    //TODO
    //TODO implement the bi-infinite edge fix
    //TODO
    
    cells[i] = BG::boost_clip<CoordType>(boundary, cell, generator, orphans);
    POLY_ASSERT(not BG::boost_intersects<CoordType>(cells[i]));
  }

  if (not orphans.empty()) {
    cerr << "Orphans detected. Taking no actions." << endl;
  }

  // Input nodes and construct the final mesh topology
  this->constructBoundedTopology(points, geometry, cells, mesh);
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
                 vector<vector<unsigned> >& cellNodes,
                 map<int, PointType>& id2node,
                 vector<unsigned>& infNodes) const{
  const int numGenerators = points.size()/2;
  // The Boost.Polygon Voronoi Diagram object
  // VB voroBuilder;
  VD voronoi;

  // Convert point generators to Polytope integer points
  vector<IntPoint> generators(numGenerators);
  for (unsigned i = 0; i < numGenerators; ++i) {
    IntPoint ip = mCoords.quantize(&points[2*i]);
    ip.index = i;
    generators[i] = ip;
  }
  
  // Sort the input points by the first element of the generator-index pair.
  // The second element provides the pre-sort generator index. Shame on you,
  // Boost, for making us do this.
  sort(generators.begin(), generators.end());
  
  // for (unsigned i = 0; i < numGenerators; ++i) {
  //   voroBuilder.insert_point(generators[i].x, generators[i].y);
  // }
  // voroBuilder.construct(&voronoi);
  construct_voronoi(generators.begin(), generators.end(), &voronoi);
  
  // // Blago!
  // {
  //   int sortedIndex = 0;
  //   cerr << "ccc = []" << endl;
  //   for (typename VD::const_cell_iterator cellItr = voronoi.cells().begin(); 
  //        cellItr != voronoi.cells().end(); 
  //        ++cellItr, ++sortedIndex) {
  //     // cerr << endl << "Cell " << cellItr->source_index() << " of " 
  //     //      << voronoi.num_cells() << endl
  //     //      << "genrat index : " << generators[sortedIndex].index << endl
  //     //      << "sorted index : " << sortedIndex << endl
  //     //      << "int generator: " << generators[sortedIndex] << endl;
  //     cerr << "cc" << cellItr->source_index() << " = [" << endl;
  //     const typename VD::edge_type* edge = cellItr->incident_edge();
  //     do {
  //       const typename VD::vertex_type* v0 = edge->vertex0();
  //       const typename VD::vertex_type* v1 = edge->vertex1();
  //       if (v0 and v1) {
  //         RealPoint vert = RealPoint(v0->x(), v0->y());
  //         PointType node = BTT::deboost(mCoords, &vert.x);
  //         cerr << "(" << node.x << "," << node.y << ")," << endl;
  //       } else {
  //         const typename VD::vertex_type* vfin = v0 ? v0 : v1;
  //         RealPoint vert = RealPoint(vfin->x(), vfin->y());
  //         PointType node = BTT::deboost(mCoords, &vert.x);
  //         PointType endpt = node;
  //         const typename VD::cell_type* cell1 = edge->cell();
  //         const typename VD::cell_type* cell2 = edge->twin()->cell();
  //         const size_t index1 = cell1->source_index();
  //         const size_t index2 = cell2->source_index();
  //         const unsigned cellIndex1 = generators[index1].index;
  //         const unsigned cellIndex2 = generators[index2].index;
  //         RealPoint r = RealPoint(points[2*cellIndex2  ] - points[2*cellIndex1  ],
  //                                 points[2*cellIndex2+1] - points[2*cellIndex1+1]);
  //         RealPoint direction = v0 ? RealPoint(-r.y, r.x) : RealPoint(r.y, -r.x);
  //         geometry::unitVector<2, RealType>(&direction.x);
  //         PointType pinf = BTT::project(mCoords, endpt, &direction.x);
  //         RealPoint vpinf =  mCoords.projectPoint(&vert.x, &direction.x);
  //         if (v0) {
  //           cerr << "(" << node.x << "," << node.y << ")," << endl;
  //           cerr << "(" << pinf.x << "," << pinf.y << ")," << endl;
  //         }
  //         else {
  //           cerr << "(" << pinf.x << "," << pinf.y << ")," << endl;
  //         }
  //       }
  //       edge = edge->next();
  //     } while (edge != cellItr->incident_edge());
  //     cerr << "]" << endl;
  //     cerr << "ccc.append(cc" << cellItr->source_index() << ")" << endl;
  //   }
  // }
  // // Blago!

  POLY_ASSERT(voronoi.num_cells() == numGenerators);

  // Compute a bounding box for the floating point vertex positions
  RealPoint vlow  = RealPoint( numeric_limits<RealType>::max(),  
                               numeric_limits<RealType>::max());
  RealPoint vhigh = RealPoint(-numeric_limits<RealType>::max(), 
			      -numeric_limits<RealType>::max());
  for (typename VD::const_vertex_iterator itr = voronoi.vertices().begin();
       itr != voronoi.vertices().end(); 
       ++itr) {
    RealPoint rp = RealPoint(itr->x(), itr->y());
    PointType pp = BTT::deboost(mCoords, &rp.x);
    RealPoint p  = BTT::dequantize(mCoords, pp);
    vlow.x = min(vlow.x, p.x);  vhigh.x = max(vhigh.x, p.x);
    vlow.y = min(vlow.y, p.y);  vhigh.y = max(vhigh.y, p.y);
  }
  
  // Expand the bounding box to include
  mCoords.expand(&vlow.x, &vhigh.x);
  
  // Iterate over the edges. Boost has organized them CCW around each generator.
  int sortedIndex=0, cellIndex;
  PointType node, pinf, endpt;
  RealPoint direction, vert;
  map<PointType, int> node2id;
  // map<PointType, int, ThreeByThreeCompare<CoordType> > node2id;
  // map<PointType, int, ThreeByThreeTolCompare<CoordType> > 
  //    node2id(ThreeByThreeTolCompare<CoordType>(1.49e-8));
  cellNodes.resize(numGenerators);
  for (typename VD::const_cell_iterator cellItr = voronoi.cells().begin(); 
       cellItr != voronoi.cells().end(); 
       ++cellItr, ++sortedIndex) {
    const typename VD::edge_type* edge = cellItr->incident_edge();
    vector<unsigned> nodeChain;
    bool boundedCell = true;
    int infBegin = -1;
    int infEnd = -1;
    POLY_CONTRACT_VAR(infEnd);
    do {
      // Some pre-conditions
      POLY_ASSERT2(sortedIndex == cellItr->source_index(),
                   sortedIndex << " != " << cellItr->source_index());
      POLY_ASSERT(sortedIndex <  numGenerators);
      cellIndex = generators[sortedIndex].index;
      POLY_ASSERT(cellIndex   <  numGenerators);      

      // The two vertex pointers for this edge
      // NOTE: If edge is infinite, one of these pointers is null
      const typename VD::vertex_type* v0 = edge->vertex0();
      const typename VD::vertex_type* v1 = edge->vertex1();

      // Finite edge: just add vertex 0 to the cell nodes
      if (v0 and v1) {
        vert = RealPoint(v0->x(), v0->y());
        node = BTT::deboost(mCoords, &vert.x);
        const unsigned old_size = node2id.size();
        const unsigned j = internal::addKeyToMap(node, node2id);
        nodeChain.push_back(j);
        if (j == old_size) id2node[j] = node;
      }
      
      // Infinite edge: Determine the direction of the ray pointing to infinity.
      // Add the origin vertex of the ray and the projected point
      else {
        POLY_ASSERT(v0 or v1);
        boundedCell *= false;
        const typename VD::vertex_type* vfin = v0 ? v0 : v1;
        vert  = RealPoint(vfin->x(), vfin->y());
        node  = BTT::deboost(mCoords, &vert.x);
        endpt = node;
        
        // Determine the edge direction pointing to infinity
        const typename VD::cell_type* cell1 = edge->cell();
        const typename VD::cell_type* cell2 = edge->twin()->cell();

        // Assume only point-generators for the time being
        POLY_ASSERT(cell1->contains_point() and cell2->contains_point());
        const size_t index1 = cell1->source_index();
        const size_t index2 = cell2->source_index();
        POLY_ASSERT(index1 < numGenerators and index2 < numGenerators);
        const unsigned cellIndex1 = generators[index1].index;
        const unsigned cellIndex2 = generators[index2].index;
        RealPoint r = RealPoint(points[2*cellIndex2  ] - points[2*cellIndex1  ],
                                points[2*cellIndex2+1] - points[2*cellIndex1+1]);
        direction = v0 ? RealPoint(-r.y, r.x) : RealPoint(r.y, -r.x);
	geometry::unitVector<2, RealType>(&direction.x);
        
        // Project the finite vertex to the infinite shell
        pinf = BTT::project(mCoords, endpt, &direction.x);
        
        // Vertex 0 is finite, vertex 1 is the projected infNode. Add them in order
        if (v0) {
          { // Vertex 0
            const unsigned old_size = node2id.size();
            const unsigned j = internal::addKeyToMap(node, node2id);
            nodeChain.push_back(j);
            if (j == old_size)  id2node[j] = node;
          }

          { // Vertex 1
            node = pinf;
            const unsigned old_size = node2id.size();
            const unsigned j = internal::addKeyToMap(node, node2id);
            nodeChain.push_back(j);
            if (j == old_size)  id2node[j] = node;
            if (j == old_size)  infNodes.push_back(j);
            infEnd = j;
          }
        }
        
        // Vertex 0 is the projected infNode. Only add vertex 0.
        else {
          node = pinf;
          const unsigned old_size = node2id.size();
          const unsigned j = internal::addKeyToMap(node, node2id);
          nodeChain.push_back(j);
          if (j == old_size) id2node[j] = node;
          if (j == old_size) infNodes.push_back(j);
          infBegin = j;
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

    if (not boundedCell) {
      POLY_ASSERT(infBegin != -1 and infEnd != -1);
      POLY_ASSERT(std::find(nodeChain.begin(), nodeChain.end(), infBegin) != nodeChain.end());
      POLY_ASSERT(std::find(nodeChain.begin(), nodeChain.end(), infEnd  ) != nodeChain.end());
      vector<unsigned>::iterator beginIt = std::find(nodeChain.begin(), nodeChain.end(), infBegin);
      std::rotate(nodeChain.begin(), beginIt, nodeChain.end());
      POLY_ASSERT(nodeChain.front() == infBegin);
      POLY_ASSERT(nodeChain.back()  == infEnd  );
      
      // Edge 1
      const unsigned n11 = *(nodeChain.begin()  );     //Projected
      const unsigned n12 = *(nodeChain.begin()+1);
      POLY_ASSERT(n11 != n12);
      POLY_ASSERT(id2node.find(n11) != id2node.end());
      POLY_ASSERT(id2node.find(n12) != id2node.end());
      RealPoint rp11 = BTT::dequantize(mCoords, id2node[n11]);
      RealPoint rp12 = BTT::dequantize(mCoords, id2node[n12]);

      // Edge 2
      const unsigned n21 = *(nodeChain.end()-2);
      const unsigned n22 = *(nodeChain.end()-1);       //Projected
      POLY_ASSERT(n21 != n22);
      POLY_ASSERT(id2node.find(n21) != id2node.end());
      POLY_ASSERT(id2node.find(n22) != id2node.end());
      RealPoint rp21 = BTT::dequantize(mCoords, id2node[n21]);
      RealPoint rp22 = BTT::dequantize(mCoords, id2node[n22]);

      // Compute self-intersections
      RealPoint result;
      bool intersects;
      if (n12 == n21) intersects = false;
      else            intersects = geometry::segmentIntersection2D(&rp11.x, &rp12.x,
                                                                   &rp21.x, &rp22.x,
                                                                   &result.x,
                                                                   mDegeneracy);

      if (intersects) {
        // // Blago!
        // cerr << "Self-intersection detected:" << endl
        //      << "  Cell " << cellIndex << endl
	//      << "  " << rp11 << "  " << rp12 << endl
	//      << "  " << rp21 << "  " << rp22 << endl
	//      << "  Intersection pt = " << result << endl;
        // // Blago!

        if (geometry::distance<2,RealType>(&rp12.x, &result.x) > mDegeneracy and
            geometry::distance<2,RealType>(&rp21.x, &result.x) > mDegeneracy) {
          const PointType newNode = BTT::quantize(mCoords, &result.x);
          id2node[n11] = newNode;
          id2node[n22] = newNode;
        }
      }
    }

    // // Blago!
    // if (cellIndex == 53) {
    //    cerr << "Cell " << cellIndex << ":" << endl;
    //    for (unsigned j = 0; j < nodeChain.size(); ++j) {
    //       cerr << "  " << j << ": " << nodeChain[j] << " " << id2node[nodeChain[j]] << endl;
    //    }
    //    cerr << endl;
    // }
    // // Blago!

    cellNodes[cellIndex] = nodeChain;
  }

  POLY_ASSERT(cellNodes.size() == numGenerators);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
computeCellNodesCollinear(const vector<RealType>& points,
                          vector<vector<unsigned> >& cellNodes,
                          map<int, PointType>& id2node,
                          vector<unsigned>& infNodes) const{

  // A dummy expansion just in case we're doing a distributed mesh construction.
  mCoords.expand(&mCoords.low.x, &mCoords.high.x);

  // Call the 1d routine for projecting a line of points
  vector<RealPoint> nodes;
  constructCells1d(points, &(mCoords.center).x, mCoords.rinf, cellNodes, nodes);
  POLY_ASSERT(cellNodes.size() == points.size()/2);
  POLY_ASSERT(nodes.size() == points.size());

  // Quantize nodes and assign indices
  std::set<PointType> uniqueNodes;
  infNodes.resize(nodes.size());
  for (unsigned i = 0; i < nodes.size(); ++i) {
    const PointType p = BTT::quantize(mCoords, &(nodes[i]).x);
    uniqueNodes.insert(p);
    POLY_ASSERT(uniqueNodes.size() == i+1);
    id2node[i] = p;
    infNodes[i] = i;   // All nodes are projected inf nodes
  }
  POLY_ASSERT(uniqueNodes.size() == nodes.size());
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
constructBoundedTopology(const vector<RealType>& points,
			 const ReducedPLC<2, RealType>& geometry,
			 const vector<ReducedPLC<2, CoordType> >& cellRings,
			 Tessellation<2, RealType>& mesh) const {
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  
  // Now build the unique mesh nodes and cell info.
  const unsigned numGenerators = cellRings.size();
  // map<PointType, int> point2node;
  map<PointType, int, ThreeByThreeTolCompare<CoordType> > 
     point2node(ThreeByThreeTolCompare<CoordType>(1.49e-8));
  // map<IntPoint, int, ThreeByThreeCompare<CoordHash> > point2node;
  map<EdgeHash, int> edgeHash2id;
  map<int, vector<int> > edgeCells;
  int i, j, k, iedge;
  mesh.cells = std::vector<std::vector<int> >(numGenerators);
  for (i = 0; i != numGenerators; ++i) { 
    POLY_ASSERT(cellRings[i].facets.size() > 2);
    for (unsigned ifacet = 0; ifacet < cellRings[i].facets.size(); ++ifacet) {
      POLY_ASSERT(cellRings[i].facets[ifacet].size() == 2);
      const int i1 = cellRings[i].facets[ifacet][0];
      const int i2 = cellRings[i].facets[ifacet][1];
      POLY_ASSERT(i1 != i2);
      const PointType p1 = PointType(cellRings[i].points[2*i1], cellRings[i].points[2*i1+1]);
      const PointType p2 = PointType(cellRings[i].points[2*i2], cellRings[i].points[2*i2+1]);
      // const PointType pp1 = PointType(cellRings[i].points[2*i1], cellRings[i].points[2*i1+1]);
      // const PointType pp2 = PointType(cellRings[i].points[2*i2], cellRings[i].points[2*i2+1]);
      // const RealPoint rp1 = BTT::dequantize(mCoords, pp1);
      // const RealPoint rp2 = BTT::dequantize(mCoords, pp2);
      // const IntPoint  p1  = mCoords.quantize(&rp1.x);
      // const IntPoint  p2  = mCoords.quantize(&rp2.x);
      POLY_ASSERT(p1 != p2);
      j = internal::addKeyToMap(p1, point2node);
      k = internal::addKeyToMap(p2, point2node);
      POLY_ASSERT(j != k);
      // if (j != k) {
        iedge = internal::addKeyToMap(internal::hashEdge(j,k), edgeHash2id);
        edgeCells[iedge].push_back(j < k ? i : ~i);
        POLY_ASSERT2(edgeCells[iedge].size() == 1 or edgeCells[iedge][0]*edgeCells[iedge][1] <= 0,
                     "BLAGO: " << iedge << " " << j << " " << k << " " << edgeCells[iedge][0] 
                     << " " << edgeCells[iedge][1]);
        mesh.cells[i].push_back(j < k ? iedge : ~iedge);
      // }
    }
    POLY_ASSERT(mesh.cells[i].size() >= 3);
  }
  POLY_ASSERT(edgeCells.size() == edgeHash2id.size());
  
  // Fill in the mesh nodes.
  RealPoint node;
  mesh.nodes = std::vector<RealType>(2*point2node.size());
  for (typename std::map<PointType, int>::const_iterator itr = point2node.begin();
  // for (typename std::map<IntPoint, int>::const_iterator itr = point2node.begin();
       itr != point2node.end(); 
       ++itr) {
    const RealPoint p = BTT::dequantize(mCoords, itr->first);
    // const RealPoint p = mCoords.dequantize(&(itr->first).x);
    i = itr->second;
    POLY_ASSERT(i < mesh.nodes.size()/2);
    // if (not BG::boost_within(p, geometry)) {
    //   RealPoint result;
    //   RealType dist = nearestPoint(&node.x, 
    // 				      geometry.points.size()/2,
    // 				      &geometry.points[0],
    // 				      geometry,
    // 				      &result.x);
    //   POLY_ASSERT(dist < 1.0e-10);
    //   node = result;
    // }
    node = p;
    mesh.nodes[2*i]   = node.x;
    mesh.nodes[2*i+1] = node.y;
  }
  
  // Fill in the mesh faces.
  mesh.faces = std::vector<std::vector<unsigned> >(edgeHash2id.size());
  for (typename std::map<EdgeHash, int>::const_iterator itr = edgeHash2id.begin();
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

  // Post-conditions
  POLY_ASSERT(mesh.faceCells.size() == mesh.faces.size());
  POLY_ASSERT(mesh.cells.size()     == numGenerators    );
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
           const ReducedPLC<2, CoordHash>& intGeometry,
           const QuantizedCoordinates<2, RealType>& coords,
           vector<ReducedPLC<2, CoordHash> >& intCells) const {
  // Pre-conditions
  POLY_ASSERT(not intGeometry.empty());
  POLY_ASSERT(not points.empty() > 0);
  POLY_ASSERT(points.size() % 2 == 0);
  POLY_ASSERT(not coords.empty());

  // The Quantized coordinates
  mCoords = coords;
  
  const unsigned numGenerators = points.size()/2;
  const bool collinear = geometry::collinear<2, RealType>(points, mDegeneracy);
  vector<vector<unsigned> > cellNodes;
  map<int, PointType> id2node;
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

  // vector<ReducedPLC<2, CoordHash> > dummy;
  // intCells.resize(numGenerators);
  // for (int i = 0; i < numGenerators; ++i) {
  //    const IntPoint intGenerator 
  // }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template class BoostTessellator<double>;


} //end polytope namespace
