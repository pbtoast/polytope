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

#include "polytope.hh" // Pulls in POLY_ASSERT
#include "within.hh"
#include "nearestPoint.hh"
#include "intersect.hh"

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
			     const RealType* endpt,
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
  RealType midpt[2] = {0.5*(p1.realx(low[0],delta) + p2.realx(low[0],delta)), 
		       0.5*(p1.realy(low[1],delta) + p2.realy(low[1],delta))};
  direction[0] = midpt[0] - endpt[0];
  direction[1] = midpt[1] - endpt[1];
  geometry::unitVector<2, double>(direction);
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
  mLow [0] = min(low [0], cen[0]-rinf);
  mHigh[0] = max(high[0], cen[0]+rinf);
  mLow [1] = min(low [1], cen[1]-rinf);
  mHigh[1] = max(high[1], cen[1]+rinf);
  mDelta = max(degeneracy, 2.0*rinf/coordMax);
  POLY_ASSERT(mLow[0] <= mHigh[0] and mLow[1] <= mHigh[1]);

  mRinf = rinf;
  mCenter.push_back(0.5*(mLow[0] + mHigh[0]));
  mCenter.push_back(0.5*(mLow[1] + mHigh[1]));

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
BoostTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const vector<RealType>& PLCpoints,
           const PLC<2, RealType>& geometry,
           Tessellation<2, RealType>& mesh) const 
{
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
computeCellRings(const vector<RealType>& points,
		 const vector<RealType>& PLCpoints,
		 const PLC<2, RealType>& geometry,
		 vector<BGring>& cellRings,
		 bool performCellAdoption) const {
  const int numGenerators = points.size()/2;
  const int numPLCpoints  = PLCpoints.size()/2;
  int i, j, k;

  // Quantize the PLC points
  vector<IntPoint> IntPLCPoints(numPLCpoints);
  for (i = 0; i != numPLCpoints; ++i) {
    IntPLCPoints[i] = IntPoint(PLCpoints[2*i], PLCpoints[2*i+1],
			       mLow[0], mLow[1], mDelta);
  }

  // Generate the quantized boundary to handle Boost.Geometry intersections
  BGpolygon boundary;
  buildBoostBoundary(IntPLCPoints, geometry, boundary);
  
  // Convert point generators to Polytope integer points
  vector<IntPoint> generators(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
    generators[i] = IntPoint(points[2*i], points[2*i+1],
			     mLow[0], mLow[1], mDelta);
  }

  // Sort the input points
  // NOTE: Boost.Polygon will do this by default.
  // TODO: Store the pre-sorted indices and store the tessellation cell
  //       info in terms of these original indices
  sort(generators.begin(), generators.end());

  // The Boost.Polygon Voronoi Diagram object
  VD voronoi;
  construct_voronoi(generators.begin(), generators.end(), &voronoi);

  // Build up the ring of integer nodes around each generator directly from 
  bool test, inside;
  int nodeIndex, cellIndex=0;
  RealPoint direction, pinf, endpt;
  IntPoint X, finiteVertex, node;
  size_t colorIndex = 1;
  map<VD::vertex_type, int, vertexCompare> vertexToNodeIndex;
  vector<BGring> orphanage;
  cellRings.resize(numGenerators);
  for (VD::const_cell_iterator cellItr = voronoi.cells().begin(); 
       cellItr != voronoi.cells().end(); ++cellItr, ++cellIndex) {
    vector<IntPoint> cellBoundary;
    const VD::edge_type* edge = cellItr->incident_edge();
    do {
      POLY_ASSERT2(colorIndex < numeric_limits<size_t>::max()/32, "BoostTessellator Error: "
                   << "overflow of maximum allowable face index");

      // Check if the edge on the opposite side of the face has been given a "color".
      // If it hasn't add the vertex info of this edge.
      //if (edge->twin()->color() == 0) {
        edge->color(colorIndex);
        ++colorIndex;
        
        // Input cell-to-face info into mesh.cells
        POLY_ASSERT(cellIndex < voronoi.num_cells());
        POLY_ASSERT(cellItr->source_index() == cellIndex);
        
	// The two vertex pointers for this edge
        // NOTE: If edge is infinite, one of these pointers is null
        const VD::vertex_type* v0 = edge->vertex0();
        const VD::vertex_type* v1 = edge->vertex1();

        // Finite edge: Add both vertices to the map
        if (edge->is_finite()) {
	  POLY_ASSERT(v0 and v1);
	  // vertex 0
//           k = vertexToNodeIndex.size();
//           nodeIndex = internal::addKeyToMap(*v0, vertexToNodeIndex);
//           if (nodeIndex == k) {
	    node = IntPoint(v0->x(), v0->y());
	    cellBoundary.push_back(node);
//           }
          
// 	  // vertex 1
//           k = vertexToNodeIndex.size();
//           nodeIndex = internal::addKeyToMap(*v1, vertexToNodeIndex);
//           if (nodeIndex == k) {
// 	    node = IntPoint(v1->x(), v1->y());
// 	    cellBoundary.push_back(node);
//           }
        }

        // Infinite edge: Determine the direction of the ray pointing to infinity.
        // Add the origin vertex of the ray and the projected point
        else {
          finiteVertex = (edge->vertex0())                       ?
	    IntPoint(edge->vertex0()->x(), edge->vertex0()->y()) :
	    IntPoint(edge->vertex1()->x(), edge->vertex1()->y());
	  endpt = RealPoint(finiteVertex.realx(mLow[0],mDelta),
			    finiteVertex.realy(mLow[1],mDelta));
	  computeInfiniteEdgeDirection(edge, generators, &endpt.x, &mLow[0], mDelta, &direction.x);
	  test = geometry::rayCircleIntersection(&endpt.x, &direction.x,
                                                 &mCenter[0], mRinf, 1.0e-6, &pinf.x);
          POLY_ASSERT(test);

          // Vertex 0 exists, vertex 1 is a projected infNode. Add them in order
          if (v0) {
	    // vertex 0
//             k = vertexToNodeIndex.size();
//             nodeIndex = internal::addKeyToMap(*v0, vertexToNodeIndex);
//             if (nodeIndex == k) {
              node = IntPoint(v0->x(), v0->y());
	      cellBoundary.push_back(node);
//             }

	    // infinite vertex 1
//             k = vertexToNodeIndex.size();
//             nodeIndex = internal::addKeyToMap(VD::vertex_type(pinf.x,pinf.y),
//                                               vertexToNodeIndex);
//             if (nodeIndex == k) {
	      node = IntPoint(pinf.x, pinf.y, mLow[0], mLow[1], mDelta);
	      cellBoundary.push_back(node);
//             }
          } 
          
          // Vertex 0 is a projected infNode, vertex 1 exists. Add them in order
          else {
	    // infinite vertex 0
//             k = vertexToNodeIndex.size();
//             nodeIndex = internal::addKeyToMap(VD::vertex_type(pinf.x,pinf.y),
//                                               vertexToNodeIndex);
//             if (nodeIndex == k) {
	      node = IntPoint(pinf.x, pinf.y, mLow[0], mLow[1], mDelta);
	      cellBoundary.push_back(node);
//             }
            
// 	    // vertex 1
// 	    k = vertexToNodeIndex.size();
//             nodeIndex = internal::addKeyToMap(*v1, vertexToNodeIndex);
//             if (nodeIndex == k) {
//               node = IntPoint(v1->x(), v1->y());
// 	      cellBoundary.push_back(node);
//             }
          }
	  //}
      }
      edge = edge->next();
    } while (edge != cellItr->incident_edge());
    POLY_ASSERT(cellBoundary.size() > 0);
    cellBoundary.push_back(cellBoundary[0]);
    boost::geometry::assign(cellRings[cellIndex], BGring(cellBoundary.begin(),
							 cellBoundary.end()));
    POLY_ASSERT(cellRings[cellIndex].size() > 0);

    // Blago!
    cerr << "\nCell " << cellIndex << endl;
    for (typename BGring::const_iterator itr = cellRings[cellIndex].begin();
         itr != cellRings[cellIndex].end(); ++itr) {
      cerr << (*itr).realx(mLow[0],mDelta) << " "
           << (*itr).realy(mLow[1],mDelta) << endl;
    }
    // Blago!
    

    // Intersect with the boundary to get the bounded cell.
    // Since for complex boundaries this may return more than one polygon, we find
    // the one that contains the generator.
    vector<BGring> cellIntersections;
    boost::geometry::intersection(boundary, cellRings[cellIndex], cellIntersections);
    if (cellIntersections.size() == 0) {
      cerr << points[2*cellIndex] << " " << points[2*cellIndex+1] << endl 
           << boost::geometry::dsv(cellRings[cellIndex]) << endl
           << boost::geometry::dsv(boundary) << endl;
    }
    POLY_ASSERT(cellIntersections.size() > 0);
    if (cellIntersections.size() == 1) {
      cellRings[cellIndex] = cellIntersections[0];
    } else {
      X = IntPoint(points[2*cellIndex], points[2*cellIndex+1], 
		   mLow[0], mLow[1], mDelta);
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
      cellRings[cellIndex] = cellIntersections[k];
    }
  }

  // If any orphaned cells exist, run the adoption algorithm
  // and modify the neighboring cell rings
  if (performCellAdoption and orphanage.size() > 0){
    POLY_ASSERT2(false, "Adoption has not been implemented for "
		 << "the BoostTessellator yet");
  }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
computeVoronoiUnbounded(const vector<RealType>& points,
			Tessellation<2, RealType>& mesh) const{
  const int numGenerators = points.size()/2;
  int i;

  // Convert point generators to Polytope integer points
  vector<IntPoint> generators(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
    generators[i] = IntPoint(points[2*i], points[2*i+1],
			     mLow[0], mLow[1], mDelta);
  }

  // Sort the input points
  // NOTE: Boost.Polygon will do this by default.
  // TODO: Store the pre-sorted indices and store the tessellation cell
  //       info in terms of these original indices
  sort(generators.begin(), generators.end());

  // The Boost.Polygon Voronoi Diagram object
  VD voronoi;
  construct_voronoi(generators.begin(), generators.end(), &voronoi);

  // Set the "color" of each edge and use it as a local index
  bool test;
  size_t colorIndex = 1;
  int k, nodeIndex, faceIndex, cellIndex=0;
  RealPoint direction, pinf, endpt;
  IntPoint finiteVertex, node;
  set<int> infCells;
  map<VD::vertex_type, int, vertexCompare> vertexToNodeIndex;
  mesh.cells.resize(voronoi.num_cells());
  for (VD::const_cell_iterator cellItr = voronoi.cells().begin(); 
       cellItr != voronoi.cells().end(); ++cellItr, ++cellIndex) {
    const VD::edge_type* edge = cellItr->incident_edge();
    do {
      POLY_ASSERT2(colorIndex < numeric_limits<size_t>::max()/32, "BoostTessellator Error: "
                   << "overflow of maximum allowable face index");

      // Store the indices of all infinite cells
      if (edge->is_infinite())  infCells.insert(cellIndex);


      // Blago!
      cout << endl << "Edge index " << colorIndex << endl;
      cout << "   Infinite edge? " << ((edge->is_infinite()) ? "yes" : "no") << endl;
      cout << "   Is it primary? " << ((edge->is_primary()) ? "yes" : "no" ) << endl;
      cout << "   Vertex 0:" << endl;
      if (edge->vertex0()) {
         IntPoint vert = IntPoint(edge->vertex0()->x(),edge->vertex0()->y());
         cout << "      position = (" << vert.realx(mLow[0],mDelta) << ","
              << vert.realy(mLow[1],mDelta) << ")" << endl;
      } else {
         cout << "      Inf node" << endl;
      }
      cout << "   Vertex 1:" << endl;
      if (edge->vertex1()) {
         IntPoint vert = IntPoint(edge->vertex1()->x(),edge->vertex1()->y());
         cout << "      position = (" << vert.realx(mLow[0],mDelta) << ","
              << vert.realy(mLow[1],mDelta) << ")" << endl;
      } else {
         cout << "      Inf node" << endl;
      }
      // Blago!

      if (edge->twin()->color() == 0) {
        faceIndex = colorIndex-1;
        edge->color(colorIndex);
        ++colorIndex;
        
        // Input cell-to-face info into mesh.cells
        POLY_ASSERT(cellIndex < voronoi.num_cells());
        POLY_ASSERT(cellItr->source_index() == cellIndex);
        mesh.cells[cellIndex].push_back(faceIndex);
        size_t oppCellIndex = edge->twin()->cell()->source_index();
        POLY_ASSERT(oppCellIndex < voronoi.num_cells());
        mesh.cells[oppCellIndex].push_back(~faceIndex);
        
        // Input face-to-cell info info mesh.faceCells
        mesh.faceCells.push_back(vector<int>(2));
        POLY_ASSERT(faceIndex < mesh.faceCells.size());
        POLY_ASSERT(mesh.faceCells[faceIndex].size() == 2);
        mesh.faceCells[faceIndex][0] = cellIndex;
        mesh.faceCells[faceIndex][1] = ~(int(oppCellIndex));
        
        // Size the face and infFace arrays
        mesh.faces.push_back(vector<unsigned>(2));
        mesh.infFaces.push_back(0);
        POLY_ASSERT(faceIndex < mesh.faces.size());
        POLY_ASSERT(mesh.faces[faceIndex].size() == 2);

        // The two vertex pointers for this edge
        // NOTE: If edge is infinite, one of these pointers is null
        const VD::vertex_type* v0 = edge->vertex0();
        const VD::vertex_type* v1 = edge->vertex1();

        // Finite edge: Add both vertices to the map
        if (edge->is_finite()) {
	  POLY_ASSERT(v0 and v1);
	  // vertex 0
          k = vertexToNodeIndex.size();
          nodeIndex = internal::addKeyToMap(*v0, vertexToNodeIndex);
          if (nodeIndex == k) {
	    node = IntPoint(v0->x(), v0->y());
            mesh.nodes.push_back(node.realx(mLow[0], mDelta));
            mesh.nodes.push_back(node.realy(mLow[1], mDelta));
            mesh.infNodes.push_back(0);
          }
          mesh.faces[faceIndex][0] = nodeIndex;
          
	  // vertex 1
          k = vertexToNodeIndex.size();
          nodeIndex = internal::addKeyToMap(*v1, vertexToNodeIndex);
          if (nodeIndex == k) {
	    node = IntPoint(v1->x(), v1->y());
            mesh.nodes.push_back(node.realx(mLow[0], mDelta));
            mesh.nodes.push_back(node.realy(mLow[1], mDelta));
            mesh.infNodes.push_back(0);
          }
          mesh.faces[faceIndex][1] = nodeIndex;
        }

        // Infinite edge: Determine the direction of the ray pointing to infinity.
        // Add the origin vertex of the ray and the projected point
        else {
          finiteVertex = (edge->vertex0())              ?
	    IntPoint(edge->vertex0()->x(), edge->vertex0()->y()) :
	    IntPoint(edge->vertex1()->x(), edge->vertex1()->y());
	  endpt = RealPoint(finiteVertex.realx(mLow[0],mDelta),
			    finiteVertex.realy(mLow[1],mDelta));
	  computeInfiniteEdgeDirection(edge, generators, &endpt.x, &mLow[0], mDelta, &direction.x);
	  test = geometry::rayCircleIntersection(&endpt.x, &direction.x,
                                                 &mCenter[0], mRinf, 1.0e-6, &pinf.x);
          POLY_ASSERT(test);

          // Vertex 0 exists, vertex 1 is a projected infNode. Add them in order
          if (v0) {
	    // vertex 0
            k = vertexToNodeIndex.size();
            nodeIndex = internal::addKeyToMap(*v0, vertexToNodeIndex);
            if (nodeIndex == k) {
              node = IntPoint(v0->x(), v0->y());
              mesh.nodes.push_back(node.realx(mLow[0], mDelta));
              mesh.nodes.push_back(node.realy(mLow[1], mDelta));
	      mesh.infNodes.push_back(0);
            }
            mesh.faces[faceIndex][0] = nodeIndex;

	    // infinite vertex 1
            k = vertexToNodeIndex.size();
            nodeIndex = internal::addKeyToMap(VD::vertex_type(pinf.x,pinf.y),
                                              vertexToNodeIndex);
            if (nodeIndex == k) {
              mesh.nodes.push_back(pinf.x);
              mesh.nodes.push_back(pinf.y);
              mesh.infNodes.push_back(1);
            }
            mesh.faces[faceIndex][1] = nodeIndex;
          } 
          
          // Vertex 0 is a projected infNode, vertex 1 exists. Add them in order
          else {
	    // infinite vertex 0
            k = vertexToNodeIndex.size();
            nodeIndex = internal::addKeyToMap(VD::vertex_type(pinf.x,pinf.y),
                                              vertexToNodeIndex);
            if (nodeIndex == k) {
              mesh.nodes.push_back(pinf.x);
              mesh.nodes.push_back(pinf.y);
              mesh.infNodes.push_back(1);
            }
            mesh.faces[faceIndex][0] = nodeIndex;
            
	    // vertex 1
	    k = vertexToNodeIndex.size();
            nodeIndex = internal::addKeyToMap(*v1, vertexToNodeIndex);
            if (nodeIndex == k) {
              node = IntPoint(v1->x(), v1->y());
              mesh.nodes.push_back(node.realx(mLow[0], mDelta));
              mesh.nodes.push_back(node.realy(mLow[1], mDelta));
	      mesh.infNodes.push_back(0);
            }
            mesh.faces[faceIndex][1] = nodeIndex;
          }
        }
      }
      edge = edge->next();
    } while (edge != cellItr->incident_edge());
  }
  
  // Use the "color" index as a counter for the number of faces
  int numFaces = colorIndex-1;
  // A few sanity checks
  POLY_ASSERT(mesh.faces.size()     == numFaces);
  POLY_ASSERT(mesh.faceCells.size() == numFaces);
  POLY_ASSERT(mesh.infNodes.size()  == mesh.nodes.size()/2);
  POLY_ASSERT(mesh.infFaces.size()  == numFaces);

  // Blago!
  cerr << "Inf cells: ";
  for (set<int>::const_iterator itr = infCells.begin();
       itr != infCells.end(); ++itr) {
     cerr << *itr << " ";
  }
  cerr << endl;
  // Blago!

  // Infinite cells by convention are bounded by inf faces. These connect
  // the inf nodes around an unbounded cell.
  for (set<int>::const_iterator it = infCells.begin(); 
       it != infCells.end(); ++it){
     int icell = *it;
     POLY_ASSERT(mesh.cells[icell].size() > 1);
     vector<int>::iterator faceItr0 = mesh.cells[icell].end()-1;
     vector<int>::iterator faceItr1 = mesh.cells[icell].begin();
     bool lookingForInfFace = true;
     unsigned infNode0, infNode1;
     while (lookingForInfFace and faceItr1 != mesh.cells[icell].end()) {
        const unsigned iface0 = *faceItr0 < 0 ? ~(*faceItr0) : *faceItr0;
        const unsigned iface1 = *faceItr1 < 0 ? ~(*faceItr1) : *faceItr1;
        POLY_ASSERT(iface0 < mesh.faces.size() and iface1 < mesh.faces.size());
        POLY_ASSERT(mesh.faces[iface0].size() == 2 and mesh.faces[iface1].size() == 2);
        const unsigned inode00 = *faceItr0 < 0 ? mesh.faces[iface0][1] : mesh.faces[iface0][0];
        const unsigned inode01 = *faceItr0 < 0 ? mesh.faces[iface0][0] : mesh.faces[iface0][1];
        const unsigned inode10 = *faceItr1 < 0 ? mesh.faces[iface1][1] : mesh.faces[iface1][0];
        const unsigned inode11 = *faceItr1 < 0 ? mesh.faces[iface1][0] : mesh.faces[iface1][1];
        POLY_ASSERT(inode00 < mesh.nodes.size()/2 and inode01 < mesh.nodes.size()/2 and
                    inode10 < mesh.nodes.size()/2 and inode11 < mesh.nodes.size()/2 );
        if (mesh.infNodes[inode00] == 0 and 
            mesh.infNodes[inode01] == 1 and
            mesh.infNodes[inode10] == 1 and
            mesh.infNodes[inode11] == 0) {
          POLY_ASSERT(inode01 != inode10);
          lookingForInfFace = false;
          infNode0 = inode01;
          infNode1 = inode10;
        }
        faceItr0 = faceItr1;
        ++faceItr1;
     }
     POLY_ASSERT(!lookingForInfFace);
     vector<unsigned> faceNodes(2);
     faceNodes[0] = infNode0;
     faceNodes[1] = infNode1;
     mesh.faces.push_back(faceNodes);
     mesh.faceCells.push_back(vector<int>(1, icell));
     mesh.infFaces.push_back(1);
     mesh.cells[icell].insert(faceItr0, numFaces);
     ++numFaces;
  }
  // A few repeat sanity checks
  POLY_ASSERT(mesh.faces.size()     == numFaces);
  POLY_ASSERT(mesh.faceCells.size() == numFaces);
  POLY_ASSERT(mesh.infFaces.size()  == numFaces);
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
  const int numGenerators = points.size()/2;
  const int numPLCpoints  = PLCpoints.size()/2;
  int i, j, k;
  
  vector<BGring> cellRings;
  this->computeCellRings(points, PLCpoints, geometry, cellRings, false);

  // Now build the unique mesh nodes and cell info.
  map<IntPoint, int> point2node;
  map<EdgeHash, int> edgeHash2id;
  map<int, vector<int> > edgeCells;
  int iedge;
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
    node[0] = p.realx(mLow[0],mDelta);
    node[1] = p.realy(mLow[1],mDelta);
    
    // Check if nodes are inside boundary (either bounding box or PLC, if defined)
    bool inside = within(node, numPLCpoints, &PLCpoints[0], geometry);
    
    if( !inside ){
      RealType result[2];
      RealType dist = nearestPoint( node, numPLCpoints, &PLCpoints[0], geometry, result );
      // Check the node has not moved more than 2.5 quantized mesh spacings. NOTE: this is
      // not a sharp estimate. Theoreticallly, the distance ought to be at most sqrt(2)*cdx, 
      // but nodes will fail this strict of a test.
      POLY_ASSERT2( dist < 5.0*mDelta,
                    dist << " " << 2.5*mDelta << " : "
                    << "(" << node[0]   << " " << node[1]   << ") "
                    << "(" << result[0] << " " << result[1] << ")\n" << geometry );
      node[0] = result[0];
      node[1] = result[1];
    }
     
    mesh.nodes[2*i]   = node[0];
    mesh.nodes[2*i+1] = node[1];
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
      for (j = 0; j != edgeCells[i].size(); ++j) {
	cerr << " --> " << edgeCells[i][j] << " " 
	     << points[2*edgeCells[i][j]] << " " 
	     << points[2*edgeCells[i][j]+1] << endl;
      }
    }
    POLY_ASSERT(edgeCells[i].size() == 1 or edgeCells[i].size() == 2);
    mesh.faceCells[i] = edgeCells[i];
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class BoostTessellator<double>;
} //end polytope namespace
