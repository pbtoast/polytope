#ifndef POLYTOPE_TESSELLATOR_UTILITIES_HH
#define POLYTOPE_TESSELLATOR_UTILITIES_HH

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/algorithms/unique.hpp>

#include "Point.hh"
#include "QuantizedCoordinates.hh"
#if HAVE_BOOST
BOOST_GEOMETRY_REGISTER_POINT_2D(polytope::Point2<int64_t>, int64_t, boost::geometry::cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_POINT_2D(polytope::Point2<double>, double, boost::geometry::cs::cartesian, x, y);
#endif


namespace polytope {

//------------------------------------------------------------------------------
// constructUnboundedMeshTopology
//
// Constructs the mesh structures for cells, faces, faceCells, and infFaces.
// INPUT: 
//   cellNodes: map from cell Index to collection of node indices
// PRE-CONDITIONS:
//   The node and infNode vectors have already been set
//------------------------------------------------------------------------------
template<typename RealType>
void
constructUnboundedMeshTopology(std::vector<std::vector<unsigned> >& cellNodes,
                               const std::vector<RealType> points,
                               Tessellation<2,RealType>& mesh) {  
  // Pre-conditions
  POLY_ASSERT(!mesh.nodes.empty());
  POLY_ASSERT(mesh.cells.empty() and mesh.faceCells.empty() and 
              mesh.faces.empty() and mesh.infFaces.empty());
  
  // Typedefs
  typedef std::pair<int, int> EdgeHash;

  const unsigned numGenerators = cellNodes.size();

  unsigned i, j, i1, i2, iface;
  EdgeHash face;
  std::vector<unsigned> faceVec(2);
  std::map<EdgeHash, int> face2id;
  internal::CounterMap<unsigned> faceCounter;
  mesh.cells.resize(numGenerators);
  for (i = 0; i != numGenerators; ++i) {
    // Check the orientation of the nodes around the generator. Reverse the order
    // if they're ordered clockwise. It is enough to check the order of the first
    // two nodes in the cellNodes list, since its node list is sequential
    POLY_ASSERT(cellNodes[i].size() > 2);
    if( orient2d(&mesh.nodes[2*cellNodes[i][0]], 
		 &mesh.nodes[2*cellNodes[i][1]], 
		 (double*)&points[2*i]) < 0 ){
      reverse(cellNodes[i].begin(), cellNodes[i].end());
    }
    
    // The ordered mesh nodes around a given generator
    for (j = 0; j != cellNodes[i].size(); ++j){
      i1 = cellNodes[i][j];
      i2 = cellNodes[i][(j+1) % cellNodes[i].size()];
      POLY_ASSERT(i1 != i2);
      face  = internal::hashEdge(i1,i2);
      iface = internal::addKeyToMap(face,face2id);
      ++faceCounter[iface];
      
      // If you're looking at this face for the first time, add its
      // nodes, classify it as an interior or inf face, and resize faceCells
      POLY_ASSERT( faceCounter[iface] > 0 );
      if( faceCounter[iface] == 1 ){
        faceVec[0] = face.first; faceVec[1] = face.second;
        mesh.faces.push_back(faceVec);
        mesh.faceCells.resize(iface+1);
	mesh.infFaces.resize(iface+1);
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

  // Post-conditions
  POLY_ASSERT(mesh.faceCells.size() == mesh.faces.size()  );
  POLY_ASSERT(mesh.infFaces.size()  == mesh.faces.size()  );
  POLY_ASSERT(mesh.infNodes.size()  == mesh.nodes.size()/2);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// constructBoundedMeshTopology
//
// Constructs the mesh structures for cells, faces, faceCells, and infFaces.
// INPUT: 
//    cellRings: Boost.Geometry ring of quantized node positions for each cell
//    PLCpoints: Input boundary points
//    geometry : PLC boundary
//    coords   : Quantized coordinate system and bounding box
// PRE-CONDITIONS: 
//    The mesh is empty
//------------------------------------------------------------------------------
template<typename RealType>
void
constructBoundedMeshTopology(const std::vector<boost::geometry::model::ring
                               <polytope::Point2<int64_t>, false> >& cellRings,
                             const std::vector<RealType>& points,
                             //const std::vector<RealType>& PLCpoints,
                             //const PLC<2,RealType>& geometry,
                             const QuantizedCoordinates<2, RealType>& coords,
                             Tessellation<2,RealType>& mesh) {
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(!coords.empty());
  
  // Typedefs
  typedef int64_t CoordHash;
  typedef std::pair<int, int> EdgeHash;
  typedef polytope::Point2<CoordHash> IntPoint;
  typedef polytope::Point2<double> RealPoint;
  typedef boost::geometry::model::ring<IntPoint, false> BGring;

  const unsigned numGenerators = cellRings.size();
  // const unsigned numPLCpoints  = PLCpoints.size();

  // Now build the unique mesh nodes and cell info.
  std::map<IntPoint, int> point2node;
  std::map<EdgeHash, int> edgeHash2id;
  std::map<int, std::vector<int> > edgeCells;
  int i, j, k, iedge;
  mesh.cells = std::vector<std::vector<int> >(numGenerators);
  for (i = 0; i != numGenerators; ++i) { 
    POLY_ASSERT(cellRings[i].size() > 2);
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
  RealPoint node;
  // RealType node[2];
  mesh.nodes = std::vector<RealType>(2*point2node.size());
  for (typename std::map<IntPoint, int>::const_iterator itr = point2node.begin();
       itr != point2node.end(); ++itr) {
    const IntPoint& p = itr->first;
    i = itr->second;
    POLY_ASSERT(i < mesh.nodes.size()/2);
    node = coords.dequantize(&p.x);
    
    // // Check if nodes are inside boundary (either bounding box or PLC, if defined)
    // bool inside = within(&node.x, numPLCpoints, &PLCpoints[0], geometry);
    // inside = true;
    // if(!inside){
    //   RealType result[2];
    //   RealType dist = nearestPoint(&node.x, numPLCpoints, &PLCpoints[0], geometry, result);
    //   // Check the node has not moved more than 2.5 quantized mesh spacings. NOTE: this is
    //   // not a sharp estimate. Theoreticallly, the distance ought to be at most sqrt(2)*dx, 
    //   // but nodes will fail this strict of a test.
    //   // POLY_ASSERT2( dist < 5.0*dx,
    //   //               dist << " " << 2.5*dx << " : "
    //   //               << "(" << node[0]   << " " << node[1]   << ") "
    //   //               << "(" << result[0] << " " << result[1] << ")\n" << geometry );
    //   node.x = result[0];
    //   node.y = result[1];
    // }

    mesh.nodes[2*i]   = node.x;
    mesh.nodes[2*i+1] = node.y;
  }
  
  // Fill in the mesh faces.
  mesh.faces = std::vector<std::vector<unsigned> >(edgeHash2id.size());
  for (typename std::map<EdgeHash, int>::const_iterator itr = edgeHash2id.begin();
       itr != edgeHash2id.end(); ++itr) {
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
// constructBoostBoundary
//
// Store PointType PLC boundary data as a Boost.Geometry polygon
// INPUT: 
//    PLCpoints: vector of PointType data for the vertices of the PLC
//    geometry : The boundary PLC
//    boundary : Ref to the output Boost.Geometry polygon
// PRE-CONDITIONS: 
//    The PointType must first be registered with Boost
//------------------------------------------------------------------------------
template<typename RealType, typename PointType>
void
constructBoostBoundary(const std::vector<PointType>& PLCPoints,
                       const PLC<2,RealType>& geometry,
                       boost::geometry::model::polygon<PointType,false>& boundary) {
  typedef boost::geometry::model::polygon<PointType,false> BGpolygon;
  int i, j, k;
  boost::geometry::append( boundary, PLCPoints[geometry.facets[0][0]] );
  for (j = 0; j != geometry.facets.size(); ++j){
    POLY_ASSERT(geometry.facets[j].size() == 2);
    i =  geometry.facets[j][1];
    boost::geometry::append( boundary, PLCPoints[i]);
  }
  POLY_ASSERT(boundary.outer().size() == geometry.facets.size() + 1);
  POLY_ASSERT(boundary.outer().front() == boundary.outer().back());
  
  // Add any interior holes.
  const unsigned numHoles = geometry.holes.size();
  if (numHoles > 0) {
    typename BGpolygon::inner_container_type& holes = boundary.inners();
    holes.resize(numHoles);
    for (k = 0; k != numHoles; ++k) {
      boost::geometry::append( holes[k], PLCPoints[geometry.holes[k][0][0]] );
      for (j = 0; j != geometry.holes[k].size(); ++j) {
	POLY_ASSERT(geometry.holes[k][j].size() == 2);
	i =  geometry.holes[k][j][1];
	boost::geometry::append( holes[k], PLCPoints[i]); 
      }
      POLY_ASSERT(holes[k].size() == geometry.holes[k].size() + 1 );
      POLY_ASSERT(holes[k].front() == holes[k].back());
    }
    POLY_ASSERT(boundary.inners().size() == numHoles );
  }
}

//------------------------------------------------------------------------------
// convertTessellationToRings
//
// Convert mesh cells to Boost.Geometry rings
// INPUT: 
//    mesh: full tessellation
//    low :  (x,y) coordinate defining the origin of the quantized grid
//    dx  :  Quantized grid spacing
// PRE-CONDITIONS: 
//    The mesh cells, faces, and nodes are full
//------------------------------------------------------------------------------
template <typename RealType>
void convertTessellationToRings(const polytope::Tessellation<2,RealType>& mesh,
                                const RealType* low,
                                const RealType dx,
                                std::vector<boost::geometry::model::ring<
                                   Point2<int64_t>, false> >& cellRings) {
  // Pre-conditions
  POLY_ASSERT(!mesh.cells.empty());
  POLY_ASSERT(!mesh.faces.empty());
  POLY_ASSERT(!mesh.nodes.empty());
  POLY_ASSERT(low != 0 and dx != 0);

  // typedefs
  typedef int64_t CoordHash;
  typedef Point2<CoordHash> IntPoint;
  typedef boost::geometry::model::ring<IntPoint, false> BGring;

  const int numCells = mesh.cells.size();
  cellRings.resize(numCells);
  for (unsigned i = 0; i != numCells; ++i) {
    std::vector<IntPoint> cellNodes;
    for (std::vector<int>::const_iterator faceItr = mesh.cells[i].begin();
         faceItr != mesh.cells[i].end(); ++faceItr){
      const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
      POLY_ASSERT(iface < mesh.faceCells.size());
      POLY_ASSERT(mesh.faces[iface].size() == 2);
      const unsigned inode = *faceItr < 0 ? mesh.faces[iface][1] : mesh.faces[iface][0];
      cellNodes.push_back(IntPoint(mesh.nodes[2*inode  ],
                                   mesh.nodes[2*inode+1],
                                   low[0], low[1], dx));
    }
    boost::geometry::assign(cellRings[i], BGring(cellNodes.begin(), cellNodes.end()));
    boost::geometry::correct(cellRings[i]);
    POLY_ASSERT(cellRings[i].front() == cellRings[i].back());
    POLY_ASSERT(!boost::geometry::intersects(cellRings[i]));
  }
  POLY_ASSERT(cellRings.size() == mesh.cells.size());
}


//------------------------------------------------------------------------------
// intersectBoundingBox
//
// Compute the intersection of a line segment and a PLC. Returns number of
// intersections, intersected facet(s) along PLC, and intersection point(s)
// INPUT: 
//    point1, point2 : Line segment endpoints
//    numVertices    : Number of vertices of bounding box
//    vertices       : List of vertex coordinates
//    facets         : Facets of the bounding box
// OUTPUT:
//    facetIntersections : Vector of facet indices where intersection(s) occur
//    result             : Vector of coordinates for intersection(s)
//------------------------------------------------------------------------------
template<typename RealType>
inline
unsigned intersectBoundingBox(const RealType* point1,
			      const RealType* point2,
			      const unsigned numVertices,
			      const RealType* vertices,
			      const std::vector<std::vector<int> >& facets,
                              const RealType tol,
			      std::vector<int>& facetIntersections,
			      std::vector<RealType>& result) {
  unsigned i, j, numIntersections=0;
  RealType intersectionPoint[2];
  bool intersects, addPoint;
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
    addPoint = false;
    if (intersects) {
      RealType p[2];
      p[(ifacet+1)%2] = vertices[2*ifacet + (ifacet+1)%2];
      p[ ifacet   %2] = intersectionPoint[ifacet%2];

      // Make sure we're not adding the same point twice (i.e. with a corner)
      if (numIntersections == 0) addPoint = true;
      else {
         if (std::abs(p[0] - (result.back()-1)) > tol and 
             std::abs(p[1] - (result.back()  )) > tol ) addPoint = true;
      }
      
      if (addPoint) {
        ++numIntersections;
        result.push_back(p[0]);
        result.push_back(p[1]);
        facetIntersections.push_back(ifacet);
      }
    }
  }
  POLY_ASSERT(numIntersections == result.size()/2);
  return numIntersections;
}


//------------------------------------------------------------------------------
// computeCellNodesCollinear
//
// Compute the quantized node positions for each cell from a collection of
// collinear generators
// INPUT: 
//    points     : Vector of collinear generators
//    coords     : Quantized coordinate system object
// OUTPUT:
//    nodeMap    : Map from quantized node position to node index
//    cellNodes  : Collection of sorted node indices around each cell
//------------------------------------------------------------------------------
template<typename RealType>
void
computeCellNodesCollinear(const std::vector<RealType>& points,
                          const QuantizedCoordinates<2, RealType>& coords,
			  std::vector<Point2<RealType> >& nodeList,
			  std::vector<std::vector<unsigned> >& cellNodes) {  
  const unsigned numGenerators = points.size()/2;
  int i;

  // typedefs
  typedef Point2<RealType> RealPoint;

  // Sort the generators but keep their original indices
  std::vector<std::pair<RealPoint,int> > pointIndexPairs;
  for (i = 0; i != numGenerators; ++i){
    pointIndexPairs.push_back(make_pair(RealPoint(points[2*i], points[2*i+1]), i));
  }
  sort( pointIndexPairs.begin(), pointIndexPairs.end(),
	internal::pairCompareFirst<RealPoint,int> );

  // Number of nodes
  const int nnodes = 2*numGenerators;

  unsigned inode, icell1, icell2;
  RealPoint p1, p2, r1, r2, node, midpt;
  nodeList.resize(nnodes);
  cellNodes.resize(numGenerators);
  
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
  node = coords.projectPoint(&p1.x, &r1.x);
  nodeList[inode] = node;

  // Node 1: endpt of first interior face
  node = coords.projectPoint(&midpt.x, &r2.x);
  nodeList[inode+1] = node;
  
  // Node 2: other endpt of first interior face
  r2 *= -1.0;
  node = coords.projectPoint(&midpt.x, &r2.x);
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
    node = coords.projectPoint(&midpt.x, &r2.x);
    nodeList[inode] = node;
    
    // Node 1: other endpt of interior face
    r2 *= -1.0;
    node = coords.projectPoint(&midpt.x, &r2.x);
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
  
  node = coords.projectPoint(&p2.x, &r1.x);
  nodeList[inode] = node;
    
  // Last node for final cell
  cellNodes[icell1].push_back(inode);

  // Post-conditions
  POLY_ASSERT(nodeList.size()  == nnodes       );
  POLY_ASSERT(cellNodes.size() == numGenerators);
}


} //end namespace polytope

#endif
