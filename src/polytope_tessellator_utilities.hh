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
#if HAVE_BOOST
BOOST_GEOMETRY_REGISTER_POINT_2D(polytope::Point2<int64_t>, int64_t, boost::geometry::cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_POINT_2D(polytope::Point2<double>, double, boost::geometry::cs::cartesian, x, y);
#endif


// Fast predicate for determining point orientation
//extern double orient2d(double* pa, double* pb, double* pc);


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
constructUnboundedMeshTopology(std::map<int, std::vector<unsigned> >& cellNodes,
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
//    geometry:  PLC boundary
//    low     :  (x,y) coordinate defining the origin of the quantized grid
//    dx      :  Quantized grid spacing
// PRE-CONDITIONS: 
//    The mesh is empty
//------------------------------------------------------------------------------
template<typename RealType>
void
constructBoundedMeshTopology(const std::vector<boost::geometry::model::ring
                               <polytope::Point2<int64_t>, false> >& cellRings,
                             const std::vector<RealType>& points,
                             const std::vector<RealType>& PLCpoints,
                             const PLC<2,RealType>& geometry,
                             const RealType* low,
                             const RealType dx,
                             Tessellation<2,RealType>& mesh) {
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(low != 0 and dx != 0);
  
  // Typedefs
  typedef int64_t CoordHash;
  typedef std::pair<int, int> EdgeHash;
  typedef polytope::Point2<CoordHash> IntPoint;
  typedef polytope::Point2<double> RealPoint;
  typedef boost::geometry::model::ring<IntPoint, false> BGring;

  const unsigned numGenerators = cellRings.size();
  const unsigned numPLCpoints  = PLCpoints.size();

  // Now build the unique mesh nodes and cell info.
  std::map<IntPoint, int> point2node;
  std::map<EdgeHash, int> edgeHash2id;
  std::map<int, std::vector<int> > edgeCells;
  int i, j, k, iedge;
  mesh.cells = std::vector<std::vector<int> >(numGenerators);
  for (i = 0; i != numGenerators; ++i) { 
    POLY_ASSERT(cellRings[i].size() > 2);
    //boost::geometry::unique(cellRings[i]);
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
  mesh.nodes = std::vector<RealType>(2*point2node.size());
  for (typename std::map<IntPoint, int>::const_iterator itr = point2node.begin();
       itr != point2node.end(); ++itr) {
    const IntPoint& p = itr->first;
    i = itr->second;
    POLY_ASSERT(i < mesh.nodes.size()/2);
    node[0] = p.realx(low[0],dx);
    node[1] = p.realy(low[1],dx);
    
    // Check if nodes are inside boundary (either bounding box or PLC, if defined)
    bool inside = within(node, numPLCpoints, &PLCpoints[0], geometry);
    if(!inside){
      RealType result[2];
      RealType dist = nearestPoint( node, numPLCpoints, &PLCpoints[0], geometry, result );
      // Check the node has not moved more than 2.5 quantized mesh spacings. NOTE: this is
      // not a sharp estimate. Theoreticallly, the distance ought to be at most sqrt(2)*dx, 
      // but nodes will fail this strict of a test.
      POLY_ASSERT2( dist < 5.0*dx,
                    dist << " " << 2.5*dx << " : "
                    << "(" << node[0]   << " " << node[1]   << ") "
                    << "(" << result[0] << " " << result[1] << ")\n" << geometry );
      node[0] = result[0];
      node[1] = result[1];
    }
    mesh.nodes[2*i]   = node[0];
    mesh.nodes[2*i+1] = node[1];
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
    // Make sure we're not adding the same point twice (i.e. with a corner)
    addPoint = false;
    if (intersects) {
      if (numIntersections == 0) addPoint = true;
      else {
        if (intersectionPoint[0] != result.back() - 1 and
            intersectionPoint[1] != result.back() ) addPoint = true;
        else addPoint = false;
      }
    }

    if (addPoint) {
      RealType p[2];
      p[(ifacet+1)%2] = vertices[2*ifacet + (ifacet+1)%2];
      p[ ifacet   %2] = intersectionPoint[ifacet%2];
      ++numIntersections;
      result.push_back(p[0]);
      result.push_back(p[1]);
      facetIntersections.push_back(ifacet);
    }
  }
  POLY_ASSERT(numIntersections == result.size()/2);
  return numIntersections;
}



 


} //end namespace polytope

#endif
