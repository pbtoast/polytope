#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>

#include "polytope.hh"
#include "polytope_test_utilities.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

// #include "Tessellator.hh"
// #include "Point.hh"

// // The Voronoi tools in Boost.Polygon
// #include <boost/polygon/voronoi.hpp>

// // Some Boost.Polygon stuff to make life easier
// #include <boost/geometry.hpp>
// #include <boost/geometry/geometries/geometries.hpp>
// #include <boost/geometry/geometries/register/point.hpp>

using namespace std;
using namespace polytope;

// //------------------------------------------------------------------------
// // Typedefs
// //------------------------------------------------------------------------
// typedef int64_t CoordHash;
// typedef Point2<CoordHash> IntPoint;
// typedef Point2<double> RealPoint;
// typedef boost::polygon::voronoi_diagram<double> VD;

// //------------------------------------------------------------------------
// // Mapping point stuff to Boost.Polygon
// //------------------------------------------------------------------------
// namespace boost{
// namespace polygon{
 
// template <>
// struct geometry_concept<IntPoint> { typedef point_concept type; };
  
// template <>
// struct point_traits<IntPoint> {
//   typedef CoordHash coordinate_type;
   
//   static inline coordinate_type get(const IntPoint& point, orientation_2d orient) {
//     return (orient == HORIZONTAL) ? point.x : point.y;
//   }
// };

// } // end boost
// } // end polygon


// //------------------------------------------------------------------------
// // Comparison operator for voronoi vertices
// //------------------------------------------------------------------------
// struct vertexCompare {
//    bool operator()(const VD::vertex_type v1, const VD::vertex_type v2) {
//       return (v1.x() < v2.x()                      ? true :
//               v1.x() == v2.x() and v1.y() < v2.y() ? true :
//               false);
//    }
// };


// //------------------------------------------------------------------------
// // computeInfiniteEdgeDirection
// //------------------------------------------------------------------------
// template<typename RealType>
// void
// computeInfiniteEdgeDirection(const VD::edge_type* edge,
//                              vector<IntPoint>& points,
// 			     RealType* endpt,
// 			     RealType* low,
// 			     RealType delta,
//                              RealType* direction) {
//   POLY_ASSERT(edge->is_infinite());
//   const VD::cell_type* cell1 = edge->cell();
//   const VD::cell_type* cell2 = edge->twin()->cell();
//   // Assume only point-generators for the time being
//   POLY_ASSERT(cell1->contains_point() and cell2->contains_point());
//   size_t index1 = cell1->source_index();
//   size_t index2 = cell2->source_index();
//   POLY_ASSERT(index1 < points.size() and index2 < points.size());
//   const IntPoint p1 = points[index1];
//   const IntPoint p2 = points[index2];
//   RealType midpt[2] = {0.5*(p1.realx(low[0],delta) + p2.realx(low[0],delta)), 
// 		       0.5*(p1.realy(low[1],delta) + p2.realy(low[1],delta))};
//   direction[0] = midpt[0] - endpt[0];
//   direction[1] = midpt[1] - endpt[1];
//   geometry::unitVector<2, double>(direction);
// }



//------------------------------------------------------------------------
// main
//------------------------------------------------------------------------
int 
main(int argc, char** argv) 
{
#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  const CoordHash coordMax = (1LL << 26);
  const double degeneracy = 1.5e-8;

  // Create the generators.
  vector<double> points;
  const int nx = 3;
  const double x1 = 0.0, y1 = 0.0;
  const double x2 = 1.0, y2 = 1.0;
  const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx;
  double low [2] = { numeric_limits<double>::max(),  numeric_limits<double>::max()};
  double high[2] = {-numeric_limits<double>::max(), -numeric_limits<double>::max()};
  unsigned ix, iy;
  double xi, yi;
  for (iy = 0; iy != nx; ++iy) {
    yi = y1 + (iy + 0.5)*dy;
    for (ix = 0; ix != nx; ++ix) {
      xi = x1 + (ix + 0.5)*dx;
      points.push_back(xi);
      points.push_back(yi);
      low [0] = min(low [0], xi);
      low [1] = min(low [1], yi);
      high[0] = max(high[0], xi);
      high[1] = max(high[1], yi);
    }
  }

  vector<double> PLCpoints;
  PLCpoints.push_back(x1);  PLCpoints.push_back(y1);
  PLCpoints.push_back(x2);  PLCpoints.push_back(y1);
  PLCpoints.push_back(x2);  PLCpoints.push_back(y2);
  PLCpoints.push_back(x1);  PLCpoints.push_back(y2);

  PLC<2,double> boundary;
  int nSides = 4;
  boundary.facets.resize(nSides, vector<int>(2));
  for (int i = 0; i != nSides; ++i) {
    boundary.facets[i][0] = i;
    boundary.facets[i][1] = (i+1)%nSides;
  }

  Tessellation<2,double> tmesh;
  BoostTessellator<double> tessellator;
  //tessellator.tessellate(points, tmesh);
  tessellator.tessellate(points, PLCpoints, boundary, tmesh);
  cout << tmesh << endl;

#if HAVE_SILO
  vector<double> index(tmesh.cells.size());
  for (int i = 0; i < tmesh.cells.size(); ++i) index[i] = double(i);
  map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
  cellFields["cell_index"] = &index[0];
  ostringstream os;
  os << "test_BoostTessellator_mesh";
  polytope::SiloWriter<2, double>::write(tmesh, nodeFields, edgeFields, 
                                         faceFields, cellFields, os.str());
#endif

//   const int numGenerators = points.size()/2;

//   const double box[2] = {high[0] - low[0], high[1] - low[1]};
//   const double rinf   = 4.0*max(box[0], box[1]);
//   const double cen[2] = {0.5*(low[0] + high[0]), 0.5*(low[1] + high[1])};

//   low [0] = min(low [0], cen[0]-rinf);
//   high[0] = max(high[0], cen[0]+rinf);
//   low [1] = min(low [1], cen[1]-rinf);
//   high[1] = max(high[1], cen[1]+rinf);

//   const double delta = max(degeneracy, 2.0*rinf/coordMax);

//   vector<IntPoint> generators(numGenerators);
//   for (int i = 0; i != numGenerators; ++i) {
//     generators[i] = IntPoint(points[2*i], points[2*i+1],
// 			     low[0], low[1], delta);
//   }

//   sort(generators.begin(), generators.end());
  
//   VD vd;
//   construct_voronoi(generators.begin(), generators.end(), &vd);
//   Tessellation<2, double> mesh;
  
//   // Set the "color" of each edge and use it as a local index
//   bool test;
//   size_t colorIndex = 1;
//   int k, nodeIndex, faceIndex, cellIndex=0;
//   RealPoint direction, pinf;
//   IntPoint node;
//   map<VD::vertex_type, int, vertexCompare> vertexToNodeIndex;
//   set<int> infCells;
//   mesh.cells.resize(vd.num_cells());
//   for (VD::const_cell_iterator cellItr = vd.cells().begin(); 
//        cellItr != vd.cells().end(); ++cellItr, ++cellIndex) {
//     const VD::edge_type* edge = cellItr->incident_edge();
//     do {
//       POLY_ASSERT2(colorIndex < numeric_limits<size_t>::max()/32, "BoostTessellator Error: "
//                    << "overflow of maximum allowable face index");

//       if (edge->is_infinite()) {
//         set<int>::iterator it = infCells.find(cellIndex);
//         if (it == infCells.end())  infCells.insert(cellIndex);
//       }

//       // // Blago!
//       // cout << endl << "Edge index " << colorIndex << endl;
//       // cout << "   Infinite edge? " << ((edge->is_infinite()) ? "yes" : "no") << endl;
//       // cout << "   Is it primary? " << ((edge->is_primary()) ? "yes" : "no" ) << endl;
//       // cout << "   Vertex 0:" << endl;
//       // if (edge->vertex0()) {
//       //    cout << "      position = (" << edge->vertex0()->x() << ","
//       //         << edge->vertex0()->y() << ")" << endl;
//       // } else {
//       //    cout << "      Inf node" << endl;
//       // }
//       // cout << "   Vertex 1:" << endl;
//       // if (edge->vertex1()) {
//       //    cout << "      position = (" << edge->vertex1()->x() << ","
//       //         << edge->vertex1()->y() << ")" << endl;
//       // } else {
//       //    cout << "      Inf node" << endl;
//       // }
//       // // Blago!

//       if (edge->twin()->color() == 0) {
//         faceIndex = colorIndex-1;
//         edge->color(colorIndex);
//         ++colorIndex;
        
//         // Input cell-to-face info into mesh.cells
//         POLY_ASSERT(cellIndex < vd.num_cells());
//         POLY_ASSERT(cellItr->source_index() == cellIndex);
//         mesh.cells[cellIndex].push_back(faceIndex);
//         size_t oppCellIndex = edge->twin()->cell()->source_index();
//         POLY_ASSERT(oppCellIndex < vd.num_cells());
//         mesh.cells[oppCellIndex].push_back(~faceIndex);
        
//         // Input face-to-cell info info mesh.faceCells
//         mesh.faceCells.push_back(vector<int>(2));
//         POLY_ASSERT(faceIndex < mesh.faceCells.size());
//         POLY_ASSERT(mesh.faceCells[faceIndex].size() == 2);
//         mesh.faceCells[faceIndex][0] = cellIndex;
//         mesh.faceCells[faceIndex][1] = ~(int(oppCellIndex));
        
//         // Size the face and infFace arrays
//         mesh.faces.push_back(vector<unsigned>(2));
//         mesh.infFaces.push_back(0);
//         POLY_ASSERT(faceIndex < mesh.faces.size());
//         POLY_ASSERT(mesh.faces[faceIndex].size() == 2);

//         // The two vertex pointers for this edge
//         // NOTE: If edge is infinite, one of them will be null
//         const VD::vertex_type* v0 = edge->vertex0();
//         const VD::vertex_type* v1 = edge->vertex1();

//         // Finite edge: Add both vertices to the map
//         if (edge->is_finite()) {
//           k = vertexToNodeIndex.size();
//           nodeIndex = internal::addKeyToMap(*v0, vertexToNodeIndex);
//           if (nodeIndex == k) {
// 	    node = IntPoint(v0->x(), v0->y());
//             mesh.nodes.push_back(node.realx(low[0], delta));
//             mesh.nodes.push_back(node.realy(low[1], delta));
//             mesh.infNodes.push_back(0);
//           }
//           mesh.faces[faceIndex][0] = nodeIndex;
          
//           k = vertexToNodeIndex.size();
//           nodeIndex = internal::addKeyToMap(*v1, vertexToNodeIndex);
//           if (nodeIndex == k) {
// 	    node = IntPoint(v1->x(), v1->y());
//             mesh.nodes.push_back(node.realx(low[0], delta));
//             mesh.nodes.push_back(node.realy(low[1], delta));
//             mesh.infNodes.push_back(0);
//           }
//           mesh.faces[faceIndex][1] = nodeIndex;
//         }

//         // Infinite edge: Determine the direction of the ray pointing to infinity.
//         // Add the origin vertex of the ray and the projected point
//         else {
//           // Project the non-null vertex to infinity along the edge direction
//           //computeInfiniteEdgeDirection(edge, generators, &direction.x);
          
// 	  IntPoint finiteVertex = (edge->vertex0())              ?
// 	    IntPoint(edge->vertex0()->x(), edge->vertex0()->y()) :
// 	    IntPoint(edge->vertex1()->x(), edge->vertex1()->y());
// 	  RealPoint endpt = RealPoint(finiteVertex.realx(low[0],delta),
// 				      finiteVertex.realy(low[1],delta));
// 	  computeInfiniteEdgeDirection(edge, generators, &endpt.x, low, delta, &direction.x);
	  
// // 	  // Blago!
// // 	  const VD::cell_type* cell1 = edge->cell();
// // 	  const VD::cell_type* cell2 = edge->twin()->cell();
// // 	  // Assume only point-generators for the time being
// // 	  POLY_ASSERT(cell1->contains_point() and cell2->contains_point());
// // 	  size_t index1 = cell1->source_index();
// // 	  size_t index2 = cell2->source_index();
// // 	  POLY_ASSERT(index1 < points.size() and index2 < points.size());
// // 	  const IntPoint p1 = generators[index1];
// // 	  const IntPoint p2 = generators[index2];
// // 	  RealPoint midpt = RealPoint(0.5*(p1.realx(low[0],delta) + p2.realx(low[0],delta)), 
// // 				      0.5*(p1.realy(low[1],delta) + p2.realy(low[1],delta)));
// // 	  IntPoint vertex = (edge->vertex0())                    ?
// // 	    IntPoint(edge->vertex0()->x(), edge->vertex0()->y()) :
// // 	    IntPoint(edge->vertex1()->x(), edge->vertex1()->y());
// // 	  RealPoint endpt = RealPoint(vertex.realx(low[0],delta),
// // 				      vertex.realy(low[1],delta));
// // 	  direction.x = midpt.x - endpt.x;
// // 	  direction.y = midpt.y - endpt.y;
// // 	  geometry::unitVector<2, double>(&direction.x);
// // 	  // Blago!
	  
// 	  //RealPoint node = v0           ? 
// 	  //  RealPoint(v0->x(), v0->y()) :
// 	  //  RealPoint(v1->x(), v1->y());
// 	  //test = geometry::rayCircleIntersection(&node.x, &direction.x,
//           //                                       cen, rinf, 1.0e-6, &pinf.x);
//           test = geometry::rayCircleIntersection(&endpt.x, &direction.x,
//                                                  cen, rinf, 1.0e-6, &pinf.x);
//           POLY_ASSERT(test);

//           // Vertex 0 exists, vertex 1 is a projected infNode. Add them in order
//           if (v0) {
//             k = vertexToNodeIndex.size();
//             nodeIndex = internal::addKeyToMap(*v0, vertexToNodeIndex);
//             if (nodeIndex == k) {
//               node = IntPoint(v1->x(), v1->y());
//               mesh.nodes.push_back(node.realx(low[0], delta));
//               mesh.nodes.push_back(node.realy(low[1], delta));
// 	      mesh.infNodes.push_back(0);
//             }
//             mesh.faces[faceIndex][0] = nodeIndex;

//             k = vertexToNodeIndex.size();
//             nodeIndex = internal::addKeyToMap(VD::vertex_type(pinf.x,pinf.y),
//                                               vertexToNodeIndex);
//             if (nodeIndex == k) {
//               mesh.nodes.push_back(pinf.x);
//               mesh.nodes.push_back(pinf.y);
//               mesh.infNodes.push_back(1);
//             }
//             mesh.faces[faceIndex][1] = nodeIndex;
//           } 
          
//           // Vertex 0 is a projected infNode, vertex 1 exists. Add them in order
//           else {
//             k = vertexToNodeIndex.size();
//             nodeIndex = internal::addKeyToMap(VD::vertex_type(pinf.x,pinf.y),
//                                               vertexToNodeIndex);
//             if (nodeIndex == k) {
//               mesh.nodes.push_back(pinf.x);
//               mesh.nodes.push_back(pinf.y);
//               mesh.infNodes.push_back(1);
//             }
//             mesh.faces[faceIndex][0] = nodeIndex;
            
// 	    k = vertexToNodeIndex.size();
//             nodeIndex = internal::addKeyToMap(*v1, vertexToNodeIndex);
//             if (nodeIndex == k) {
//               node = IntPoint(v1->x(), v1->y());
//               mesh.nodes.push_back(node.realx(low[0], delta));
//               mesh.nodes.push_back(node.realy(low[1], delta));
// 	      mesh.infNodes.push_back(0);
//             }
//             mesh.faces[faceIndex][1] = nodeIndex;
//           }
//         }
//       }
//       edge = edge->next();
//     } while (edge != cellItr->incident_edge());
//   }
//   int numFaces = colorIndex-1;
//   // A few sanity checks
//   POLY_ASSERT(mesh.faces.size()     == numFaces);
//   POLY_ASSERT(mesh.faceCells.size() == numFaces);
//   POLY_ASSERT(mesh.infNodes.size()  == mesh.nodes.size()/2);
//   POLY_ASSERT(mesh.infFaces.size()  == numFaces);

//   // If the next edge in the list is also infinite, we add an inf face to
//   // complete the faces around a cell (even in the unbounded case)
//   for (set<int>::const_iterator it = infCells.begin(); 
//        it != infCells.end(); ++it){
//      int icell = *it;
//      POLY_ASSERT(mesh.cells[icell].size() > 1);
//      vector<int>::iterator faceItr0 = mesh.cells[icell].end()-1;
//      vector<int>::iterator faceItr1 = mesh.cells[icell].begin();
//      bool lookingForInfFace = true;
//      unsigned infNode0, infNode1;
//      while (lookingForInfFace and faceItr1 != mesh.cells[icell].end()) {
//         const unsigned iface0 = *faceItr0 < 0 ? ~(*faceItr0) : *faceItr0;
//         const unsigned iface1 = *faceItr1 < 0 ? ~(*faceItr1) : *faceItr1;
//         POLY_ASSERT(iface0 < mesh.faces.size() and iface1 < mesh.faces.size());
//         POLY_ASSERT(mesh.faces[iface0].size() == 2 and mesh.faces[iface1].size() == 2);
//         const unsigned inode00 = *faceItr0 < 0 ? mesh.faces[iface0][1] : mesh.faces[iface0][0];
//         const unsigned inode01 = *faceItr0 < 0 ? mesh.faces[iface0][0] : mesh.faces[iface0][1];
//         const unsigned inode10 = *faceItr1 < 0 ? mesh.faces[iface1][1] : mesh.faces[iface1][0];
//         const unsigned inode11 = *faceItr1 < 0 ? mesh.faces[iface1][0] : mesh.faces[iface1][1];
//         POLY_ASSERT(inode00 < mesh.nodes.size()/2 and inode01 < mesh.nodes.size()/2 and
//                     inode10 < mesh.nodes.size()/2 and inode11 < mesh.nodes.size()/2 );
//         if (mesh.infNodes[inode00] == 0 and 
//             mesh.infNodes[inode01] == 1 and
//             mesh.infNodes[inode10] == 1 and
//             mesh.infNodes[inode11] == 0) {
//           POLY_ASSERT(inode01 != inode10);
//           lookingForInfFace = false;
//           infNode0 = inode01;
//           infNode1 = inode10;
//         }
//         faceItr0 = faceItr1;
//         ++faceItr1;
//      }
//      POLY_ASSERT(!lookingForInfFace);
//      vector<unsigned> faceNodes(2);
//      faceNodes[0] = infNode0;
//      faceNodes[1] = infNode1;
//      mesh.faces.push_back(faceNodes);
//      mesh.faceCells.push_back(vector<int>(1, icell));
//      mesh.infFaces.push_back(1);
//      mesh.cells[icell].insert(faceItr0, numFaces);
//      ++numFaces;
//   }
//   // A few repeat sanity checks
//   POLY_ASSERT(mesh.faces.size()     == numFaces);
//   POLY_ASSERT(mesh.faceCells.size() == numFaces);
//   POLY_ASSERT(mesh.infFaces.size()  == numFaces);

//   cout << mesh << endl;
 
// #if HAVE_SILO
//   vector<double> index(mesh.cells.size());
//   for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
//   map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
//   cellFields["cell_index"] = &index[0];
//   ostringstream os;
//   os << "test_BoostTessellator_mesh";
//   polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
//                                          faceFields, cellFields, os.str());
// #endif


#if HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
//------------------------------------------------------------------------
