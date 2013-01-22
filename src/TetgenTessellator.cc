//------------------------------------------------------------------------
// TetgenTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>

#include "polytope.hh" // Pulls in POLY_ASSERT and TetgenTessellator.hh.
#include "Point.hh"

// Pull in tetgen stuff.
#define TETLIBRARY
#include "tetgen.h"

namespace polytope {

using namespace std;
using std::min;
using std::max;
using std::abs;

namespace {

//------------------------------------------------------------------------
// This function computes the determinant of the 3x3 matrix A with 
// the given 9 coefficients.
//------------------------------------------------------------------------
double 
det3(double A11, double A12, double A13,
     double A21, double A22, double A23,
     double A31, double A32, double A33) {
  return A11*(A22*A33-A32*A23) - A12*(A21*A33-A31*A23) + A13*(A21*A31-A31*A22);
}

//------------------------------------------------------------------------
// This function computes the determinant of the 4x4 matrix A with 
// the given 16 coefficients.
//------------------------------------------------------------------------
double 
det4(double A11, double A12, double A13, double A14,
     double A21, double A22, double A23, double A24,
     double A31, double A32, double A33, double A34,
     double A41, double A42, double A43, double A44) {
  double det31 = det3(A22, A23, A24,
                      A32, A33, A34,
                      A42, A43, A44);
  double det32 = det3(A21, A23, A24,
                      A31, A33, A34,
                      A41, A43, A44);
  double det33 = det3(A21, A22, A24,
                      A31, A32, A34,
                      A41, A42, A44);
  double det34 = det3(A21, A22, A23,
                      A31, A32, A33,
                      A41, A42, A43);
  return A11*det31 - A12*det32 + A13*det33 - A14*det34;
}

//------------------------------------------------------------------------
// This function computes the circumcenter of a tetrahedron with vertices
// A = (Ax, Ay, Az), B = (Bx, By, Bz), C = (Cx, Cy, Cz), and D = (Dx, Dy, Dz),
// and places the result in X.
//------------------------------------------------------------------------
void 
computeCircumcenter(double* A, double* B, double* C, double* D, double* X) {
  // This solution was taken from 
  // http://mathworld.wolfram.com/Circumsphere.html.
  double r12 = A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
  double r22 = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
  double r32 = C[0]*C[0] + C[1]*C[1] + C[2]*C[2];
  double r42 = D[0]*D[0] + D[1]*D[1] + D[2]*D[2];
  double a = det4(A[0], A[1], A[2], 1.0,
                  B[0], B[1], B[2], 1.0,
                  C[0], C[1], C[2], 1.0,
                  D[0], D[1], D[2], 1.0);
  double Dx = det4(r12, A[1], A[2], 1.0,
                   r22, B[1], B[2], 1.0,
                   r32, C[1], C[2], 1.0,
                   r42, D[1], D[2], 1.0);
  double Dy = det4(r12, A[0], A[2], 1.0,
                   r22, B[0], B[2], 1.0,
                   r32, C[0], C[2], 1.0,
                   r42, D[0], D[2], 1.0);
  double Dz = det4(r12, A[0], A[1], 1.0,
                   r22, B[0], B[1], 1.0,
                   r32, C[0], C[1], 1.0,
                   r42, D[0], D[1], 1.0);
  X[0] = 0.5*Dx/a;
  X[1] = 0.5*Dy/a;
  X[2] = 0.5*Dz/a;
}

//------------------------------------------------------------------------------
// An implementation of the map specialized to help constructing counters.
// This thing just overloads the index operator to start the count at zero
// for new key values.
//------------------------------------------------------------------------------
template<typename Key, 
         typename Comparator = std::less<Key> >
class CounterMap: public std::map<Key, unsigned> {
public:
  CounterMap(): std::map<Key, unsigned>() {}
  virtual ~CounterMap() {}
  unsigned& operator[](const Key& key) {
    typename std::map<Key, unsigned>::iterator itr = this->find(key);
    if (itr == this->end()) {
      std::map<Key, unsigned>::operator[](key) = 0U;
      itr = this->find(key);
    }
    POLY_ASSERT(itr != this->end());
    return itr->second;
  }
};

//------------------------------------------------------------------------------
// Increment the given position by 
//------------------------------------------------------------------------------
void
incrementPosition(Point3<double>& val,
                  const vector<double>& points,
                  const vector<int>& indices) {
  Point3<double> xi(0.0, 0.0, 0.0);
  for (unsigned k = 0; k != indices.size(); ++k) {
    const unsigned i = indices[k];
    POLY_ASSERT(i < points.size()/3);
    xi.x += points[3*i];
    xi.y += points[3*i + 1];
    xi.z += points[3*i + 2];
  }
  xi /= indices.size();
  val += xi;
}

} // end anonymous namespace

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
TetgenTessellator::
TetgenTessellator(const bool directComputation):
  Tessellator<3, double>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
TetgenTessellator::
~TetgenTessellator() {
}

//------------------------------------------------------------------------------
// Unbounded tessellation.
//------------------------------------------------------------------------------
void
TetgenTessellator::
tessellate(const vector<double>& points,
           Tessellation<3, double>& mesh) const {

  // Pre-conditions.
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 3 == 0);

  typedef int64_t CoordHash;
  typedef pair<int, int> EdgeHash;
  typedef set<int> FaceHash;
  typedef Point3<CoordHash> IntPoint;
  typedef Point3<RealType> RealPoint;
  const CoordHash coordMax = (1LL << 34); // numeric_limits<CoordHash>::max() >> 32U;
  const double degeneracy = 1.0e-12;

  // Compute the normalized generators.
  RealType low[3], high[3], box[3];
  const unsigned numGenerators = points.size() / 3;
  vector<double> generators = this->computeNormalizedPoints(points, points, true, low, high);
  unsigned i, j, k;
  for (j = 0; j != 3; ++j) {
    box[j] = high[j] - low[j];
    POLY_ASSERT(box[j] > 0.0);
  }

  // Build the input to tetgen.
  tetgenio in;
  in.firstnumber = 0;
  in.mesh_dim = 3;
  in.pointlist = &generators.front();
  in.pointattributelist = 0;
  in.pointmtrlist = 0;
  in.pointmarkerlist = 0;
  in.numberofpoints = generators.size() / 3;
  in.numberofpointattributes = 0;
  in.numberofpointmtrs = 0;

  // Do the tetrahedralization.
  tetgenio out;
  tetrahedralize((char*)"vV", &in, &out);

  // Make sure we got something.
  if (out.numberoftetrahedra == 0)
    error("TetgenTessellator: Delauney tetrahedralization produced 0 tetrahedra!");
  if (out.numberofpoints != generators.size()/3) {
    char err[1024];
    snprintf(err, 1024, "TetgenTessellator: Delauney tetrahedralization produced %d tetrahedra\n(%d generating points given)", 
             out.numberofpoints, (int)numGenerators);
    error(err);
  }

  cerr << "Finished tetgen." << endl;

  // Find the maximum extent of the Voronoi points.  We will use this to create the
  // unbounded (or "infinite") points.
  RealType vlow[3], vhigh[3];
  geometry::computeBoundingBox<3, RealType>(out.vpointlist,
                                            3*out.numberofvpoints,
                                            true,
                                            vlow,
                                            vhigh);
  const RealType rinf = 2.0*max(abs(vlow[0]), 
                                max(abs(vlow[1]), 
                                    max(abs(vlow[2]),
                                        max(vhigh[0],
                                            max(vhigh[1], vhigh[2])))));

  // Read out the vertices, converting to quantized (hashed) positions and
  // building the map of hash -> node ID.
  map<IntPoint, int> point2node;
  for (i = 0; i != out.numberofvpoints; ++i) {
    IntPoint p(out.vpointlist[3*i], out.vpointlist[3*i + 1], out.vpointlist[3*i + 2], 0.0, 0.0, 0.0, degeneracy);
    internal::addKeyToMap(p, point2node);
  }
  const unsigned nvertices = point2node.size();

  // Read the normalized mesh node positions.
  mesh.nodes.resize(3*point2node.size());
  for (map<IntPoint, int>::const_iterator itr = point2node.begin();
       itr != point2node.end();
       ++itr) {
    i = itr->second;
    mesh.nodes[3*i  ] = itr->first.realx(0.0, degeneracy);
    mesh.nodes[3*i+1] = itr->first.realy(0.0, degeneracy);
    mesh.nodes[3*i+2] = itr->first.realz(0.0, degeneracy);
  }

  // Read the tetgen voro edges, creating our "infinte" surface nodes for
  // unbounded rays.
  int v1, v2;
  unsigned n1, n2;
  RealType origin[3] = {0.0, 0.0, 0.0}, ray_sph_int[3];
  vector<EdgeHash> edges;
  edges.reserve(out.numberofvedges);
  mesh.infNodes = vector<unsigned>();
  for (i = 0; i != out.numberofvedges; ++i) {
    tetgenio::voroedge& vedge = out.vedgelist[i];
    v1 = vedge.v1;
    v2 = vedge.v2;
    POLY_ASSERT2(v1 >= 0 and v1 < nvertices, i << " " << v1 << " " << v2);
    POLY_ASSERT2(v2 == -1 or (v2 >= 0 and v2 < nvertices), i << " " << v1 << " " << v2);
    IntPoint p1(out.vpointlist[3*v1], out.vpointlist[3*v1+1], out.vpointlist[3*v1+2], 0.0, 0.0, 0.0, degeneracy);
    POLY_ASSERT(point2node.find(p1) != point2node.end());
    n1 = point2node[p1];
    if (v2 == -1) {
      // This edge is a ray, so we construct a node on the rinf spherical surface
      // to hook to.
      geometry::rayPlaneIntersection(&mesh.nodes[3*n1],
                                     out.vedgelist[i].vnormal,
                                     origin,
                                     rinf,
                                     degeneracy,
                                     ray_sph_int);
      IntPoint p1(ray_sph_int[0], ray_sph_int[1], ray_sph_int[2], 0.0, 0.0, 0.0, degeneracy);
      k = point2node.size();
      point2node[p1] = k;
      mesh.infNodes.push_back(k);
      n2 = k;
    } else {
      IntPoint p2(out.vpointlist[3*v2], out.vpointlist[3*v2+1], out.vpointlist[3*v2+2], 0.0, 0.0, 0.0, degeneracy);
      POLY_ASSERT(point2node.find(p2) != point2node.end());
      n2 = point2node[p2];
    }
    POLY_ASSERT(n1 != n2);
    edges.push_back(internal::hashEdge(n1, n2));
  }
  POLY_ASSERT(edges.size() == out.numberofvedges);

  // Create the faces.
  const unsigned nfaces = out.numberofvfacets;
  mesh.faces.resize(nfaces);
  mesh.faceCells.resize(nfaces);
  int ne, ie, ie_last;
  unsigned m1, m2;
  for (i = 0; i != nfaces; ++i) {
    const tetgenio::vorofacet& vfacet = out.vfacetlist[i];
    ne = vfacet.elist[0];
    POLY_ASSERT(ne >= 3);
    mesh.faces[i].reserve(ne);
    for (k = 0; k != ne; ++k) {
      ie = vfacet.elist[k + 1];
      POLY_ASSERT(ie < edges.size());
      n1 = edges[ie].first;
      n2 = edges[ie].second;

      if (k == 0) {
        mesh.faces[i].push_back(n1);
        mesh.faces[i].push_back(n2);
      } else if (k < ne - 1) {
        if (n1 == mesh.faces[i].back()) {
          mesh.faces[i].push_back(n2);
        } else if (n2 == mesh.faces[i].back()) {
          mesh.faces[i].push_back(n1);
        } else {
          // If we couldn't find any continuity in the nodes from the last
          // edge, both this and the previous edge should be unbound rays.
          ie_last = (ie - 1) % ne;
          POLY_ASSERT(out.vedgelist[ie_last].v2 == -1 and
                      out.vedgelist[ie].v2 == -1);
          mesh.faces[i].push_back(n2);
          mesh.faces[i].push_back(n1);
        }
      } else {
        ie_last = (ie - 1) % ne;
        if (out.vedgelist[ie_last].v2 == -1 and
            out.vedgelist[ie].v2 == -1) {
          POLY_ASSERT(mesh.faces[i].front() == n1);
          mesh.faces[i].push_back(n2);
        }
      }
    }

    // Now the cells that touch this face.
    POLY_ASSERT(vfacet.c1 >= 0 and vfacet.c1 < numGenerators);
    POLY_ASSERT(vfacet.c2 >= 0 and vfacet.c2 < numGenerators);
    mesh.faceCells[i].push_back(vfacet.c1);
    mesh.faceCells[i].push_back(vfacet.c2);
  }
  POLY_ASSERT(mesh.faces.size() == nfaces);

  // Read out the cell structure as collections of faces.
  POLY_ASSERT(out.numberofvcells == numGenerators);
  mesh.cells = vector<vector<int> >(numGenerators);
  unsigned nf;
  int fi;
  RealType ccent[3], fcent[3], fhat[3];
  for (i = 0; i != numGenerators; ++i) {
    nf = out.vcelllist[i][0];
    POLY_ASSERT(nf >= 4);
    std::copy(&out.vcelllist[i][1], &out.vcelllist[i][1] + nf, std::back_inserter(mesh.cells[i]));
    POLY_ASSERT(mesh.cells[i].size() == nf);

    // Now we have to walk the face indices again and orient them with respect to this cell.
    geometry::computeCellCentroid(mesh, i, ccent);
    for (k = 0; k != nf; ++k) {
      fi = mesh.cells[i][k];
      POLY_ASSERT(fi >= 0 and fi < mesh.faces.size());

      // We have to determine the orientation of this face with respect to this cell.
      geometry::computeFaceCentroidAndNormal(mesh, fi, fcent, fhat);
      for (j = 0; j != 3; ++j) fcent[j] -= ccent[j];
      if (geometry::dot<3, RealType>(fcent, fhat) < 0.0) mesh.cells[i][k] = ~fi;
    }
  }

  // Rescale the mesh node positions for the input geometry.
  for (i = 0; i != mesh.nodes.size(); ++i) {
    j = i % 3;
    mesh.nodes[i] = low[j] + mesh.nodes[i]*box[j];
  }
}

//------------------------------------------------------------------------------
// Tessellate within a bounding box.
//------------------------------------------------------------------------------
void
TetgenTessellator::
tessellate(const vector<double>& points,
           double* low,
           double* high,
           Tessellation<3, double>& mesh) const {

  // Pre-conditions.
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 3 == 0);
  POLY_ASSERT(low[0] < high[0] and
              low[1] < high[1] and
              low[2] < high[2]);

  typedef int64_t CoordHash;
  typedef pair<int, int> EdgeHash;
  typedef set<int> FaceHash;
  typedef Point3<CoordHash> IntPoint;
  typedef Point3<RealType> RealPoint;
  const CoordHash coordMax = (1LL << 34); // numeric_limits<CoordHash>::max() >> 32U;
  const double degeneracy = 1.0e-12;

  // Compute the normalized generators.
  RealType box[3];
  const unsigned numGenerators = points.size() / 3;
  vector<double> generators = this->computeNormalizedPoints(points, points, false, low, high);
  for (unsigned j = 0; j != 3; ++j) {
    box[j] = high[j] - low[j];
    POLY_ASSERT(box[j] > 0.0);
  }

  // Build the input to tetgen.
  tetgenio in;
  in.firstnumber = 0;
  in.mesh_dim = 3;
  in.pointlist = &generators.front();
  in.pointattributelist = 0;
  in.pointmtrlist = 0;
  in.pointmarkerlist = 0;
  in.numberofpoints = generators.size() / 3;
  in.numberofpointattributes = 0;
  in.numberofpointmtrs = 0;

  // // Copy the PLC boundary info to the tetgen input.
  // vector<tetgenio::polygon> plcFacetPolygons(geometry.facets.size());
  // for (unsigned ifacet = 0; ifacet != geometry.facets.size(); ++ifacet) {
  //   POLY_ASSERT(geometry.facets[ifacet].size() >= 3);
  //   plcFacetPolygons[ifacet].numberofvertices = geometry.facets[ifacet].size();
  //   plcFacetPolygons[ifacet].vertexlist = const_cast<int*>(&geometry.facets[ifacet].front());
  // }
  // unsigned nfacets = plcFacetPolygons.size();
  // vector<RealPoint> holeCentroids(geometry.holes.size(), RealPoint(0.0, 0.0, 0.0));
  // for (unsigned ihole = 0; ihole != geometry.holes.size(); ++ihole) {
  //   plcFacetPolygons.resize(nfacets + geometry.holes[ihole].size());
  //   for (unsigned ifacet = 0; ifacet != geometry.holes[ihole].size(); ++ifacet) {
  //     POLY_ASSERT(geometry.holes[ihole][ifacet].size() >= 3);
  //     plcFacetPolygons[nfacets + ifacet].numberofvertices = geometry.holes[ihole][ifacet].size();
  //     plcFacetPolygons[nfacets + ifacet].vertexlist = const_cast<int*>(&geometry.holes[ihole][ifacet].front());
  //     incrementPosition(holeCentroids[ihole], PLCpoints, geometry.holes[ihole][ifacet]);
  //   }
  //   holeCentroids[ihole] /= RealType(geometry.holes[ihole].size());
  //   nfacets = plcFacetPolygons.size();
  // }
  // POLY_ASSERT(nfacets == plcFacetPolygons.size());
  // vector<tetgenio::facet> plcFacets(nfacets);
  // for (unsigned ifacet = 0; ifacet != nfacets; ++ifacet) {
  //   in.init(&plcFacets[ifacet]);
  //   plcFacets[ifacet].polygonlist = &plcFacetPolygons[ifacet];
  //   plcFacets[ifacet].numberofpolygons = 1;
  // }
  // in.numberoffacets = nfacets;
  // in.facetlist = &plcFacets.front();

  // Do the tetrahedralization.
  tetgenio out;
  tetrahedralize((char*)"vV", &in, &out);

  // Make sure we got something.
  if (out.numberoftetrahedra == 0)
    error("TetgenTessellator: Delauney tetrahedralization produced 0 tetrahedra!");
  if (out.numberofpoints != generators.size()/3) {
    char err[1024];
    snprintf(err, 1024, "TetgenTessellator: Delauney tetrahedralization produced %d tetrahedra\n(%d generating points given)", 
             out.numberofpoints, (int)numGenerators);
    error(err);
  }

  cerr << "Finished tetgen." << endl;

  // Read out the vertices, converting to quantized (hashed) positions and
  // buiding the map of hash -> node ID.
  map<IntPoint, int> point2node;
  for (unsigned i = 0; i != out.numberofvpoints; ++i) {
    IntPoint p(out.vpointlist[3*i], out.vpointlist[3*i + 1], out.vpointlist[3*i + 2], 0.0, 0.0, 0.0, degeneracy);
    internal::addKeyToMap(p, point2node);
  }
  const unsigned nvertices = point2node.size();

  // Read the tetgen voro edges, converting the rays to edges terminating on our bounding box.
  vector<EdgeHash> edges;
  edges.reserve(out.numberofvedges);
  for (unsigned i = 0; i != out.numberofvedges; ++i) {
    tetgenio::voroedge& vedge = out.vedgelist[i];
    int n1 = vedge.v1, n2 = vedge.v2;
    POLY_ASSERT2(n1 >= 0 and n1 < nvertices, i << " " << n1 << " " << vedge.v2);
    POLY_ASSERT2(n2 == -1 or (n2 >= 0 and n2 < nvertices), i << " " << n1 << " " << n2);
    if (n2 == -1) {
      // This edge is actually a ray, so we will introduce a new point where the ray
      // intersects the boundary.
      RealType pnew[3];
      const bool test = geometry::rayBoxIntersection<RealType>(&out.vpointlist[3*n1], vedge.vnormal, low, high, degeneracy, pnew);
      POLY_ASSERT(test);
      IntPoint ipnew(pnew[0], pnew[1], pnew[2], 0.0, 0.0, 0.0, degeneracy);
      n2 = internal::addKeyToMap(ipnew, point2node);
    }
    edges.push_back(internal::hashEdge(n1, n2));
  }
  POLY_ASSERT(edges.size() == out.numberofvedges);

  // Read the face info.
  const unsigned nfaces = out.numberofvfacets;
  mesh.faces.resize(nfaces);
  mesh.faceCells.resize(nfaces);
  for (unsigned i = 0; i != nfaces; ++i) {
    const tetgenio::vorofacet& vfacet = out.vfacetlist[i];
    const unsigned ne = vfacet.elist[0];
    POLY_ASSERT(ne >= 3);
    mesh.faces[i].reserve(ne);
    for (unsigned k = 0; k != ne; ++k) {
      const unsigned ie = vfacet.elist[k + 1];
      POLY_ASSERT(ie < edges.size());
      int n1 = edges[ie].first, n2 = edges[ie].second;
      if (k == 0) {
        mesh.faces[i].push_back(n1);
        mesh.faces[i].push_back(n2);
      } else if (k < ne - 1) {
        POLY_ASSERT(n1 == mesh.faces[i].back() or
                    n2 == mesh.faces[i].back());
        mesh.faces[i].push_back(n1 == mesh.faces[i].back() ?
                                n2 :
                                n1);
      } else {
        POLY_ASSERT((n1 == mesh.faces[i].back() and n2 == mesh.faces[i].front()) or
                    (n2 == mesh.faces[i].back() and n1 == mesh.faces[i].front()));
      }
    }
    POLY_ASSERT(mesh.faces[i].size() == ne);

    // Now the cells that touch this face.
    POLY_ASSERT(vfacet.c1 >= 0 and vfacet.c1 < numGenerators);
    POLY_ASSERT(vfacet.c2 >= 0 and vfacet.c2 < numGenerators);
    mesh.faceCells[i].push_back(vfacet.c1);
    mesh.faceCells[i].push_back(vfacet.c2);
  }

}

// //------------------------------------------------------------------------------
// // Tessellate within a PLC boundary.
// //------------------------------------------------------------------------------
// void
// TetgenTessellator::
// tessellate(const vector<double>& points,
//            const vector<double>& PLCpoints,
//            const PLC<3, double>& geometry,
//            Tessellation<3, double>& mesh) const {
//   if (mDirectComputation) {
//     this->computeDirectVoronoi(points, PLCpoints, geometry, mesh);
//   } else {
//     this->computeDualOfTetrahedralization(points, PLCpoints, geometry, mesh);
//   }
// }

// //------------------------------------------------------------------------------
// // Internal method to build the tessellation directly using Tetgen's native 
// // Voronoi capability.
// //------------------------------------------------------------------------------
// void
// TetgenTessellator::
// computeDirectVoronoi(const vector<double>& points,
//                      const vector<double>& PLCpoints,
//                      const PLC<3, double>& geometry,
//                      Tessellation<3, double>& mesh) const {

//   // Pre-conditions.
//   POLY_ASSERT(not points.empty());
//   POLY_ASSERT(points.size() % 3 == 0);

//   cout << "TetgenTessellator recieved PLC: " << endl << geometry << endl;

//   typedef int64_t CoordHash;
//   typedef set<int> FaceHash;
//   typedef Point3<CoordHash> IntPoint;
//   typedef Point3<RealType> RealPoint;
//   const CoordHash coordMax = (1LL << 34); // numeric_limits<CoordHash>::max() >> 32U;
//   const double degeneracy = 1.0e-12;

//   // Compute the normalized generators.
//   RealType low[3], high[3], box[3];
//   const unsigned numGenerators = points.size() / 3;
//   vector<double> generators = this->computeNormalizedPoints(points, PLCpoints, low, high);
//   for (unsigned j = 0; j != 3; ++j) {
//     box[j] = high[j] - low[j];
//     POLY_ASSERT(box[j] > 0.0);
//   }

//   // Build the input to tetgen.
//   tetgenio in;
//   in.firstnumber = 0;
//   in.mesh_dim = 3;
//   in.pointlist = &generators.front();
//   in.pointattributelist = 0;
//   in.pointmtrlist = 0;
//   in.pointmarkerlist = 0;
//   in.numberofpoints = generators.size() / 3;
//   in.numberofpointattributes = 0;
//   in.numberofpointmtrs = 0;

//   // // Copy the PLC boundary info to the tetgen input.
//   // vector<tetgenio::polygon> plcFacetPolygons(geometry.facets.size());
//   // for (unsigned ifacet = 0; ifacet != geometry.facets.size(); ++ifacet) {
//   //   POLY_ASSERT(geometry.facets[ifacet].size() >= 3);
//   //   plcFacetPolygons[ifacet].numberofvertices = geometry.facets[ifacet].size();
//   //   plcFacetPolygons[ifacet].vertexlist = const_cast<int*>(&geometry.facets[ifacet].front());
//   // }
//   // unsigned nfacets = plcFacetPolygons.size();
//   // vector<RealPoint> holeCentroids(geometry.holes.size(), RealPoint(0.0, 0.0, 0.0));
//   // for (unsigned ihole = 0; ihole != geometry.holes.size(); ++ihole) {
//   //   plcFacetPolygons.resize(nfacets + geometry.holes[ihole].size());
//   //   for (unsigned ifacet = 0; ifacet != geometry.holes[ihole].size(); ++ifacet) {
//   //     POLY_ASSERT(geometry.holes[ihole][ifacet].size() >= 3);
//   //     plcFacetPolygons[nfacets + ifacet].numberofvertices = geometry.holes[ihole][ifacet].size();
//   //     plcFacetPolygons[nfacets + ifacet].vertexlist = const_cast<int*>(&geometry.holes[ihole][ifacet].front());
//   //     incrementPosition(holeCentroids[ihole], PLCpoints, geometry.holes[ihole][ifacet]);
//   //   }
//   //   holeCentroids[ihole] /= RealType(geometry.holes[ihole].size());
//   //   nfacets = plcFacetPolygons.size();
//   // }
//   // POLY_ASSERT(nfacets == plcFacetPolygons.size());
//   // vector<tetgenio::facet> plcFacets(nfacets);
//   // for (unsigned ifacet = 0; ifacet != nfacets; ++ifacet) {
//   //   in.init(&plcFacets[ifacet]);
//   //   plcFacets[ifacet].polygonlist = &plcFacetPolygons[ifacet];
//   //   plcFacets[ifacet].numberofpolygons = 1;
//   // }
//   // in.numberoffacets = nfacets;
//   // in.facetlist = &plcFacets.front();

//   // Do the tetrahedralization.
//   tetgenio out;
//   tetrahedralize((char*)"vV", &in, &out);

//   // Make sure we got something.
//   if (out.numberoftetrahedra == 0)
//     error("TetgenTessellator: Delauney tetrahedralization produced 0 tetrahedra!");
//   if (out.numberofpoints != generators.size()/3) {
//     char err[1024];
//     snprintf(err, 1024, "TetgenTessellator: Delauney tetrahedralization produced %d tetrahedra\n(%d generating points given)", 
//              out.numberofpoints, (int)numGenerators);
//     error(err);
//   }

//   cerr << "Finished tetgen." << endl;

//   // Read out the vertices, converting to quantized (hashed) positions and
//   // buiding the map of hash -> node ID.
//   map<IntPoint, int> point2node;
//   for (unsigned i = 0; i != out.numberofvpoints; ++i) {
//     IntPoint p(out.vpointlist[3*i], out.vpointlist[3*i + 1], out.vpointlist[3*i + 2], 0.0, 0.0, 0.0, degeneracy);
//     internal::addKeyToMap(p, point2node);
//   }

//   // Read the face info.
//   const unsigned nfaces = out.numberofvfacets;
//   mesh.faces.resize(nfaces);
//   mesh.faceCells.resize(nfaces);
//   for (unsigned i = 0; i != nfaces; ++i) {
//     const tetgenio::vorofacet& vfacet = out.vfacetlist[i];
//     const unsigned ne = vfacet.elist[0];
//     POLY_ASSERT(ne >= 3);
//     mesh.faces[i].reserve(ne);
//     for (unsigned k = 0; k != ne; ++k) {
//       const unsigned ie = vfacet.elist[k];
//       POLY_ASSERT(ie < out.numberofvedges);
//       const tetgenio::voroedge& vedge = out.vedgelist[ie];
//       int n1 = vedge.v1;
//       int n2 = vedge.v2;
//       POLY_ASSERT2(n1 >= 0 and n1 < nvertices, ie << " " << n1 << " " << vedge.v2);
//       POLY_ASSERT2(n2 == -1 or (n2 >= 0 and n2 < nvertices), ie << " " << n1 << " " << n2);
//       if (k == 0) {
//         mesh.faces[i].push_back(n1);
//         mesh.faces[i].push_back(n2);
//       } else if (k < ne - 1) {
//         POLY_ASSERT(n1 == mesh.faces[i].back() or
//                     n2 == mesh.faces[i].back());
//         mesh.faces[i].push_back(n1 == mesh.faces[i].back() ?
//                                 n2 :
//                                 n1);
//       } else {
//         POLY_ASSERT((n1 == mesh.faces[i].back() and n2 == mesh.faces[i].front()) or
//                     (n2 == mesh.faces[i].back() and n1 == mesh.faces[i].front()));
//       }
//     }
//     POLY_ASSERT(mesh.faces[i].size() == ne);

//     // Now the cells that touch this face.
//     POLY_ASSERT(vfacet.c1 >= 0 and vfacet.c1 < numGenerators);
//     POLY_ASSERT(vfacet.c2 >= 0 and vfacet.c2 < numGenerators);
//     mesh.faceCells[i].push_back(vfacet.c1);
//     mesh.faceCells[i].push_back(vfacet.c2);
//   }
// }

// //------------------------------------------------------------------------------
// // Internal method to build the tessellation by first doing the 
// // tetrahedralization, and then computing the dual.
// //------------------------------------------------------------------------------
// void
// TetgenTessellator::
// computeDualOfTetrahedralization(const vector<double>& points,
//                                 const vector<double>& PLCpoints,
//                                 const PLC<3, double>& geometry,
//                                 Tessellation<3, double>& mesh) const {

//   // Pre-conditions.
//   POLY_ASSERT(not points.empty());
//   POLY_ASSERT(points.size() % 3 == 0);

//   typedef int64_t CoordHash;
//   typedef set<int> FaceHash;
//   typedef Point3<CoordHash> IntPoint;
//   typedef Point3<RealType> RealPoint;

//   // // Normalize coordinates in the input box.
//   // vector<double> generators;
//   // const unsigned numGenerators = points.size() / 3;
//   // generators.reserve(3*(numGenerators + 8));
//   // double box[3] = {high[0] - low[0],
//   //                  high[1] - low[1],
//   //                  high[2] - low[2]};
//   // POLY_ASSERT(box[0] > 0.0 and box[1] > 0.0 and box[2] > 0.0);
//   // for (unsigned i = 0; i != points.size(); ++i) {
//   //   const unsigned j = i % 3;
//   //   generators.push_back((points[i] - low[j])/box[j]);
//   //   POLY_ASSERT(generators.back() >= 0.0 and generators.back() <= 1.0);
//   // }
//   // POLY_ASSERT(generators.size() == points.size());

//   // // Add the boundary points to the generator list.
//   // generators.push_back(0.0); generators.push_back(0.0); generators.push_back(0.0);
//   // generators.push_back(1.0); generators.push_back(0.0); generators.push_back(0.0);
//   // generators.push_back(1.0); generators.push_back(1.0); generators.push_back(0.0);
//   // generators.push_back(0.0); generators.push_back(1.0); generators.push_back(0.0);
//   // generators.push_back(0.0); generators.push_back(0.0); generators.push_back(1.0);
//   // generators.push_back(1.0); generators.push_back(0.0); generators.push_back(1.0);
//   // generators.push_back(1.0); generators.push_back(1.0); generators.push_back(1.0);
//   // generators.push_back(0.0); generators.push_back(1.0); generators.push_back(1.0);
//   // POLY_ASSERT(generators.size() == 3*(numGenerators + 8));

//   // // Build the input to tetgen.
//   // tetgenio in;
//   // in.firstnumber = 0;
//   // in.mesh_dim = 3;
//   // in.pointlist = &generators.front();
//   // in.pointattributelist = 0;
//   // in.pointmtrlist = 0;
//   // in.pointmarkerlist = 0;
//   // in.numberofpoints = generators.size() / 3;
//   // in.numberofpointattributes = 0;
//   // in.numberofpointmtrs = 0;

//   // // Do the tetrahedralization.
//   // tetgenio out;
//   // tetrahedralize((char*)"qV", &in, &out);

//   // // Make sure we got something.
//   // if (out.numberoftetrahedra == 0)
//   //   error("TetgenTessellator: Delauney tetrahedralization produced 0 tetrahedra!");
//   // if (out.numberofpoints != generators.size()/3) {
//   //   char err[1024];
//   //   snprintf(err, 1024, "TetgenTessellator: Delauney tetrahedralization produced %d tetrahedra\n(%d generating points given)", 
//   //            out.numberofpoints, (int)numGenerators);
//   //   error(err);
//   // }

//   // // Find the circumcenters of each tetrahedron.
//   // RealType  clow[3] = { numeric_limits<RealType>::max(),  
//   //                       numeric_limits<RealType>::max(),
//   //                       numeric_limits<RealType>::max()};
//   // RealType chigh[3] = {-numeric_limits<RealType>::max(), 
//   //                      -numeric_limits<RealType>::max(),
//   //                      -numeric_limits<RealType>::max()};
//   // CounterMap<FaceHash> faceCounter;
//   // vector<RealPoint> circumcenters(out.numberoftetrahedra);
//   // map<int, set<int> > gen2tet;
//   // map<int, set<int> > neighbors;
//   // unsigned p1, p2, p3, p4;
//   // for (unsigned i = 0; i != out.numberoftetrahedra; ++i) {
//   // }
// }

// //------------------------------------------------------------------------------
// // Access the internal attribute to compute directly or not.
// //------------------------------------------------------------------------------
// bool
// TetgenTessellator::
// directComputation() const {
//   return mDirectComputation;
// }

// void
// TetgenTessellator::
// directComputation(const bool x) {
//   mDirectComputation = x;
// }

}
