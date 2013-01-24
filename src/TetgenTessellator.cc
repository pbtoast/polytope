//------------------------------------------------------------------------
// TetgenTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>
#include <set>

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

    // BLAGO!
    cerr << "Input edges :";
    for (unsigned i = 0; i != edges.size(); ++i) cerr << " (" << edges[i].first << " " << edges[i].second << ")";
    cerr << endl << "nodes2edges: " << endl;
    for (std::map<int, std::set<EdgeHash> >::const_iterator itr = nodes2edges.begin();
         itr != nodes2edges.end();
         ++itr) {
      cerr << itr->first << " : ";
      for (std::set<EdgeHash>::const_iterator eitr = itr->second.begin();
           eitr != itr->second.end();
           ++eitr) cerr << " (" << eitr->first << " " << eitr->second << ")";
      cerr << endl;
    }
    cerr << "nodeUseCount: " << endl;
    for (internal::CounterMap<int>::const_iterator itr = nodeUseCount.begin();
         itr != nodeUseCount.end();
         ++itr) {
      cerr << "   " << itr->first << " : " << itr->second << endl;
    }
    // BLAGO!

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
        lastNode = (nodeUseCount[edges[i].first] == 1 ? edges[i].first : edges[i].second);
        nodes2edges[lastNode].erase(orderedEdges.back());
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
      orderedEdges.push_back(edges[0]);
      nodes2edges[lastNode].erase(edges[0]);
    }

    // Walk the remaining edges
    while (orderedEdges.size() != nedges) {
      POLY_ASSERT(nodes2edges[lastNode].size() > 0);
      orderedEdges.push_back(*nodes2edges[lastNode].begin());
      nodes2edges[lastNode].erase(orderedEdges.back());
      lastNode = (orderedEdges.back().first == lastNode ? orderedEdges.back().second : orderedEdges.back().first);
    }
    
    // BLAGO!
    cerr << "Sorted edges : ";
    for (i = 0; i != nedges; ++i) cerr << " (" << orderedEdges[i].first << " " << orderedEdges[i].second << ")";
    cerr << endl;
    // BLAGO!

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
      cerr << "Looking at " << i << " " << j << endl;
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
// Default constructor.
//------------------------------------------------------------------------------
TetgenTessellator::
TetgenTessellator():
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
  in.pointlist = new double[generators.size()];
  copy(&generators.front(), &generators.front() + generators.size(), in.pointlist);
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
  vlow[0] = min(vlow[0], low[0]); vlow[1] = min(vlow[1], low[1]); vlow[2] = min(vlow[2], low[2]);
  vhigh[0] = max(vhigh[0], high[0]); vhigh[1] = max(vhigh[1], high[1]); vhigh[2] = max(vhigh[2], high[2]);
  const RealType vbox[3] = {vhigh[0] - vlow[0], vhigh[1] - vlow[1], vhigh[2] - vlow[2]};
  const RealType rinf = 2.0*max(vbox[0], max(vbox[1], vbox[2]));
  const RealType vboxc[3] = {0.5*(vlow[0] + vhigh[0]), 0.5*(vlow[1] + vhigh[1]), 0.5*(vlow[2] + vhigh[2])};
  vlow[0] = vboxc[0] - rinf; vlow[1] = vboxc[1] - rinf; vlow[2] = vboxc[2] - rinf;
  vhigh[0] = vboxc[0] + rinf; vhigh[1] = vboxc[1] + rinf; vhigh[2] = vboxc[2] + rinf;

  // Read out the vertices, converting to quantized (hashed) positions and
  // building the map of hash -> node ID.
  // Note the Tetgen will produce degenerate node positions, so this collapse of those
  // degeneracies is necessary!
  map<IntPoint, int> point2node;
  for (i = 0; i != out.numberofvpoints; ++i) {
    POLY_ASSERT(out.vpointlist[3*i  ] >= vlow[0] and out.vpointlist[3*i  ] <= vhigh[0]);
    POLY_ASSERT(out.vpointlist[3*i+1] >= vlow[1] and out.vpointlist[3*i+1] <= vhigh[1]);
    POLY_ASSERT(out.vpointlist[3*i+2] >= vlow[2] and out.vpointlist[3*i+2] <= vhigh[2]);
    IntPoint p(out.vpointlist[3*i], out.vpointlist[3*i+1], out.vpointlist[3*i+2], vlow[0], vlow[1], vlow[2], degeneracy);
    internal::addKeyToMap(p, point2node);
    // cerr << "   ----> (" << out.vpointlist[3*i] << " " << out.vpointlist[3*i + 1] << " " << out.vpointlist[3*i + 2] << ") " << p << " " << point2node[p] << endl;
  }

  // Read the normalized mesh node positions.
  mesh.nodes.resize(3*point2node.size());
  for (map<IntPoint, int>::const_iterator itr = point2node.begin();
       itr != point2node.end();
       ++itr) {
    i = itr->second;
    mesh.nodes[3*i  ] = itr->first.realx(vlow[0], degeneracy);
    mesh.nodes[3*i+1] = itr->first.realy(vlow[1], degeneracy);
    mesh.nodes[3*i+2] = itr->first.realz(vlow[2], degeneracy);
  }

  // Read the Tetgen voro edges, creating our "infinte" surface nodes for
  // unbounded rays.
  // Note that Tetgen will also produce degenerate edges between the degenerate
  // nodes!
  int v1, v2;
  unsigned n1, n2;
  RealType ray_sph_int[3];
  vector<EdgeHash> edges;
  edges.reserve(out.numberofvedges);
  mesh.infNodes = vector<unsigned>();
  for (i = 0; i != out.numberofvedges; ++i) {
    tetgenio::voroedge& vedge = out.vedgelist[i];
    v1 = vedge.v1;
    v2 = vedge.v2;
    POLY_ASSERT2(v1 >= 0 and v1 < out.numberofvpoints, i << " " << v1 << " " << v2 << " " << out.numberofvpoints);
    POLY_ASSERT2(v2 == -1 or (v2 >= 0 and v2 < out.numberofvpoints), i << " " << v1 << " " << v2 << " " << out.numberofvpoints);
    IntPoint p1(out.vpointlist[3*v1], out.vpointlist[3*v1+1], out.vpointlist[3*v1+2], vlow[0], vlow[1], vlow[2], degeneracy);
    POLY_ASSERT(point2node.find(p1) != point2node.end());
    n1 = point2node[p1];
    if (v2 == -1) {
      // This edge is a ray, so we construct a node on the rinf spherical surface
      // to hook to.
      geometry::raySphereIntersection(&mesh.nodes[3*n1],
                                      out.vedgelist[i].vnormal,
                                      vboxc,
                                      rinf,
                                      degeneracy,
                                      ray_sph_int);
      POLY_ASSERT(ray_sph_int[0] >= vlow[0] and ray_sph_int[0] <= vhigh[0]);
      POLY_ASSERT(ray_sph_int[1] >= vlow[1] and ray_sph_int[1] <= vhigh[1]);
      POLY_ASSERT(ray_sph_int[2] >= vlow[2] and ray_sph_int[2] <= vhigh[2]);
      IntPoint p1(ray_sph_int[0], ray_sph_int[1], ray_sph_int[2], vlow[0], vlow[1], vlow[2], degeneracy);
      k = point2node.size();
      n2 = internal::addKeyToMap(p1, point2node);
      if (k != point2node.size()) {
        mesh.infNodes.push_back(n2);
        mesh.nodes.push_back(p1.realx(vlow[0], degeneracy));
        mesh.nodes.push_back(p1.realy(vlow[1], degeneracy));
        mesh.nodes.push_back(p1.realz(vlow[2], degeneracy));
        // cerr << " ---> new inf point (" << ray_sph_int[0] << " " << ray_sph_int[1] << " " << ray_sph_int[2] << ") " << p1 << " " << n2 << endl;
      }
    } else {
      IntPoint p2(out.vpointlist[3*v2], out.vpointlist[3*v2+1], out.vpointlist[3*v2+2], vlow[0], vlow[1], vlow[2], degeneracy);
      POLY_ASSERT(point2node.find(p2) != point2node.end());
      n2 = point2node[p2];
    }
    if (n1 == n2) {   // Flag degenerate edges with negative values.
      n1 = -n1;
      n2 = -n1 - 1;
    }
    edges.push_back(internal::hashEdge(n1, n2));
  }
  POLY_ASSERT(edges.size() == out.numberofvedges);

  // Add any of our new "inf" nodes to the mesh node list.
  const unsigned nold = mesh.nodes.size()/3;
  mesh.nodes.resize(3*point2node.size());
  for (map<IntPoint, int>::const_iterator itr = point2node.begin();
       itr != point2node.end();
       ++itr) {
    i = itr->second;
    if (i >= nold) {
      mesh.nodes[3*i  ] = itr->first.realx(vlow[0], degeneracy);
      mesh.nodes[3*i+1] = itr->first.realy(vlow[1], degeneracy);
      mesh.nodes[3*i+2] = itr->first.realz(vlow[2], degeneracy);
    }
  }

  // Create the faces, checking for degeneracies here too.
  map<unsigned, int> faceMap;
  mesh.faces.reserve(out.numberofvfacets);
  mesh.faceCells.resize(out.numberofvfacets);
  int ne, ie, ie_last;
  for (i = 0; i != out.numberofvfacets; ++i) {
    const tetgenio::vorofacet& vfacet = out.vfacetlist[i];
    ne = vfacet.elist[0];
    POLY_ASSERT(ne >= 2);
    vector<EdgeHash> faceEdges;
    faceEdges.reserve(ne);
    for (k = 0; k != ne; ++k) {
      ie = vfacet.elist[k + 1];
      POLY_ASSERT(ie < edges.size());
      if (edges[ie].first >= 0) faceEdges.push_back(edges[ie]); // Exclude degenerate edges.
    }

    // Reduce to the unique edges, and extrac the nodes ringing the face.
    sort(faceEdges.begin(), faceEdges.end());
    faceEdges.erase(unique(faceEdges.begin(), faceEdges.end()), faceEdges.end());
    const vector<unsigned> faceNodes = computeSortedFaceNodes(faceEdges);

    // Are there enough nodes to constitute a non-degenerate face?
    if (faceNodes.size() > 2) {
      faceMap[i] = mesh.faces.size();
      mesh.faces.push_back(faceNodes);

      // Now the cells that touch this face.
      POLY_ASSERT(vfacet.c1 >= 0 and vfacet.c1 < numGenerators);
      POLY_ASSERT(vfacet.c2 >= 0 and vfacet.c2 < numGenerators);
      mesh.faceCells[i].push_back(vfacet.c1);
      mesh.faceCells[i].push_back(vfacet.c2);
    } else {
      faceMap[i] = -1;
    }
  }
  const unsigned nfaces = mesh.faces.size();
  POLY_ASSERT(faceMap.size() == out.numberofvfacets);

  // Read out the cell structure as collections of faces.
  POLY_ASSERT(out.numberofvcells == numGenerators);
  mesh.cells = vector<vector<int> >(numGenerators);
  unsigned nf;
  int fi;
  RealType ccent[3], fcent[3], fhat[3];
  for (i = 0; i != numGenerators; ++i) {
    nf = out.vcelllist[i][0];
    POLY_ASSERT(nf >= 4);
    for (k = 0; k != nf; ++k) {
      POLY_ASSERT(faceMap.find(out.vcelllist[i][k+1]) != faceMap.end());
      fi = faceMap[out.vcelllist[i][k+1]];
      if (fi >= 0) mesh.cells[i].push_back(fi);
    }
    POLY_ASSERT(mesh.cells[i].size() <= nf);
    nf = mesh.cells[i].size();

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

  // Start by creating the unbounded tessellation.
  this->tessellate(points, mesh);
}

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

}
