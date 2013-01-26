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

// Returns (positive, 0.0, negative) if pd is (below, coplanar, above) the plane
// (pa, pb, pc), where above is defined such that (pa, pb, pc) is counter-clockwise.
// This puppy is defined in predicates.cc
extern double orient3d(double* pa, double* pb, double* pc, double* pd);

namespace polytope {

using namespace std;
using std::min;
using std::max;
using std::abs;

namespace {

//------------------------------------------------------------------------------
// Borrow the Point3 type as a tuple to create 3 node facets hashes.
//------------------------------------------------------------------------------
Point3<unsigned>
hashFacet(const unsigned i, const unsigned j, const unsigned k) {
  typedef Point3<unsigned> Tuple3;
  POLY_ASSERT(i != j and i != k and j != k);
  if (i < j and i < k) {
    if (j < k) {
      return Tuple3(i, j, k);
    } else {
      return Tuple3(i, k, j);
    }
  } else if (j < i and j < k) {
    if (i < k) {
      return Tuple3(j, i, k);
    } else {
      return Tuple3(j, k, i);
    }
  } else {
    if (i < j) {
      return Tuple3(k, i, j);
    } else {
      return Tuple3(k, j, i);
    }
  }
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

    // // BLAGO!
    // cerr << "Input edges :";
    // for (unsigned i = 0; i != edges.size(); ++i) cerr << " (" << edges[i].first << " " << edges[i].second << ")";
    // cerr << endl << "nodes2edges: " << endl;
    // for (std::map<int, std::set<EdgeHash> >::const_iterator itr = nodes2edges.begin();
    //      itr != nodes2edges.end();
    //      ++itr) {
    //   cerr << itr->first << " : ";
    //   for (std::set<EdgeHash>::const_iterator eitr = itr->second.begin();
    //        eitr != itr->second.end();
    //        ++eitr) cerr << " (" << eitr->first << " " << eitr->second << ")";
    //   cerr << endl;
    // }
    // cerr << "nodeUseCount: " << endl;
    // for (internal::CounterMap<int>::const_iterator itr = nodeUseCount.begin();
    //      itr != nodeUseCount.end();
    //      ++itr) {
    //   cerr << "   " << itr->first << " : " << itr->second << endl;
    // }
    // // BLAGO!

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
        nodes2edges[edges[i].first].erase(edges[i]);
        nodes2edges[edges[i].second].erase(edges[i]);
        lastNode = (nodeUseCount[edges[i].first] == 1 ? edges[i].first : edges[i].second);
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
      // orderedEdges.push_back(edges[0]);
      // nodes2edges[lastNode].erase(edges[0]);
    }

    // Walk the remaining edges
    while (orderedEdges.size() != nedges) {
      POLY_ASSERT(nodes2edges[lastNode].size() > 0);
      orderedEdges.push_back(*nodes2edges[lastNode].begin());
      nodes2edges[orderedEdges.back().first].erase(orderedEdges.back());
      nodes2edges[orderedEdges.back().second].erase(orderedEdges.back());
      lastNode = (orderedEdges.back().first == lastNode ? orderedEdges.back().second : orderedEdges.back().first);
    }
    
    // // BLAGO!
    // cerr << "Sorted edges : ";
    // for (i = 0; i != nedges; ++i) cerr << " (" << orderedEdges[i].first << " " << orderedEdges[i].second << ")";
    // cerr << endl;
    // // BLAGO!

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
// Constructor.
//------------------------------------------------------------------------------
TetgenTessellator::
TetgenTessellator(const bool directComputation):
  Tessellator<3, double>(),
  mDirectComputation(directComputation) {
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
  if (mDirectComputation) {
    this->computeVoronoiNatively(points, mesh);
  } else {
    this->computeVoronoiThroughTetrahedralization(points, mesh);
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

//------------------------------------------------------------------------------
// Compute the unbouded tessellation using Tetgen's directo Voronoi computation.
//------------------------------------------------------------------------------
void
TetgenTessellator::
computeVoronoiNatively(const vector<double>& points,
                       Tessellation<3, double>& mesh) const {

  // Pre-conditions.
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 3 == 0);

  typedef int64_t CoordHash;
  typedef pair<int, int> EdgeHash;
  typedef Point3<CoordHash> IntPoint;
  typedef Point3<RealType> RealPoint;
  const CoordHash coordMax = (1LL << 34); // numeric_limits<CoordHash>::max() >> 32U;
  const double degeneracy = 1.0/coordMax;

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

  // Do the tetgen tessellation.
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
  int ne, ie;
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
// Use tetgen to compute the Delaunay tetrathedralization, and do the dual 
// ourselves.  This is the way Tetgen is more commonly used, so maybe safer.
//------------------------------------------------------------------------------
void
TetgenTessellator::
computeVoronoiThroughTetrahedralization(const vector<double>& points,
                                        Tessellation<3, double>& mesh) const {

  // Pre-conditions.
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 3 == 0);

  typedef int64_t CoordHash;
  typedef pair<int, int> EdgeHash;
  typedef Point3<CoordHash> IntPoint;
  typedef Point3<RealType> RealPoint;
  typedef Point3<unsigned> TetFacetHash;  // kind of nefarious!
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
  tetrahedralize((char*)"V", &in, &out);

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

  // Compute the circumcenters of the tetrahedra, and the set of tets associated
  // with each generator.
  RealType clow[3], chigh[3], X[3];
  vector<RealPoint> circumcenters(out.numberoftetrahedra);
  unsigned a, b, c, d;
  EdgeHash ab, ac, ad, bc, bd, cd;
  TetFacetHash abc, abd, bcd, acd;
  map<TetFacetHash, vector<unsigned> > facet2tets;      // Tets which share a facet.
  vector<vector<EdgeHash> > gen2edges(numGenerators);   // Edges from each generator.
  map<EdgeHash, set<unsigned> > edge2tets;              // Tets which share a given edge.
  for (i = 0; i != out.numberoftetrahedra; ++i) {
    a = out.tetrahedronlist[4*i];
    b = out.tetrahedronlist[4*i+1];
    c = out.tetrahedronlist[4*i+2];
    d = out.tetrahedronlist[4*i+3];
    POLY_ASSERT(a < numGenerators);
    POLY_ASSERT(b < numGenerators);
    POLY_ASSERT(c < numGenerators);
    POLY_ASSERT(d < numGenerators);
    geometry::computeCircumcenter3d(&out.pointlist[3*a],
                                    &out.pointlist[3*b],
                                    &out.pointlist[3*c],
                                    &out.pointlist[3*d],
                                    &circumcenters[i].x);
    ab = internal::hashEdge(a, b);
    ac = internal::hashEdge(a, c);
    ad = internal::hashEdge(a, d);
    bc = internal::hashEdge(b, c);
    bd = internal::hashEdge(b, d);
    cd = internal::hashEdge(c, d);
    abc = hashFacet(a, b, c);
    abd = hashFacet(a, b, d);
    bcd = hashFacet(b, c, d);
    acd = hashFacet(a, c, d);
    facet2tets[abc].push_back(i);
    facet2tets[abd].push_back(i);
    facet2tets[bcd].push_back(i);
    facet2tets[acd].push_back(i);
    gen2edges[i].push_back(ab);
    gen2edges[i].push_back(ac);
    gen2edges[i].push_back(ad);
    gen2edges[i].push_back(bc);
    gen2edges[i].push_back(bd);
    gen2edges[i].push_back(cd);
    edge2tets[ab].insert(i);
    edge2tets[ac].insert(i);
    edge2tets[ad].insert(i);
    edge2tets[bc].insert(i);
    edge2tets[bd].insert(i);
    edge2tets[cd].insert(i);
    clow[0] = min(clow[0], circumcenters.back().x);
    clow[1] = min(clow[1], circumcenters.back().y);
    clow[2] = min(clow[2], circumcenters.back().z);
    chigh[0] = max(chigh[0], circumcenters.back().x);
    chigh[1] = max(chigh[1], circumcenters.back().y);
    chigh[2] = max(chigh[2], circumcenters.back().z);
  }
  POLY_BEGIN_CONTRACT_SCOPE;
  {
    for (map<TetFacetHash, vector<unsigned> >::const_iterator itr = facet2tets.begin();
         itr != facet2tets.end();
         ++itr) POLY_ASSERT(itr->second.size() == 1 or itr->second.size() == 2);
    POLY_ASSERT(gen2edges.size() == numGenerators);
    for (unsigned i = 0; i != numGenerators; ++i) POLY_ASSERT(gen2edges[i].size() >= 3);
    for (map<EdgeHash, set<unsigned> >::const_iterator itr = edge2tets.begin();
         itr != edge2tets.end();
         ++itr) POLY_ASSERT(itr->second.size() >= 1);
    POLY_ASSERT(clow[0] < chigh[0] and clow[1] < chigh[1] and clow[2] < chigh[2]);
  }
  POLY_END_CONTRACT_SCOPE;

  // The bounding box for the circumcenters, useful for quantizing those pups.
  clow[0] = min(clow[0], low[0]); clow[1] = min(clow[1], low[1]); clow[2] = min(clow[2], low[2]);
  chigh[0] = max(chigh[0], high[0]); chigh[1] = max(chigh[1], high[1]); chigh[2] = max(chigh[2], high[2]);
  const RealType cbox[3]  = {chigh[0] - clow[0], chigh[1] - clow[1], chigh[2] - clow[2]};
  const RealType cboxsize = 2.0*max(cbox[0], max(cbox[1], cbox[2]));
  const RealType cdx = max(degeneracy, cboxsize/coordMax);

  // Any surface facets create new "infinite" or "unbounded" rays, which originate at
  // the tet circumcenter and pass through the circumcenter of the triangular facet.
  const RealType rinf = 2.0*cboxsize;
  const RealType cboxc[3] = {0.5*(clow[0] + chigh[0]), 0.5*(clow[1] + chigh[1]), 0.5*(clow[2] + chigh[2])};

  // Create the quantized circumcenters, and the map from the (possibly) degenerate
  // circumcenters to their unique IDs.
  map<IntPoint, int> circ2id;
  map<int, vector<unsigned> > tet2ids;
  for (i = 0; i != out.numberoftetrahedra; ++i) {
    IntPoint ip(circumcenters[i].x, circumcenters[i].y, circumcenters[i].z,
                clow[0], clow[1], clow[2], cdx);
    j = internal::addKeyToMap(ip, circ2id);
    tet2ids[i].push_back(j);
  }

  // Look for any surface facets we need to project unbounded rays through.
  bool test;
  RealPoint fhat, tetcent, test_point, a_b, a_c, pinf;
  for (map<TetFacetHash, vector<unsigned> >::const_iterator facetItr = facet2tets.begin();
       facetItr != facet2tets.end();
       ++facetItr) {
    const TetFacetHash& facet = facetItr->first;
    const vector<unsigned>& tets = facetItr->second;
    if (tets.size() == 1) {
      i = tets[0];
      POLY_ASSERT(i < out.numberoftetrahedra);
      a = out.tetrahedronlist[4*i];
      b = out.tetrahedronlist[4*i+1];
      c = out.tetrahedronlist[4*i+2];
      d = out.tetrahedronlist[4*i+3];
      POLY_ASSERT(a < numGenerators);
      POLY_ASSERT(b < numGenerators);
      POLY_ASSERT(c < numGenerators);
      POLY_ASSERT(d < numGenerators);
      geometry::computeTetCentroid(&out.pointlist[3*a],
                                   &out.pointlist[3*b],
                                   &out.pointlist[3*c],
                                   &out.pointlist[3*d],
                                   &tetcent.x);

      // We need the ray unit vector.
      test = geometry::computeTriangleCircumcenter3d(&out.pointlist[3*facet.x],
                                                     &out.pointlist[3*facet.y],
                                                     &out.pointlist[3*facet.z],
                                                     &fhat.x);
      POLY_ASSERT(test);
      fhat -= circumcenters[i];

      // Check for the special case of the tet circumcenter coplanar with the facet.
      if (abs(geometry::dot<3, RealType>(&fhat.x, &fhat.x)) < degeneracy) {
        // Yep, it's in the plane.  Just project the ray out orthogonally to the facet.
        a_b.x = out.pointlist[3*facet.y]   - out.pointlist[3*facet.x];
        a_b.y = out.pointlist[3*facet.y+1] - out.pointlist[3*facet.x+1];
        a_b.z = out.pointlist[3*facet.y+2] - out.pointlist[3*facet.x+2];
        a_c.x = out.pointlist[3*facet.z]   - out.pointlist[3*facet.x];
        a_c.y = out.pointlist[3*facet.z+1] - out.pointlist[3*facet.x+1];
        a_c.z = out.pointlist[3*facet.z+2] - out.pointlist[3*facet.x+2];
        geometry::cross<3, RealType>(&a_b.x, &a_c.x, &fhat.x);
      }
      geometry::unitVector<3, RealType>(&fhat.x);

      // The ray unit vector should point in the opposite direction from the facet as the tet centroid.
      POLY_ASSERT(abs(orient3d(&out.pointlist[3*facet.x], &out.pointlist[3*facet.y], &out.pointlist[3*facet.z], &tetcent.x)) > degeneracy);
      copy(&out.pointlist[3*facet.x], &out.pointlist[3*facet.x] + 3, &test_point.x);
      test_point += fhat;
      if (orient3d(&out.pointlist[3*facet.x], &out.pointlist[3*facet.y], &out.pointlist[3*facet.z], &tetcent.x)*
          orient3d(&out.pointlist[3*facet.x], &out.pointlist[3*facet.y], &out.pointlist[3*facet.z], &test_point.x) > 0.0) fhat *= -1.0;

      // Now we can compute the point where this ray intersects the surrounding "inf" sphere.
      test = geometry::raySphereIntersection(&circumcenters[i].x,
                                             &fhat.x,
                                             cboxc,
                                             rinf,
                                             1.0e-10,
                                             &pinf.x);
      POLY_ASSERT(test);
      IntPoint ip(pinf.x, pinf.y, pinf.z, clow[0], clow[1], clow[2], cdx);
      j = internal::addKeyToMap(ip, circ2id);
      tet2ids[i].push_back(j);
    }
  }

  // Copy the quantized nodes to the final tessellation.
  const unsigned numNodes = circ2id.size();
  mesh.nodes.resize(3*numNodes);
  for (map<IntPoint, int>::const_iterator itr = circ2id.begin();
       itr != circ2id.end();
       ++itr) {
    POLY_ASSERT(itr->second >= 0 and itr->second < numNodes);
    i = itr->second;
    mesh.nodes[3*i  ] = itr->first.realx(clow[0], cdx);
    mesh.nodes[3*i+1] = itr->first.realy(clow[1], cdx);
    mesh.nodes[3*i+2] = itr->first.realz(clow[2], cdx);
  }

  // Build the faces corresponding to each tet edge.

  // Rescale the mesh node positions for the input geometry.
  for (i = 0; i != mesh.nodes.size(); ++i) {
    j = i % 3;
    mesh.nodes[i] = low[j] + mesh.nodes[i]*box[j];
  }

}

}
