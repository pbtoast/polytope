//------------------------------------------------------------------------
// TetgenTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <limits>

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
void
computeSortedFaceNodes(const std::vector<std::pair<int, int> >& edges,
                       std::vector<unsigned>& result,
                       std::vector<unsigned>& hangingNodes) {
  typedef std::pair<int, int> EdgeHash;

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
    //   cerr << "   " << itr->first << " : ";
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

    // Look for any edges with one node in the set.  There can be at most
    // two such edges, representing the two ends of the chain.  We will put 
    // the edges with those nodes first in the ordering, so that all remaining
    // edges should naturally hook together.
    std::vector<EdgeHash> orderedEdges;
    orderedEdges.reserve(nedges);
    int lastNode;
    for (i = 0; i != nedges; ++i) {
      if (nodeUseCount[edges[i].first] == 1 or
          nodeUseCount[edges[i].second] == 1) {
        POLY_ASSERT((nodeUseCount[edges[i].first] == 1 and nodeUseCount[edges[i].second] == 2) or
                    (nodeUseCount[edges[i].first] == 2 and nodeUseCount[edges[i].second] == 1));
        orderedEdges.push_back(edges[i]);
        nodes2edges[edges[i].first].erase(edges[i]);
        nodes2edges[edges[i].second].erase(edges[i]);
        lastNode = (nodeUseCount[edges[i].first] == 1 ? edges[i].first : edges[i].second);
        hangingNodes.push_back(lastNode);
      }
    }
    POLY_ASSERT(orderedEdges.size() == 0 or orderedEdges.size() == 2);
    POLY_ASSERT(hangingNodes.size() == 0 or hangingNodes.size() == 2);

    // Pick a node to start the chain.
    if (hangingNodes.size() == 2) {
      POLY_ASSERT(nodeUseCount[orderedEdges.back().first] == 2 or
                  nodeUseCount[orderedEdges.back().second] == 2);
      lastNode = (nodeUseCount[orderedEdges.back().first] == 2 ? 
                  orderedEdges.back().first :
                  orderedEdges.back().second);
    } else {
      lastNode = edges[0].first;
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
    if (hangingNodes.size() == 2) {
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
    POLY_ASSERT2((hangingNodes.size() > 0 and result.size() == nedges + 1) or
                 result.size() == nedges, result.size());
      

  } else {

    // There are either one or no edges, so the solution is pretty simple.
    if (nedges == 1) {
      result.push_back(edges[0].first);
      result.push_back(edges[0].second);
    }

  }

  // That's it.
}

//------------------------------------------------------------------------------
// Sort a set of edges around a face so that sequential edges share nodes.
// We allow for one break in the chain (representing on unbounded surface).
// In such a situation we insert the new edge at the beginning of the chain, and
// return "true" indicating that a new edge was created.
//------------------------------------------------------------------------------
bool
computeSortedFaceEdges(std::vector<std::pair<int, int> >& edges,
                       std::vector<int>& result) {
  typedef std::pair<int, int> EdgeHash;

  unsigned nedges = edges.size();
  POLY_ASSERT(nedges >= 2);

  // Invert the mapping, from nodes to edges.
  std::map<int, std::set<unsigned> > nodes2edges;
  internal::CounterMap<int> nodeUseCount;
  unsigned i;
  for (i = 0; i != nedges; ++i) {
    nodes2edges[edges[i].first].insert(i);
    nodes2edges[edges[i].second].insert(i);
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
  //   cerr << "   " << itr->first << " : ";
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

  // Look for any edges with one node in the set.  There can be at most
  // two such edges, representing the two ends of the chain.  We introduce a
  // new edge hooking those hanging nodes together, and off we go.
  int lastNode;
  vector<int> hangingNodes;
  for (i = 0; i != nedges; ++i) {
    if (nodeUseCount[edges[i].first] == 1 or
        nodeUseCount[edges[i].second] == 1) {
      POLY_ASSERT((nodeUseCount[edges[i].first] == 1 and nodeUseCount[edges[i].second] == 2) or
                  (nodeUseCount[edges[i].first] == 2 and nodeUseCount[edges[i].second] == 1));
      result.push_back(i);
      nodes2edges[edges[i].first].erase(i);
      nodes2edges[edges[i].second].erase(i);
      lastNode = (nodeUseCount[edges[i].first] == 1 ? edges[i].first : edges[i].second);
      hangingNodes.push_back(lastNode);
    }
  }
  POLY_ASSERT(result.size() == 0 or (hangingNodes.size() == 2 and result.size() == 2));

  // If needed create that new edge and put it in the set.
  if (hangingNodes.size() == 2) {
    result.insert(result.begin() + 1, edges.size());
    edges.push_back(internal::hashEdge(hangingNodes[0], hangingNodes[1]));
    ++nedges;
    POLY_ASSERT(result.size() == 3);
  }
  POLY_ASSERT(edges.size() == nedges);

  // Pick a node to start the chain.
  if (hangingNodes.size() == 2) {
    POLY_ASSERT(nodeUseCount[edges[result.back()].first] == 2 or
                nodeUseCount[edges[result.back()].second] == 2);
    lastNode = (nodeUseCount[edges[result.back()].first] == 2 ? 
                edges[result.back()].first :
                edges[result.back()].second);
  } else {
    lastNode = edges[0].first;
  }

  // Walk the remaining edges
  EdgeHash ehash;
  while (result.size() != nedges) {
    POLY_ASSERT(nodes2edges[lastNode].size() > 0);
    result.push_back(*nodes2edges[lastNode].begin());
    ehash = edges[result.back()];
    nodes2edges[ehash.first].erase(result.back());
    nodes2edges[ehash.second].erase(result.back());
    lastNode = (ehash.first == lastNode ? ehash.second : ehash.first);
  }
  
  // // BLAGO!
  // cerr << "Sorted edges : ";
  // for (i = 0; i != nedges; ++i) cerr << " (" << orderedEdges[i].first << " " << orderedEdges[i].second << ")";
  // cerr << endl;
  // // BLAGO!

  // Set the orientation for the ordered edges.
  lastNode = (edges[result[0]].first == edges[result[1]].first ? edges[result[0]].first : edges[result[0]].second);
  for (i = 1; i != nedges; ++i) {
    POLY_ASSERT(edges[result[i]].first == lastNode or edges[result[i]].second == lastNode);
    if (edges[result[i]].first == lastNode) {
      lastNode = edges[result[i]].second;
    } else {
      lastNode = edges[result[i]].first;
      result[i] = ~result[i];
    }
  }

  // That's it.
  POLY_BEGIN_CONTRACT_SCOPE;
  {
    POLY_ASSERT(edges.size() == result.size());
    for (int i = 0; i != edges.size(); ++i) {
      const int j = (i + 1) % edges.size();
      const int ii = result[i];
      const int jj = result[j];
      POLY_ASSERT((ii >= 0 ? ii : ~ii) < edges.size());
      POLY_ASSERT((jj >= 0 ? jj : ~jj) < edges.size());
      POLY_ASSERT(((ii >= 0 and jj >= 0) and edges[ii].second == edges[jj].first) or
                  ((ii >= 0 and jj <  0) and edges[ii].second == edges[~jj].second) or
                  ((ii <  0 and jj >= 0) and edges[~ii].first == edges[jj].first) or
                  ((ii <  0 and jj <  0) and edges[~ii].first == edges[~jj].second));
    }
  }
  POLY_END_CONTRACT_SCOPE;
  return !(hangingNodes.empty());
}

//------------------------------------------------------------------------------
// Given an array of 4 integers and 2 unique values, find the other two.
//------------------------------------------------------------------------------
void
findOtherTetIndices(const int* indices,
                    const int a,
                    const int b,
                    int& c,
                    int& d) {
  POLY_ASSERT(a == indices[0] or a == indices[1] or a == indices[2] or a == indices[3]);
  POLY_ASSERT(b == indices[0] or b == indices[1] or b == indices[2] or b == indices[3]);
  POLY_ASSERT(indices[0] != indices[1] and indices[0] != indices[2] and indices[0] != indices[3] and
                                           indices[1] != indices[2] and indices[1] != indices[3] and
                                                                        indices[2] != indices[3]);
  if (a != indices[0] and b != indices[0]) {
    c = indices[0];
    d = (a != indices[1] and b != indices[1] ? indices[1] :
         a != indices[2] and b != indices[2] ? indices[2] :
         indices[3]);
  } else if (a != indices[1] and b != indices[1]) {
    c = indices[1];
    d = (a != indices[2] and b != indices[2] ? indices[2] :
         indices[3]);
  } else {
    c = indices[2];
    d = indices[3];
  }
}

} // end anonymous namespace

//------------------------------------------------------------------------------
// Constructor.
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

  // First generate our internal quantized tessellation representation.
  internal::QuantTessellation<3, double> qmesh;
  this->computeUnboundedQuantizedTessellation(points, qmesh);

  // Convert to the output tessellation and we're done.
  qmesh.tessellation(mesh);
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

  // Flag all the inf nodes in a fast lookup array.
  unsigned nnodes = mesh.nodes.size()/3;
  int i, j, k, j0, k0, jplane, kplane;
  //vector<unsigned> infNodeFlags(nnodes, 0);
  //for (i = 0; i != mesh.infNodes.size(); ++i) infNodeFlags[mesh.infNodes[i]] = 1;

  // Find the faces with "inf" nodes, and close them with edges along the
  // bounding box.
  RealPoint ray_hat;
  unsigned nfaces = mesh.faces.size();
  for (i = 0; i != nfaces; ++i) {
    vector<unsigned>& faceNodes = mesh.faces[i];
    j = 0;
    while (j < faceNodes.size() and mesh.infNodes[faceNodes[j]] == 0) ++j;
    if (j < faceNodes.size()) {

      // Yep, this face has inf nodes.  There should be two of 'em.
      k = (j + 1) % faceNodes.size();
      j0 = (j - 1) % faceNodes.size(); 
      k0 = (k - 1) % faceNodes.size(); 
      POLY_ASSERT(mesh.infNodes[faceNodes[j]] == 1 and
                  mesh.infNodes[faceNodes[k]] == 1 and
                  mesh.infNodes[faceNodes[j0]] == 0 and
                  mesh.infNodes[faceNodes[k0]] == 0);
      ray_hat.x = mesh.nodes[3*j]   - mesh.nodes[3*j0];
      ray_hat.y = mesh.nodes[3*j+1] - mesh.nodes[3*j0+1];
      ray_hat.z = mesh.nodes[3*j+2] - mesh.nodes[3*j0+2];
      geometry::unitVector<3, RealType>(&ray_hat.x);
      jplane = geometry::rayBoxIntersection(&mesh.nodes[3*j0],
                                            &ray_hat.x,
                                            low,
                                            high,
                                            mDegeneracy,
                                            &mesh.nodes[3*j]);
      ray_hat.x = mesh.nodes[3*k]   - mesh.nodes[3*k0];
      ray_hat.y = mesh.nodes[3*k+1] - mesh.nodes[3*k0+1];
      ray_hat.z = mesh.nodes[3*k+2] - mesh.nodes[3*k0+2];
      geometry::unitVector<3, RealType>(&ray_hat.x);
      kplane = geometry::rayBoxIntersection(&mesh.nodes[3*k0],
                                            &ray_hat.x,
                                            low,
                                            high,
                                            mDegeneracy,
                                            &mesh.nodes[3*k]);

      // // Do we need to introduce a new node for this face?
      // if (jplane != kplane) {
      //   POLY_ASSERT(
    }
  }

}

// //------------------------------------------------------------------------------
// // Use tetgen to compute the Delaunay tetrahedralization, and do the dual 
// // ourselves.  This is the way Tetgen is more commonly used, so maybe safer.
// //------------------------------------------------------------------------------
// void
// TetgenTessellator::
// computeVoronoiThroughTetrahedralization(const vector<double>& points,
//                                         Tessellation<3, double>& mesh) const {

//   // Pre-conditions.
//   POLY_ASSERT(not points.empty());
//   POLY_ASSERT(points.size() % 3 == 0);

//   typedef uint64_t PointHash;
//   typedef pair<int, int> EdgeHash;
//   typedef Point3<RealType> RealPoint;
//   typedef Point3<unsigned> TetFacetHash;  // kind of nefarious!
//   const double mDegeneracy = 1.0e-12;

//   // Compute the normalized generators.
//   RealType low[3], high[3], box[3];
//   const unsigned numGenerators = points.size() / 3;
//   vector<double> generators = this->computeNormalizedPoints(points, points, true, low, high);
//   unsigned i, j, k;
//   for (j = 0; j != 3; ++j) {
//     box[j] = high[j] - low[j];
//     POLY_ASSERT(box[j] > 0.0);
//   }

//   // Build the input to tetgen.
//   tetgenio in;
//   in.firstnumber = 0;
//   in.mesh_dim = 3;
//   in.pointlist = new double[generators.size()];
//   copy(&generators.front(), &generators.front() + generators.size(), in.pointlist);
//   in.pointattributelist = 0;
//   in.pointmtrlist = 0;
//   in.pointmarkerlist = 0;
//   in.numberofpoints = generators.size() / 3;
//   in.numberofpointattributes = 0;
//   in.numberofpointmtrs = 0;

//   // Do the tetrahedralization.
//   tetgenio out;
//   tetrahedralize((char*)"Q", &in, &out);
//   // tetrahedralize((char*)"V", &in, &out);

//   // Make sure we got something.
//   if (out.numberoftetrahedra == 0)
//     error("TetgenTessellator: Delauney tetrahedralization produced 0 tetrahedra!");
//   if (out.numberofpoints != generators.size()/3) {
//     char err[1024];
//     snprintf(err, 1024, "TetgenTessellator: Delauney tetrahedralization produced %d tetrahedra\n(%d generating points given)", 
//              out.numberofpoints, (int)numGenerators);
//     error(err);
//   }

//   // Compute the circumcenters of the tetrahedra, and the set of tets associated
//   // with each generator.
//   RealType clow[3] = {numeric_limits<RealType>::max(), 
//                       numeric_limits<RealType>::max(), 
//                       numeric_limits<RealType>::max()},
//           chigh[3] = {-numeric_limits<RealType>::max(), 
//                       -numeric_limits<RealType>::max(), 
//                       -numeric_limits<RealType>::max()};
//   vector<RealPoint> circumcenters(out.numberoftetrahedra);
//   int a, b, c, d;
//   EdgeHash ab, ac, ad, bc, bd, cd;
//   TetFacetHash abc, abd, bcd, acd;
//   map<TetFacetHash, vector<unsigned> > facet2tets;      // Tets which share a facet.
//   map<EdgeHash, set<unsigned> > edge2tets;              // Tets which share a given edge.
//   for (i = 0; i != out.numberoftetrahedra; ++i) {
//     a = out.tetrahedronlist[4*i];
//     b = out.tetrahedronlist[4*i+1];
//     c = out.tetrahedronlist[4*i+2];
//     d = out.tetrahedronlist[4*i+3];
//     POLY_ASSERT(a < numGenerators);
//     POLY_ASSERT(b < numGenerators);
//     POLY_ASSERT(c < numGenerators);
//     POLY_ASSERT(d < numGenerators);
//     geometry::computeCircumcenter3d(&out.pointlist[3*a],
//                                     &out.pointlist[3*b],
//                                     &out.pointlist[3*c],
//                                     &out.pointlist[3*d],
//                                     &circumcenters[i].x);
//     ab = internal::hashEdge(a, b);
//     ac = internal::hashEdge(a, c);
//     ad = internal::hashEdge(a, d);
//     bc = internal::hashEdge(b, c);
//     bd = internal::hashEdge(b, d);
//     cd = internal::hashEdge(c, d);
//     abc = hashFacet(a, b, c);
//     abd = hashFacet(a, b, d);
//     bcd = hashFacet(b, c, d);
//     acd = hashFacet(a, c, d);
//     facet2tets[abc].push_back(i);
//     facet2tets[abd].push_back(i);
//     facet2tets[bcd].push_back(i);
//     facet2tets[acd].push_back(i);
//     edge2tets[ab].insert(i);
//     edge2tets[ac].insert(i);
//     edge2tets[ad].insert(i);
//     edge2tets[bc].insert(i);
//     edge2tets[bd].insert(i);
//     edge2tets[cd].insert(i);
//     clow[0] = min(clow[0], circumcenters[i].x);
//     clow[1] = min(clow[1], circumcenters[i].y);
//     clow[2] = min(clow[2], circumcenters[i].z);
//     chigh[0] = max(chigh[0], circumcenters[i].x);
//     chigh[1] = max(chigh[1], circumcenters[i].y);
//     chigh[2] = max(chigh[2], circumcenters[i].z);
//   }
//   POLY_BEGIN_CONTRACT_SCOPE;
//   {
//     for (map<TetFacetHash, vector<unsigned> >::const_iterator itr = facet2tets.begin();
//          itr != facet2tets.end();
//          ++itr) POLY_ASSERT(itr->second.size() == 1 or itr->second.size() == 2);
//     for (map<EdgeHash, set<unsigned> >::const_iterator itr = edge2tets.begin();
//          itr != edge2tets.end();
//          ++itr) POLY_ASSERT(itr->second.size() >= 1);
//     POLY_ASSERT(clow[0] < chigh[0] and clow[1] < chigh[1] and clow[2] < chigh[2]);
//   }
//   POLY_END_CONTRACT_SCOPE;

//   // The bounding box for the circumcenters, useful for quantizing those pups.
//   clow[0] = min(clow[0], low[0]);
//   clow[1] = min(clow[1], low[1]);
//   clow[2] = min(clow[2], low[2]);
//   chigh[0] = max(chigh[0], high[0]); 
//   chigh[1] = max(chigh[1], high[1]); 
//   chigh[2] = max(chigh[2], high[2]);
//   RealType cboxsize = max(    chigh[0] - clow[0], 
//                           max(chigh[1] - clow[1], 
//                               chigh[2] - clow[2]));
//   const RealType cboxc[3] = {clow[0],
//                              clow[1], 
//                              clow[2]};
//   clow[0] = cboxc[0] - cboxsize;
//   clow[1] = cboxc[1] - cboxsize;
//   clow[2] = cboxc[2] - cboxsize;
//   chigh[0] = cboxc[0] + 2*cboxsize;
//   chigh[1] = cboxc[1] + 2*cboxsize;
//   chigh[2] = cboxc[2] + 2*cboxsize;
//   const RealType rinf = cboxsize;
//   cboxsize *= 3.0;

//   // Create the quantized circumcenters, and the map from the (possibly) degenerate
//   // circumcenters to their unique IDs.
//   map<PointHash, int> circ2id;
//   map<int, unsigned> tet2id;
//   for (i = 0; i != out.numberoftetrahedra; ++i) {
//     PointHash ip = geometry::Hasher<3, RealType>::hashPosition(&(circumcenters[i].x), low, high, clow, chigh);
//     j = internal::addKeyToMap(ip, circ2id);
//     tet2id[i] = j;
//   }

//   // Any surface facets create new "infinite" or "unbounded" rays, which originate at
//   // the tet circumcenter and pass through the circumcenter of the triangular facet.
//   // cerr << "Centroid bounding box: (" 
//   //      << clow[0] << " " << clow[1] << " " << clow[2] << ") ("
//   //      << chigh[0] << " " << chigh[1] << " " << chigh[2] << ")" << endl
//   //      <<    " rinf = " << rinf << endl;

//   // Look for any surface facets we need to project unbounded rays through.
//   bool test;
//   RealPoint fhat, tetcent, test_point, a_b, a_c, pinf;
//   map<TetFacetHash, unsigned> facet2id;
//   mesh.infNodes = vector<unsigned>();
//   for (map<TetFacetHash, vector<unsigned> >::const_iterator facetItr = facet2tets.begin();
//        facetItr != facet2tets.end();
//        ++facetItr) {
//     const TetFacetHash& facet = facetItr->first;
//     const vector<unsigned>& tets = facetItr->second;
//     if (tets.size() == 1) {
//       i = tets[0];
//       POLY_ASSERT(i < out.numberoftetrahedra);
//       a = out.tetrahedronlist[4*i];
//       b = out.tetrahedronlist[4*i+1];
//       c = out.tetrahedronlist[4*i+2];
//       d = out.tetrahedronlist[4*i+3];
//       POLY_ASSERT(a < numGenerators);
//       POLY_ASSERT(b < numGenerators);
//       POLY_ASSERT(c < numGenerators);
//       POLY_ASSERT(d < numGenerators);
//       geometry::computeTetCentroid(&out.pointlist[3*a],
//                                    &out.pointlist[3*b],
//                                    &out.pointlist[3*c],
//                                    &out.pointlist[3*d],
//                                    &tetcent.x);

//       // We need the ray unit vector.
//       test = geometry::computeTriangleCircumcenter3d(&out.pointlist[3*facet.x],
//                                                      &out.pointlist[3*facet.y],
//                                                      &out.pointlist[3*facet.z],
//                                                      &fhat.x);
//       POLY_ASSERT(test);
//       fhat -= circumcenters[i];

//       // Check for the special case of the tet circumcenter coplanar with the facet.
//       if (abs(geometry::dot<3, RealType>(&fhat.x, &fhat.x)) < mDegeneracy) {
//         // Yep, it's in the plane.  Just project the ray out orthogonally to the facet.
//         a_b.x = out.pointlist[3*facet.y]   - out.pointlist[3*facet.x];
//         a_b.y = out.pointlist[3*facet.y+1] - out.pointlist[3*facet.x+1];
//         a_b.z = out.pointlist[3*facet.y+2] - out.pointlist[3*facet.x+2];
//         a_c.x = out.pointlist[3*facet.z]   - out.pointlist[3*facet.x];
//         a_c.y = out.pointlist[3*facet.z+1] - out.pointlist[3*facet.x+1];
//         a_c.z = out.pointlist[3*facet.z+2] - out.pointlist[3*facet.x+2];
//         geometry::cross<3, RealType>(&a_b.x, &a_c.x, &fhat.x);
//       }
//       geometry::unitVector<3, RealType>(&fhat.x);

//       // The ray unit vector should point in the opposite direction from the facet as the tet centroid.
//       POLY_ASSERT(abs(orient3d(&out.pointlist[3*facet.x], &out.pointlist[3*facet.y], &out.pointlist[3*facet.z], &tetcent.x)) > mDegeneracy);
//       copy(&out.pointlist[3*facet.x], &out.pointlist[3*facet.x] + 3, &test_point.x);
//       test_point += fhat;
//       if (orient3d(&out.pointlist[3*facet.x], &out.pointlist[3*facet.y], &out.pointlist[3*facet.z], &tetcent.x)*
//           orient3d(&out.pointlist[3*facet.x], &out.pointlist[3*facet.y], &out.pointlist[3*facet.z], &test_point.x) > 0.0) fhat *= -1.0;

//       // Now we can compute the point where this ray intersects the surrounding "inf" sphere.
//       test = geometry::raySphereIntersection(&circumcenters[i].x,
//                                              &fhat.x,
//                                              cboxc,
//                                              rinf,
//                                              1.0e-10,
//                                              &pinf.x);
//       POLY_ASSERT(test);
//       PointHash ip = geometry::Hasher<3, RealType>::hashPosition(&(pinf.x), low, high, clow, chigh);
//       k = circ2id.size();
//       j = internal::addKeyToMap(ip, circ2id);
//       POLY_ASSERT(facet2id.find(facet) == facet2id.end());
//       facet2id[facet] = j;
//       if (k != circ2id.size()) mesh.infNodes.push_back(j);
//     }
//   }

//   // Copy the quantized nodes to the final tessellation.
//   const unsigned numNodes = circ2id.size();
//   mesh.nodes.resize(3*numNodes);
//   for (map<PointHash, int>::const_iterator itr = circ2id.begin();
//        itr != circ2id.end();
//        ++itr) {
//     POLY_ASSERT(itr->second >= 0 and itr->second < numNodes);
//     i = itr->second;
//     geometry::Hasher<3, RealType>::unhashPosition(&(mesh.nodes[3*i]), low, high, clow, chigh, itr->first);
//     // cerr << "   --> " << mesh.nodes[3*i] << " " << mesh.nodes[3*i+1] << " " << mesh.nodes[3*i+2] << endl;
//   }

//   // Build the faces corresponding to each tet edge.
//   int iface;
//   RealPoint ehat, f1, f2;
//   map<EdgeHash, int> faceMap;
//   mesh.faces.reserve(edge2tets.size());
//   mesh.cells = vector<vector<int> >(numGenerators);
//   TetFacetHash lastFacet;
//   unsigned ii, jj;
//   vector<vector<EdgeHash> > cellInfEdges(numGenerators);
//   for (map<EdgeHash, set<unsigned> >::const_iterator edgeItr = edge2tets.begin();
//        edgeItr != edge2tets.end();
//        ++edgeItr) {
//     const EdgeHash& ehash = edgeItr->first;
//     const set<unsigned>& tets = edgeItr->second;
//     a = ehash.first;
//     b = ehash.second;
//     POLY_ASSERT(a < numGenerators);
//     POLY_ASSERT(b < numGenerators);

//     set<EdgeHash> meshEdges;
//     for (set<unsigned>::const_iterator tetItr = tets.begin();
//          tetItr != tets.end();
//          ++tetItr) {
//       i = *tetItr;
//       POLY_ASSERT(i < out.numberoftetrahedra);
//       POLY_ASSERT(tet2id.find(i) != tet2id.end());
//       ii = tet2id[i];

//       // Look for edges with adjacent tets.
//       findOtherTetIndices(&out.tetrahedronlist[4*i], a, b, c, d);
//       abc = hashFacet(a, b, c);
//       abd = hashFacet(a, b, d);

//       // Is abc a surface facet?
//       if (facet2tets[abc].size() == 1) {
//         POLY_ASSERT(facet2tets[abc][0] == i);
//         POLY_ASSERT(facet2id.find(abc) != facet2id.end());
//         jj = facet2id[abc];
//         POLY_ASSERT(jj != ii);
//         meshEdges.insert(internal::hashEdge(ii, jj));
//       } else {
//         POLY_ASSERT((facet2tets[abc].size() == 2 and facet2tets[abc][0] == i) or facet2tets[abc][1] == i);
//         k = (facet2tets[abc][0] == i ? facet2tets[abc][1] : facet2tets[abc][0]);
//         jj = tet2id[k];
//         if (jj != ii) meshEdges.insert(internal::hashEdge(ii, jj));
//       }

//       // Is abd a surface facet?
//       if (facet2tets[abd].size() == 1) {
//         POLY_ASSERT(facet2tets[abd][0] == i);
//         POLY_ASSERT(facet2id.find(abd) != facet2id.end());
//         jj = facet2id[abd];
//         POLY_ASSERT(jj != ii);
//         meshEdges.insert(internal::hashEdge(ii, jj));
//       } else {
//         POLY_ASSERT(facet2tets[abd].size() == 2 and facet2tets[abd][0] == i or facet2tets[abd][1] == i);
//         k = (facet2tets[abd][0] == i ? facet2tets[abd][1] : facet2tets[abd][0]);
//         jj = tet2id[k];
//         if (jj != ii) meshEdges.insert(internal::hashEdge(ii, jj));
//       }
//     }

//     // Get the face sorted nodes.  We also keep track of which cells have "inf" nodes.
//     vector<unsigned> faceNodes, hangingNodes;
//     computeSortedFaceNodes(vector<EdgeHash>(meshEdges.begin(), meshEdges.end()), faceNodes, hangingNodes);
//     if (faceNodes.size() > 2) {
//       iface = mesh.faces.size();
//       mesh.faces.push_back(faceNodes);
//       mesh.faceCells.push_back(vector<int>());

//       // Check the orientation of the face -- we want counter-clockwise
//       // ordering viewed from outside the cell.
//       // POLY_ASSERT((not geometry::collinear<3, RealType>(&mesh.nodes[3*faceNodes[0]],
//       //                                                   &mesh.nodes[3*faceNodes[1]],
//       //                                                   &mesh.nodes[3*faceNodes[2]],
//       //                                                   mDegeneracy)));
//       ehat.x = generators[3*b]   - generators[3*a];
//       ehat.y = generators[3*b+1] - generators[3*a+1];
//       ehat.z = generators[3*b+2] - generators[3*a+2];
//       f1.x = mesh.nodes[3*faceNodes[1]]   - mesh.nodes[3*faceNodes[0]];
//       f1.y = mesh.nodes[3*faceNodes[1]+1] - mesh.nodes[3*faceNodes[0]+1];
//       f1.z = mesh.nodes[3*faceNodes[1]+2] - mesh.nodes[3*faceNodes[0]+2];
//       f2.x = mesh.nodes[3*faceNodes[2]]   - mesh.nodes[3*faceNodes[0]];
//       f2.y = mesh.nodes[3*faceNodes[2]+1] - mesh.nodes[3*faceNodes[0]+1];
//       f2.z = mesh.nodes[3*faceNodes[2]+2] - mesh.nodes[3*faceNodes[0]+2];
//       geometry::cross<3, RealType>(&f1.x, &f2.x, &fhat.x);
//       if (geometry::dot<3, RealType>(&ehat.x, &fhat.x) > 0.0) {
//         mesh.cells[a].push_back(iface);
//         mesh.cells[b].push_back(~iface);
//         mesh.faceCells[iface].push_back(a);
//         mesh.faceCells[iface].push_back(~int(b));
//       } else {
//         mesh.cells[a].push_back(~iface);
//         mesh.cells[b].push_back(iface);
//         mesh.faceCells[iface].push_back(~int(a));
//         mesh.faceCells[iface].push_back(b);
//       }
//       if (!hangingNodes.empty()) {
//         // These cells need to create a new edge on the infinite sphere.
//         POLY_ASSERT(hangingNodes.size() == 2);
//         cellInfEdges[a].push_back(internal::hashEdge(hangingNodes[0], hangingNodes[1]));
//         cellInfEdges[b].push_back(internal::hashEdge(hangingNodes[0], hangingNodes[1]));
//       }
//     }
//   }

//   // Build any infFaces we need.
//   // For now we assume there is at most one infFace per cell, which is not true for all
//   // degenerate cases!  Fix at some point...
//   for (i = 0; i != numGenerators; ++i) {
//     if (!cellInfEdges[i].empty()) {
//       POLY_ASSERT(cellInfEdges[i].size() >= 3);
//       vector<unsigned> faceNodes, hangingNodes;
//       computeSortedFaceNodes(vector<EdgeHash>(cellInfEdges[i].begin(), cellInfEdges[i].end()), faceNodes, hangingNodes);
//       POLY_ASSERT(hangingNodes.size() == 0);
//       POLY_ASSERT(faceNodes.size() >= 2);
//       if (faceNodes.size() > 2) {
//         // Check if we need to reverse the face node order.
//         ehat.x = mesh.nodes[3*faceNodes[0]]   - generators[3*i];
//         ehat.y = mesh.nodes[3*faceNodes[0]+1] - generators[3*i+1];
//         ehat.z = mesh.nodes[3*faceNodes[0]+2] - generators[3*i+2];
//         f1.x = mesh.nodes[3*faceNodes[1]]   - mesh.nodes[3*faceNodes[0]];  
//         f1.y = mesh.nodes[3*faceNodes[1]+1] - mesh.nodes[3*faceNodes[0]+1];
//         f1.z = mesh.nodes[3*faceNodes[1]+2] - mesh.nodes[3*faceNodes[0]+2];
//         f2.x = mesh.nodes[3*faceNodes[2]]   - mesh.nodes[3*faceNodes[0]];
//         f2.y = mesh.nodes[3*faceNodes[2]+1] - mesh.nodes[3*faceNodes[0]+1];
//         f2.z = mesh.nodes[3*faceNodes[2]+2] - mesh.nodes[3*faceNodes[0]+2];
//         geometry::cross<3, RealType>(&f1.x, &f2.x, &fhat.x);
//         if (geometry::dot<3, RealType>(&ehat.x, &fhat.x) < 0.0) {
//           reverse(faceNodes.begin(), faceNodes.end());
//         }
//         iface = mesh.faces.size();
//         mesh.faces.push_back(faceNodes);
//         mesh.faceCells.push_back(vector<int>(1, i));
//         mesh.infFaces.push_back(iface);
//       }
//     }
//   }

//   // Rescale the mesh node positions for the input geometry.
//   for (i = 0; i != mesh.nodes.size(); ++i) {
//     j = i % 3;
//     mesh.nodes[i] = low[j] + mesh.nodes[i]*box[j];
//   }

// }

//------------------------------------------------------------------------------
// Internal method that returns an intermediated quantized representation
// of the unbounded tessellation.
//------------------------------------------------------------------------------
void
TetgenTessellator::
computeUnboundedQuantizedTessellation(const vector<double>& points,
                                      internal::QuantTessellation<3, double>& qmesh) const {

  // Pre-conditions.
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 3 == 0);

  typedef internal::QuantTessellation<3, double>::PointHash PointHash;
  typedef internal::QuantTessellation<3, double>::EdgeHash EdgeHash;
  typedef Point3<RealType> RealPoint;
  typedef Point3<unsigned> TetFacetHash;  // kind of nefarious!

  qmesh.degeneracy = mDegeneracy;

  // Compute the normalized generators.
  const unsigned numGenerators = points.size() / 3;
  qmesh.generators = this->computeNormalizedPoints(points, points, true, &qmesh.low_labframe.x, &qmesh.high_labframe.x);
  unsigned i, j, k;

  // Build the input to tetgen.
  tetgenio in;
  in.firstnumber = 0;
  in.mesh_dim = 3;
  in.pointlist = new double[qmesh.generators.size()];
  copy(&qmesh.generators.front(), &qmesh.generators.front() + qmesh.generators.size(), in.pointlist);
  in.pointattributelist = 0;
  in.pointmtrlist = 0;
  in.pointmarkerlist = 0;
  in.numberofpoints = qmesh.generators.size() / 3;
  in.numberofpointattributes = 0;
  in.numberofpointmtrs = 0;

  // Do the tetrahedralization.
  tetgenio out;
  tetrahedralize((char*)"Q", &in, &out);
  // tetrahedralize((char*)"V", &in, &out);

  // Make sure we got something.
  if (out.numberoftetrahedra == 0)
    error("TetgenTessellator: Delauney tetrahedralization produced 0 tetrahedra!");
  if (out.numberofpoints != numGenerators) {
    char err[1024];
    snprintf(err, 1024, "TetgenTessellator: Delauney tetrahedralization produced %d tetrahedra\n(%d generating points given)", 
             out.numberofpoints, (int)numGenerators);
    error(err);
  }

  // Compute the circumcenters of the tetrahedra, and the set of tets associated
  // with each generator.
  qmesh.low_inner = RealPoint(0, 0, 0);
  qmesh.high_inner = RealPoint(1, 1, 1);
  qmesh.low_outer = RealPoint(numeric_limits<RealType>::max(),
                              numeric_limits<RealType>::max(),
                              numeric_limits<RealType>::max());
  qmesh.high_outer = RealPoint(-numeric_limits<RealType>::max(),
                               -numeric_limits<RealType>::max(),
                               -numeric_limits<RealType>::max());
  vector<RealPoint> circumcenters(out.numberoftetrahedra);
  int a, b, c, d;
  EdgeHash ab, ac, ad, bc, bd, cd;
  TetFacetHash abc, abd, bcd, acd;
  map<TetFacetHash, vector<unsigned> > facet2tets;      // Tets which share a facet.
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
    edge2tets[ab].insert(i);
    edge2tets[ac].insert(i);
    edge2tets[ad].insert(i);
    edge2tets[bc].insert(i);
    edge2tets[bd].insert(i);
    edge2tets[cd].insert(i);
    qmesh.low_outer.x = min(qmesh.low_outer.x, circumcenters[i].x);
    qmesh.low_outer.y = min(qmesh.low_outer.y, circumcenters[i].y);
    qmesh.low_outer.z = min(qmesh.low_outer.z, circumcenters[i].z);
    qmesh.high_outer.x = max(qmesh.high_outer.x, circumcenters[i].x);
    qmesh.high_outer.y = max(qmesh.high_outer.y, circumcenters[i].y);
    qmesh.high_outer.z = max(qmesh.high_outer.z, circumcenters[i].z);
  }
  POLY_BEGIN_CONTRACT_SCOPE;
  {
    for (map<TetFacetHash, vector<unsigned> >::const_iterator itr = facet2tets.begin();
         itr != facet2tets.end();
         ++itr) POLY_ASSERT(itr->second.size() == 1 or itr->second.size() == 2);
    for (map<EdgeHash, set<unsigned> >::const_iterator itr = edge2tets.begin();
         itr != edge2tets.end();
         ++itr) POLY_ASSERT(itr->second.size() >= 1);
    POLY_ASSERT(qmesh.low_outer.x < qmesh.high_outer.x and
                qmesh.low_outer.y < qmesh.high_outer.y and
                qmesh.low_outer.z < qmesh.high_outer.z);
  }
  POLY_END_CONTRACT_SCOPE;

  // Expand the outer bounding box, and choose our infSphere radius.
  qmesh.low_outer.x = min(qmesh.low_outer.x, qmesh.low_inner.x);
  qmesh.low_outer.y = min(qmesh.low_outer.y, qmesh.low_inner.y);
  qmesh.low_outer.z = min(qmesh.low_outer.z, qmesh.low_inner.z);
  qmesh.high_outer.x = max(qmesh.high_outer.x, qmesh.high_inner.x);
  qmesh.high_outer.y = max(qmesh.high_outer.y, qmesh.high_inner.y);
  qmesh.high_outer.z = max(qmesh.high_outer.z, qmesh.high_inner.z);
  RealType rinf = 4.0*max(    qmesh.high_outer.x - qmesh.low_outer.x,
                          max(qmesh.high_outer.y - qmesh.low_outer.y,
                              qmesh.high_outer.z - qmesh.low_outer.z));
  const RealPoint centroid_outer = (qmesh.low_outer + qmesh.high_outer)/2;
  qmesh.low_outer.x = centroid_outer.x - 1.05*rinf;
  qmesh.low_outer.y = centroid_outer.y - 1.05*rinf;
  qmesh.low_outer.z = centroid_outer.z - 1.05*rinf;
  qmesh.high_outer.x = centroid_outer.x + 1.05*rinf;
  qmesh.high_outer.y = centroid_outer.y + 1.05*rinf;
  qmesh.high_outer.z = centroid_outer.z + 1.05*rinf;

  // Create the quantized circumcenters, and the map from the (possibly) degenerate
  // circumcenters to their unique IDs.
  map<int, unsigned> tet2id;
  for (i = 0; i != out.numberoftetrahedra; ++i) {
    tet2id[i] = qmesh.addNewNode(circumcenters[i]);
  }

  // Any surface facets create new "infinite" or "unbounded" rays, which originate at
  // the tet circumcenter and pass through the circumcenter of the triangular facet.
  // Look for any surface facets we need to project unbounded rays through.
  bool test;
  RealPoint fhat, tetcent, test_point, a_b, a_c, pinf;
  map<TetFacetHash, unsigned> facet2id;
  qmesh.infNodes = vector<unsigned>();
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
      if (abs(geometry::dot<3, RealType>(&fhat.x, &fhat.x)) < mDegeneracy) {
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
      POLY_ASSERT(abs(orient3d(&out.pointlist[3*facet.x], &out.pointlist[3*facet.y], &out.pointlist[3*facet.z], &tetcent.x)) > mDegeneracy);
      copy(&out.pointlist[3*facet.x], &out.pointlist[3*facet.x] + 3, &test_point.x);
      test_point += fhat*rinf;
      if (orient3d(&out.pointlist[3*facet.x], &out.pointlist[3*facet.y], &out.pointlist[3*facet.z], &tetcent.x)*
          orient3d(&out.pointlist[3*facet.x], &out.pointlist[3*facet.y], &out.pointlist[3*facet.z], &test_point.x) > 0.0) fhat *= -1.0;

      // Now we can compute the point where this ray intersects the surrounding "inf" sphere.
      test = geometry::raySphereIntersection(&circumcenters[i].x,
                                             &fhat.x,
                                             &centroid_outer.x,
                                             rinf,
                                             1.0e-10,
                                             &pinf.x);
      POLY_ASSERT(test);
      
      // Add this infPoint to the quantized tessellation.
      k = qmesh.point2id.size();
      j = qmesh.addNewNode(pinf);
      POLY_ASSERT(facet2id.find(facet) == facet2id.end());
      facet2id[facet] = j;
      if (k != qmesh.point2id.size()) qmesh.infNodes.push_back(j);
    }
  }

  // Build the edges and faces corresponding to each tet edge.  Recall here that a tet edge is 
  // actualy the line connecting two generators, so not the edge of the mesh we want.
  int iedge, iface;
  RealPoint ghat, e0, e1, e2, f1, f2;
  RealType vol;
  map<EdgeHash, int> faceMap;
  qmesh.faces.reserve(edge2tets.size());
  qmesh.cells = vector<vector<int> >(numGenerators);
  TetFacetHash lastFacet;
  unsigned ii, jj;
  vector<vector<EdgeHash> > cellInfEdges(numGenerators);
  for (map<EdgeHash, set<unsigned> >::const_iterator edgeItr = edge2tets.begin();
       edgeItr != edge2tets.end();
       ++edgeItr) {
    const EdgeHash& ehash = edgeItr->first;
    const set<unsigned>& tets = edgeItr->second;
    a = ehash.first;
    b = ehash.second;
    POLY_ASSERT(a < numGenerators);
    POLY_ASSERT(b < numGenerators);

    vector<EdgeHash> meshEdges;
    for (set<unsigned>::const_iterator tetItr = tets.begin();
         tetItr != tets.end();
         ++tetItr) {
      i = *tetItr;
      POLY_ASSERT(i < out.numberoftetrahedra);
      POLY_ASSERT(tet2id.find(i) != tet2id.end());
      ii = tet2id[i];

      // Look for edges with adjacent tets.
      findOtherTetIndices(&out.tetrahedronlist[4*i], a, b, c, d);
      abc = hashFacet(a, b, c);
      abd = hashFacet(a, b, d);

      // Is abc a surface facet?
      if (facet2tets[abc].size() == 1) {
        POLY_ASSERT(facet2tets[abc][0] == i);
        POLY_ASSERT(facet2id.find(abc) != facet2id.end());
        jj = facet2id[abc];
        POLY_ASSERT(jj != ii);
        meshEdges.push_back(internal::hashEdge(ii, jj));
      } else {
        POLY_ASSERT((facet2tets[abc].size() == 2 and facet2tets[abc][0] == i) or facet2tets[abc][1] == i);
        k = (facet2tets[abc][0] == i ? facet2tets[abc][1] : facet2tets[abc][0]);
        jj = tet2id[k];
        if (jj != ii) meshEdges.push_back(internal::hashEdge(ii, jj));
      }

      // Is abd a surface facet?
      if (facet2tets[abd].size() == 1) {
        POLY_ASSERT(facet2tets[abd][0] == i);
        POLY_ASSERT(facet2id.find(abd) != facet2id.end());
        jj = facet2id[abd];
        POLY_ASSERT(jj != ii);
        meshEdges.push_back(internal::hashEdge(ii, jj));
      } else {
        POLY_ASSERT((facet2tets[abd].size() == 2) and (facet2tets[abd][0] == i or facet2tets[abd][1] == i));
        k = (facet2tets[abd][0] == i ? facet2tets[abd][1] : facet2tets[abd][0]);
        jj = tet2id[k];
        if (jj != ii) meshEdges.push_back(internal::hashEdge(ii, jj));
      }
    }

    // Arrange the edges in the correctly sorted and sign oriented order
    // to construct our face.
    sort(meshEdges.begin(), meshEdges.end());
    meshEdges.erase(unique(meshEdges.begin(), meshEdges.end()), meshEdges.end());
    if (meshEdges.size() > 1) {
      vector<int> edgeOrder;
      const bool infEdge = computeSortedFaceEdges(meshEdges, edgeOrder);
      if (meshEdges.size() > 2) {

        // Add the edges and face to the quantized mesh.
        iface = qmesh.faces.size();
        qmesh.faces.push_back(vector<int>());
        for (vector<int>::const_iterator itr = edgeOrder.begin();
             itr != edgeOrder.end();
             ++itr) {
          const bool flip = (*itr < 0);
          k = (flip ? ~(*itr) : *itr);
          iedge = qmesh.addNewEdge(meshEdges[k]);
          qmesh.nodeEdges[meshEdges[k].first].push_back(iedge);
          qmesh.nodeEdges[meshEdges[k].second].push_back(iedge);
          qmesh.edgeFaces[iedge].push_back(flip ? ~iface : iface);
          qmesh.faces[iface].push_back(flip ? ~iedge : iedge);
        }

        // Add the face to its cells.
        POLY_ASSERT(qmesh.faceCells.size() == iface);
        qmesh.faceCells.push_back(vector<int>());
        e0 = qmesh.edgePosition(meshEdges[internal::positiveID(edgeOrder[0])]);
        e1 = qmesh.edgePosition(meshEdges[internal::positiveID(edgeOrder[1])]);
        e2 = qmesh.edgePosition(meshEdges[internal::positiveID(edgeOrder[2])]);
        vol = geometry::tetrahedralVolume6(&qmesh.generators[3*a], &e2.x, &e1.x, &e0.x);
        POLY_ASSERT(vol != 0.0);
        if (vol > 0.0) {
          qmesh.cells[a].push_back(iface);
          qmesh.cells[b].push_back(~iface);
          qmesh.faceCells[iface].push_back(a);
          qmesh.faceCells[iface].push_back(~int(b));
        } else {
          qmesh.cells[a].push_back(~iface);
          qmesh.cells[b].push_back(iface);
          qmesh.faceCells[iface].push_back(~int(a));
          qmesh.faceCells[iface].push_back(b);
        }

        // Did we create a new infEdge?  If so we know it was the second one in the ordered list.
        if (infEdge) {
          j = internal::positiveID(edgeOrder[1]);
          k = qmesh.edge2id[meshEdges[j]];
          qmesh.infEdges.push_back(k);
          cellInfEdges[a].push_back(meshEdges[j]);
          cellInfEdges[b].push_back(meshEdges[j]);
        }
      }
    }
  }

  // Build any infFaces we need.
  // For now we assume there is at most one infFace per cell, which is not true for all
  // degenerate cases!  Fix at some point...
  for (i = 0; i != numGenerators; ++i) {
    if (cellInfEdges[i].size() > 2) {
      vector<int> edgeOrder;
      computeSortedFaceEdges(cellInfEdges[i], edgeOrder);

      // Check if we need to reverse the face node order.
      e0 = qmesh.edgePosition(cellInfEdges[i][internal::positiveID(edgeOrder[0])]);
      e1 = qmesh.edgePosition(cellInfEdges[i][internal::positiveID(edgeOrder[1])]);
      e2 = qmesh.edgePosition(cellInfEdges[i][internal::positiveID(edgeOrder[2])]);
      vol = geometry::tetrahedralVolume6(&qmesh.generators[3*a], &e2.x, &e1.x, &e0.x);
      POLY_ASSERT(vol != 0.0);
      if (vol < 0.0) {
        reverse(edgeOrder.begin(), edgeOrder.end());
        for (j = 0; j != edgeOrder.size(); ++j) edgeOrder[j] = ~edgeOrder[j];
      }
      iface = qmesh.faces.size();
      qmesh.faces.push_back(vector<int>());
      for (vector<int>::const_iterator itr = edgeOrder.begin();
           itr != edgeOrder.end();
           ++itr) {
        const bool flip = (*itr < 0);
        k = (flip ? ~(*itr) : *itr);
        iedge = qmesh.addNewEdge(cellInfEdges[i][k]);
        qmesh.faces[iface].push_back(flip ? ~iedge : iedge);
        qmesh.edgeFaces[iedge].push_back(flip ? ~iface : iface);
      }
      qmesh.cells[i].push_back(iface);
      qmesh.faceCells.push_back(vector<int>(1, i));
      qmesh.infFaces.push_back(iface);
    }
  }

  // Reduce the sets nodeEdges to their unique values.
  for (i = 0; i != qmesh.nodeEdges.size(); ++i) {
    sort(qmesh.nodeEdges[i].begin(), qmesh.nodeEdges[i].end());
    qmesh.nodeEdges[i].erase(unique(qmesh.nodeEdges[i].begin(), qmesh.nodeEdges[i].end()), qmesh.nodeEdges[i].end());
  }

  // Post-conditions.
  POLY_BEGIN_CONTRACT_SCOPE;
  {
    POLY_ASSERT(qmesh.points.size() == qmesh.point2id.size());
    POLY_ASSERT(qmesh.edges.size() == qmesh.edge2id.size());
    POLY_ASSERT(qmesh.cells.size() == numGenerators);
    POLY_ASSERT(qmesh.nodeEdges.size() == qmesh.point2id.size());
    POLY_ASSERT(qmesh.edgeFaces.size() == qmesh.edges.size());
    POLY_ASSERT(qmesh.faceCells.size() == qmesh.faces.size());
    for (i = 0; i != qmesh.points.size(); ++i) {
      for (j = 0; j != qmesh.nodeEdges[i].size(); ++j) {
        POLY_ASSERT(qmesh.edges[qmesh.nodeEdges[i][j]].first == i or qmesh.edges[qmesh.nodeEdges[i][j]].second == i);
      }
    }
    for (i = 0; i != qmesh.edges.size(); ++i) {
      POLY_ASSERT(qmesh.edges[i].first < qmesh.points.size());
      POLY_ASSERT(qmesh.edges[i].second < qmesh.points.size());
      POLY_ASSERT(count(qmesh.nodeEdges[qmesh.edges[i].first].begin(),
                        qmesh.nodeEdges[qmesh.edges[i].first].end(), 
                        i) == 1);
      POLY_ASSERT(count(qmesh.nodeEdges[qmesh.edges[i].second].begin(),
                        qmesh.nodeEdges[qmesh.edges[i].second].end(), 
                        i) == 1);
      for (j = 0; j != qmesh.edgeFaces[i].size(); ++j) {
        iface = qmesh.edgeFaces[i][j];
        POLY_ASSERT(internal::positiveID(iface) < qmesh.faces.size());
        if (iface < 0) {
          POLY_ASSERT(count(qmesh.faces[~iface].begin(), qmesh.faces[~iface].end(), ~i) == 1);
        } else {
          POLY_ASSERT(count(qmesh.faces[iface].begin(), qmesh.faces[iface].end(), i) == 1);
        }
      }
    }
    for (i = 0; i != qmesh.faces.size(); ++i) {
      const unsigned nedges = qmesh.faces[i].size();
      POLY_ASSERT(nedges >= 3);
      for (j = 0; j != nedges; ++j) {
        k = (j + 1) % nedges;
        const int iedge1 = qmesh.faces[i][j];
        const int iedge2 = qmesh.faces[i][k];
        POLY_ASSERT(internal::positiveID(iedge1) < qmesh.edges.size());
        POLY_ASSERT(internal::positiveID(iedge2) < qmesh.edges.size());
        if (iedge1 >= 0 and iedge2 >= 0) {
          POLY_ASSERT(qmesh.edges[iedge1].second == qmesh.edges[iedge2].first);
        } else if (iedge1 >= 0 and iedge2 < 0) {
          POLY_ASSERT(qmesh.edges[iedge1].second == qmesh.edges[~iedge2].second);
        } else if (iedge1 < 0 and iedge2 >= 0) {
          POLY_ASSERT(qmesh.edges[~iedge1].first == qmesh.edges[iedge2].first);
        } else {
          POLY_ASSERT(qmesh.edges[~iedge1].first == qmesh.edges[~iedge2].second);
        }
        if (iedge1 < 0) {
          POLY_ASSERT(count(qmesh.edgeFaces[~iedge1].begin(), qmesh.edgeFaces[~iedge1].end(), ~i) == 1);
        } else {
          POLY_ASSERT(count(qmesh.edgeFaces[iedge1].begin(), qmesh.edgeFaces[iedge1].end(), i) == 1);
        }
      }
      POLY_ASSERT(qmesh.faceCells[i].size() == 1 or qmesh.faceCells[i].size() == 2);
      for (j = 0; j != qmesh.faceCells[i].size(); ++ j) {
        const int icell = qmesh.faceCells[i][j];
        POLY_ASSERT(internal::positiveID(icell) < numGenerators);
        if (icell < 0) {
          POLY_ASSERT(count(qmesh.cells[~icell].begin(), qmesh.cells[~icell].begin(), ~i) == 1);
        } else {
          POLY_ASSERT(count(qmesh.cells[icell].begin(), qmesh.cells[icell].begin(), i) == 1);
        }
      }
    }
    for (i = 0; i != numGenerators; ++i) {
      const unsigned nfaces = qmesh.cells[i].size();
      POLY_ASSERT(nfaces >= 4);
      for (j = 0; j != nfaces; ++j) {
        iface = qmesh.cells[i][j];
        POLY_ASSERT(internal::positiveID(iface) < qmesh.faces.size());
        if (iface < 0) {
          POLY_ASSERT(count(qmesh.faceCells[~iface].begin(), qmesh.faceCells[~iface].end(), ~i) == 1);
        } else {
          POLY_ASSERT(count(qmesh.faceCells[iface].begin(), qmesh.faceCells[iface].end(), i) == 1);
        }
      }
    }
  }
  POLY_END_CONTRACT_SCOPE;
}



//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
int64_t TetgenTessellator::coordMax = (1LL << 34);
double TetgenTessellator::mDegeneracy = 1.0/TetgenTessellator::coordMax;

}
