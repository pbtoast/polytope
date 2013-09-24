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

//------------------------------------------------------------------------------
// Build a ReducedPLC representation of a cell.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<3, typename internal::QuantTessellation<3, RealType>::PointHash>
plcOfCell(const internal::QuantTessellation<3, RealType>& qmesh,
          const unsigned icell) {
  typedef typename internal::QuantTessellation<3, RealType>::PointHash PointHash;
  ReducedPLC<3, PointHash> result;
  std::map<int, int> old2new;
  for (unsigned i = 0; i != qmesh.cells[icell].size(); ++i) {
    result.facets.push_back(vector<int>());
    if (qmesh.cells[icell][i] < 0) {
      const unsigned iface = ~qmesh.cells[icell][i];
      const unsigned nnodes = qmesh.faces[iface].size();
      for (int j = nnodes - 1; j != -1; --j) {
        const int iedge = qmesh.faces[iface][j];
        const int ip = iedge < 0 ? qmesh.edges[~iedge].first : qmesh.edges[~iedge].second;
        if (old2new.find(ip) == old2new.end()) {
          old2new[ip] = result.points.size();
          result.points.push_back(qmesh.points[ip]);
        }
        result.facets.back().push_back(old2new[ip]);
      }
      POLY_ASSERT(result.facets.back().size() == nnodes);
    } else {
      const unsigned iface = qmesh.cells[icell][i];
      const unsigned nnodes = qmesh.faces[iface].size();
      for (int j = 0; j != nnodes; ++j) {
        const int iedge = qmesh.faces[iface][j];
        const int ip = iedge < 0 ? qmesh.edges[iedge].second : qmesh.edges[iedge].first;
        if (old2new.find(ip) == old2new.end()) {
          old2new[ip] = result.points.size();
          result.points.push_back(qmesh.points[ip]);
        }
        result.facets.back().push_back(old2new[ip]);
      }
      POLY_ASSERT(result.facets.back().size() == nnodes);
    }
  }
  POLY_ASSERT(result.facets.size() == qmesh.cells[icell].size());
  return result;
}

// //------------------------------------------------------------------------------
// // Compare a point with plane: returns (-1,0,1) for the point 
// // (below,coplanar,above) the plane
// //------------------------------------------------------------------------------
// int
// compare(const uint64_t point,
//         const uint64_t pointPlane,
//         const uint64_t normalPlane) {
//   typedef Point3<int64_t> Point;
//   typedef geometry::Hasher<3, double> HasherType;
  
//   Point p(
// }

//------------------------------------------------------------------------------
// Clip a ReducedPLC with a plane.  The plane is specified in (point, normal)
// form in the arguments.  For now we require the plane be aligned in the x, y,
// or z direction: i.e., normals (1,0,0), (0,1,0), or (0,0,1).
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<3, typename internal::QuantTessellation<3, RealType>::PointHash>
clipReducedPLC(const ReducedPLC<3, typename internal::QuantTessellation<3, RealType>::PointHash>& cell,
               const PointHash pointPlane,
               const PointHash normalPlane) {
  typedef typename internal::QuantTessellation<3, RealType>::PointHash PointHash;
  typedef geometry::Hasher<3, double> HasherType;
  const PointHash xnorm = HasherType::qxval(normalPlane),
                  ynorm = HasherType::qyval(normalPlane),
                  znorm = HasherType::qzval(normalPlane);
  POLY_ASSERT(xnorm + ynorm + znorm == 1);

  

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
  vector<double> nonGeneratingPoints;
  this->computeUnboundedQuantizedTessellation(points, nonGeneratingPoints, qmesh);

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

  typedef geometry::Hasher<3, double> HasherType;
  typedef internal::QuantTessellation<3, double>::PointHash PointHash;
  typedef internal::QuantTessellation<3, double>::EdgeHash EdgeHash;
  typedef internal::QuantTessellation<3, double>::IntPoint IntPoint;
  typedef internal::QuantTessellation<3, double>::RealPoint RealPoint;

  const PointHash outerFlag = HasherType::outerFlag();

  // Pre-conditions.
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 3 == 0);
  POLY_ASSERT(low[0] < high[0] and
              low[1] < high[1] and
              low[2] < high[2]);
  const unsigned numGenerators = points.size()/3;
  for (unsigned i = 0; i != numGenerators; ++i) {
    for (unsigned j = 0; j != 3; ++j) {
      POLY_ASSERT(low[j] <= points[3*i+j] and points[3*i+j] <= high[j]);
    }
  }

  // Create the unbounded QuantTessellation.
  internal::QuantTessellation<3, double> qmesh0;
  vector<double> nonGeneratingPoints(6);
  copy(low, low + 3, &nonGeneratingPoints[0]);
  copy(high, high + 3, &nonGeneratingPoints[3]);
  this->computeUnboundedQuantizedTessellation(points, nonGeneratingPoints, qmesh0);

  // Find the quantized box boundaries.  Both points should be in the inner bounding
  // box.
  const PointHash plow = qmesh0.hashPosition(RealPoint(low[0], low[1], low[2])),
                 phigh = qmesh0.hashPosition(RealPoint(high[0], high[1], high[2]));
  POLY_ASSERT(plow ^ outerFlag);
  POLY_ASSERT(phigh ^ outerFlag);
  const PointHash qxlow = HasherType::qxval(plow), 
                  qylow = HasherType::qyval(plow), 
                  qzlow = HasherType::qzval(plow),
                 qxhigh = HasherType::qxval(phigh), 
                 qyhigh = HasherType::qyval(phigh), 
                 qzhigh = HasherType::qzval(phigh);

  // Create a new QuantTessellation.
  internal::QuantTessellation<3, double> qmesh1;
  qmesh1.generators = qmesh0.generators;
  qmesh1.low_labframe = qmesh0.low_labframe;
  qmesh1.high_labframe = qmesh0.high_labframe;
  qmesh1.low_inner = qmesh0.low_inner;
  qmesh1.high_inner = qmesh0.high_inner;
  qmesh1.degeneracy = qmesh0.degeneracy;

  // Walk each of the cells in the unbounded tessellation.
  for (unsigned icell = 0; icell != numGenerators; ++icell) {

    // Build a PLC to represent just this cell.
    ReducedPLC<3, PointHash> cell = plcOfCell(qmesh0, icell);

    // Clip this PLC by each plane of our bounding box.
    

  }

  // Copy over all the surviving points from the unbounded mesh.
  vector<unsigned> nodes2kill(qmesh0.points.size(), 0);
  for (unsigned i = 0; i != qmesh0.points.size(); ++i) {
    if ((qmesh0.points[i] & outerFlag) or
        (HasherType::qxval(qmesh0.points[i]) < qxlow or HasherType::qxval(qmesh0.points[i]) > qxhigh) or
        (HasherType::qyval(qmesh0.points[i]) < qylow or HasherType::qyval(qmesh0.points[i]) > qyhigh) or
        (HasherType::qzval(qmesh0.points[i]) < qzlow or HasherType::qzval(qmesh0.points[i]) > qzhigh)) {
      nodes2kill[i] = 1;
    } else {
      qmesh1.point2id[qmesh0.points[i]] = qmesh1.points.size();
      qmesh1.points.push_back(qmesh0.points[i]);
    }
  }


}

//------------------------------------------------------------------------------
// Internal method that returns an intermediated quantized representation
// of the unbounded tessellation.
//------------------------------------------------------------------------------
void
TetgenTessellator::
computeUnboundedQuantizedTessellation(const vector<double>& points,
                                      const vector<double>& nonGeneratingPoints,
                                      internal::QuantTessellation<3, double>& qmesh) const {

  // Pre-conditions.
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 3 == 0);
  POLY_ASSERT(nonGeneratingPoints.size() % 3 == 0);

  typedef internal::QuantTessellation<3, double>::PointHash PointHash;
  typedef internal::QuantTessellation<3, double>::EdgeHash EdgeHash;
  typedef Point3<RealType> RealPoint;
  typedef Point3<unsigned> TetFacetHash;  // kind of nefarious!

  qmesh.degeneracy = mDegeneracy;

  // Compute the normalized generators.
  const unsigned numGenerators = points.size() / 3;
  qmesh.generators = this->computeNormalizedPoints(points, nonGeneratingPoints, true, &qmesh.low_labframe.x, &qmesh.high_labframe.x);
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
          qmesh.faces[iface].push_back(flip ? ~iedge : iedge);
        }

        // Add the face to its cells.
        e0 = qmesh.edgePosition(meshEdges[internal::positiveID(edgeOrder[0])]);
        e1 = qmesh.edgePosition(meshEdges[internal::positiveID(edgeOrder[1])]);
        e2 = qmesh.edgePosition(meshEdges[internal::positiveID(edgeOrder[2])]);
        vol = geometry::tetrahedralVolume6(&qmesh.generators[3*a], &e2.x, &e1.x, &e0.x);
        POLY_ASSERT(vol != 0.0);
        if (vol > 0.0) {
          qmesh.cells[a].push_back(iface);
          qmesh.cells[b].push_back(~iface);
        } else {
          qmesh.cells[a].push_back(~iface);
          qmesh.cells[b].push_back(iface);
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
      }
      qmesh.cells[i].push_back(iface);
      qmesh.infFaces.push_back(iface);
    }
  }

  // Post-conditions.
  POLY_BEGIN_CONTRACT_SCOPE;
  {
    const vector<vector<unsigned> > nodeEdges = qmesh.nodeEdges();
    const vector<vector<int> > edgeFaces = qmesh.edgeFaces();
    const vector<vector<int> > faceCells = qmesh.faceCells();
    POLY_ASSERT(qmesh.points.size() == qmesh.point2id.size());
    POLY_ASSERT(qmesh.edges.size() == qmesh.edge2id.size());
    POLY_ASSERT(qmesh.cells.size() == numGenerators);
    POLY_ASSERT(nodeEdges.size() == qmesh.point2id.size());
    POLY_ASSERT(edgeFaces.size() == qmesh.edges.size());
    POLY_ASSERT(faceCells.size() == qmesh.faces.size());
    for (i = 0; i != qmesh.points.size(); ++i) {
      for (j = 0; j != nodeEdges[i].size(); ++j) {
        POLY_ASSERT(qmesh.edges[nodeEdges[i][j]].first == i or qmesh.edges[nodeEdges[i][j]].second == i);
      }
    }
    for (i = 0; i != qmesh.edges.size(); ++i) {
      POLY_ASSERT(qmesh.edges[i].first < qmesh.points.size());
      POLY_ASSERT(qmesh.edges[i].second < qmesh.points.size());
      POLY_ASSERT(count(nodeEdges[qmesh.edges[i].first].begin(),
                        nodeEdges[qmesh.edges[i].first].end(), 
                        i) == 1);
      POLY_ASSERT(count(nodeEdges[qmesh.edges[i].second].begin(),
                        nodeEdges[qmesh.edges[i].second].end(), 
                        i) == 1);
      for (j = 0; j != edgeFaces[i].size(); ++j) {
        iface = edgeFaces[i][j];
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
          POLY_ASSERT(count(edgeFaces[~iedge1].begin(), edgeFaces[~iedge1].end(), ~i) == 1);
        } else {
          POLY_ASSERT(count(edgeFaces[iedge1].begin(), edgeFaces[iedge1].end(), i) == 1);
        }
      }
      POLY_ASSERT(faceCells[i].size() == 1 or faceCells[i].size() == 2);
      for (j = 0; j != faceCells[i].size(); ++ j) {
        const int icell = faceCells[i][j];
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
          POLY_ASSERT(count(faceCells[~iface].begin(), faceCells[~iface].end(), ~i) == 1);
        } else {
          POLY_ASSERT(count(faceCells[iface].begin(), faceCells[iface].end(), i) == 1);
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
