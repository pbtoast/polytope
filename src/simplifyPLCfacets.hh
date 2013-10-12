//------------------------------------------------------------------------------
// Simplify the facets of a PLC -- combine coplanar facets and remove collinear 
// points around the facets.
// We allow degenerate points between faces, but assume each input facet does
// not have degenerate points in that facet.
//------------------------------------------------------------------------------
#ifndef __Polytope_simplifyPLCfacets__
#define __Polytope_simplifyPLCfacets__

#include <vector>
#include <set>
#include <map>

#include "PLC.hh"
#include "ReducedPLC.hh"
#include "polytope_internal.hh"
#include "polytope_geometric_utilities.hh"

namespace polytope {

template<typename RealType>
ReducedPLC<3, RealType>
simplifyPLCfacets(const PLC<3, RealType>& plc,
                  const std::vector<RealType>& points,
                  const RealType tol) {
  typedef Point3<RealType> PointType;
  typedef std::pair<int, int> EdgeHash;
  ReducedPLC<3, RealType> result;
  
  // Reduce to the set of unique points.
  std::vector<RealType> newpoints;
  std::vector<unsigned> old2new;
  geometry::uniquePoints<3, RealType>(points, result.points, old2new);
  const unsigned nnewpoints = result.points.size();

  // Build a map of the old facets attached to each unique point.
  const unsigned noldfacets = plc.facets.size();
  std::vector<std::vector<unsigned> > nodeFacets(nnewpoints);
  for (unsigned i = 0; i != noldfacets; ++i) {
    for (unsigned j = 0; j != plc.facets[i].size(); ++j) {
      nodeFacets[old2new[plc.facets[i][j]]].push_back(i);
    }
  }
  for (unsigned i = 0; i != nnewpoints; ++i) {
    std::sort(nodeFacets[i].begin(), nodeFacets[i].end());
    nodeFacets[i].erase(std::unique(nodeFacets[i].begin(), nodeFacets[i].end()), nodeFacets[i].end());
  }

  // Find the normals for each starting facet.
  std::vector<PointType> facetNormals(noldfacets);
  for (unsigned i = 0; i != noldfacets; ++i) {
    POLY_ASSERT(plc.facets[i].size() >= 3);
    unsigned k1 = plc.facets[i][0], k2 = plc.facets[i][1], k3;
    POLY_ASSERT((geometry::distance<3, RealType>(&points[3*k1], &points[3*k2]) > tol));
    unsigned j = 2;
    while (j < plc.facets[i].size() and 
           geometry::collinear<3, RealType>(&points[3*k1], &points[3*k2], &points[3*plc.facets[i][j]], tol)) ++j;
    POLY_ASSERT(j < plc.facets.size());
    k3 = plc.facets[i][j];
    RealType a[3] = {points[3*k3  ] - points[3*k1  ],
                     points[3*k3+1] - points[3*k1+1],
                     points[3*k3+2] - points[3*k1+2]},
             b[3] = {points[3*k2  ] - points[3*k1  ],
                     points[3*k2+1] - points[3*k1+1],
                     points[3*k2+2] - points[3*k1+2]};
    geometry::cross<3, RealType>(a, b, &facetNormals[i].x);
    geometry::unitVector<3, RealType>(&facetNormals[i].x);
  }

  // Walk each facet of the input PLC.
  std::vector<unsigned> facetUsed(noldfacets, 0);
  for (unsigned i = 0; i != noldfacets; ++i) {
    if (facetUsed[i] == 0) {
      facetUsed[i] = 1;

      // Build the set of old facets that will contribute to this new one
      // by checking neighbor facets (by node connectivity) to see if they're
      // coplanar.
      std::set<unsigned> oldfacets, nodes2check;
      std::map<unsigned, unsigned> nodeFacetUsedCount;
      oldfacets.insert(i);
      for (std::vector<int>::const_iterator nodeItr = plc.facets[i].begin();
           nodeItr != plc.facets[i].end();
           ++nodeItr) nodes2check.insert(old2new[*nodeItr]);
      for (std::set<unsigned>::iterator nodeItr = nodes2check.begin();
           nodeItr != nodes2check.end();
           ++nodeItr) {
        for (std::vector<unsigned>::const_iterator otherFacetItr = nodeFacets[*nodeItr].begin();
             otherFacetItr != nodeFacets[*nodeItr].end();
             ++otherFacetItr) {
          if ((*otherFacetItr != i) and (facetUsed[*otherFacetItr] == 0) and
              (std::abs(geometry::dot<3, RealType>(&facetNormals[i].x, &facetNormals[*otherFacetItr].x) - 1.0) < tol)) {
            facetUsed[*otherFacetItr] = 1;
            oldfacets.insert(*otherFacetItr);
            for (std::vector<int>::const_iterator otherNodeItr = plc.facets[*otherFacetItr].begin();
                 otherNodeItr != plc.facets[*otherFacetItr].end();
                 ++otherNodeItr) nodes2check.insert(old2new[*otherNodeItr]);
          }
        }
      }

      // oldfacets now contains all the facets we're going to glue together into a single new one.
      // Walk the edges of each facet.  Any edges we walk only once is one we're keeping.
      std::vector<EdgeHash> edges;
      internal::CounterMap<EdgeHash> edgeCount;
      for (std::set<unsigned>::const_iterator facetItr = oldfacets.begin();
           facetItr != oldfacets.end();
           ++facetItr) {
        const std::vector<int>& facetNodes = plc.facets[*facetItr];
        const unsigned n = facetNodes.size();
        for (unsigned j = 0; j != n; ++j) {
          const EdgeHash ehash = internal::hashEdge(old2new[facetNodes[j]],
                                                    old2new[facetNodes[(j + 1) % n]]);
          ++edgeCount[ehash];
        }
      }
      for (typename internal::CounterMap<EdgeHash>::const_iterator itr = edgeCount.begin();
           itr != edgeCount.end();
           ++itr) {
        if (itr->second == 1) edges.push_back(itr->first);
      }
      POLY_ASSERT(edges.size() >= 3);

      // Sort the edges around the face.
      std::vector<int> edgeOrder;
      bool hangingNodes = internal::computeSortedFaceEdges(edges, edgeOrder);
      POLY_ASSERT(!hangingNodes);

      // Remove any sequential edges that are collinear.
      result.facets.push_back(std::vector<int>());
      std::vector<int>& newfacet = result.facets.back();
      for (unsigned j = 0; j != edges.size(); ++j) {
        const unsigned k = (j + 1) % edges.size();
        const int jorder = edgeOrder[j], korder = edgeOrder[k];
        const EdgeHash ej = edges[internal::positiveID(jorder)], ek = edges[internal::positiveID(korder)];
        POLY_ASSERT((jorder < 0 ? ej.first : ej.second) == (korder < 0 ? ek.second : ek.first));
        PointType jhat(result.points[3*ej.second  ] - result.points[3*ej.first  ],
                       result.points[3*ej.second+1] - result.points[3*ej.first+1],
                       result.points[3*ej.second+2] - result.points[3*ej.first+2]),
                  khat(result.points[3*ek.second  ] - result.points[3*ek.first  ],
                       result.points[3*ek.second+1] - result.points[3*ek.first+1],
                       result.points[3*ek.second+2] - result.points[3*ek.first+2]);
        jhat *= geometry::sgn(jorder);
        khat *= geometry::sgn(korder);
        geometry::unitVector<3, RealType>(&jhat.x);
        geometry::unitVector<3, RealType>(&khat.x);
        if (1.0 - std::abs(geometry::dot<3, RealType>(&jhat.x, &khat.x)) > tol) {
          newfacet.push_back(korder < 0 ? ek.second : ek.first);
        }
      }
      POLY_ASSERT(newfacet.size() >= 3);

      // Do we need to reverse the ordering of this facet?
      const PointType a0(points[3*plc.facets[i][0]], points[3*plc.facets[i][0]+1], points[3*plc.facets[i][0]+2]), 
                      b0(points[3*plc.facets[i][1]], points[3*plc.facets[i][1]+1], points[3*plc.facets[i][1]+2]), 
                      c0(points[3*plc.facets[i][2]], points[3*plc.facets[i][2]+1], points[3*plc.facets[i][2]+2]), 
                      a1(result.points[3*newfacet[0]], result.points[3*newfacet[0]+1], result.points[3*newfacet[0]+2]), 
                      b1(result.points[3*newfacet[1]], result.points[3*newfacet[1]+1], result.points[3*newfacet[1]+2]), 
                      c1(result.points[3*newfacet[2]], result.points[3*newfacet[2]+1], result.points[3*newfacet[2]+2]),
                      ab0 = b0 - a0,
                      ac0 = c0 - a0,
                      ab1 = b1 - a1,
                      ac1 = c1 - a1;
      PointType oldnorm, newnorm;
      geometry::cross<3, RealType>(&ab0.x, &ac0.x, &oldnorm.x);
      geometry::cross<3, RealType>(&ab1.x, &ac1.x, &newnorm.x);
      if (geometry::dot<3, RealType>(&oldnorm.x, &newnorm.x) < 0.0) std::reverse(newfacet.begin(), newfacet.end());
    }
  }
  POLY_ASSERT(result.facets.size() >= 4);

  // Make a final pass and remove any points we're not actually using.
  std::vector<int> mask(result.points.size()/3, -1);
  for (unsigned i = 0; i != result.facets.size(); ++i) {
    for (unsigned j = 0; j != result.facets[i].size(); ++j) {
      mask[result.facets[i][j]] = 1;
    }
  }
  unsigned n = 0;
  for (unsigned i = 0; i != mask.size(); ++i) {
    if (mask[i] == 1) mask[i] = n++;
  }
  for (unsigned i = 0; i != mask.size(); ++i) {
    if (mask[i] != -1) {
      POLY_ASSERT(mask[i] <= i);
      result.points[3*mask[i]]   = result.points[3*i];
      result.points[3*mask[i]+1] = result.points[3*i+1];
      result.points[3*mask[i]+2] = result.points[3*i+2];
    }
  }
  result.points.resize(3*n);
  for (unsigned i = 0; i != result.facets.size(); ++i) {
    for (unsigned j = 0; j != result.facets[i].size(); ++j) {
      const unsigned k = result.facets[i][j];
      POLY_ASSERT(mask[k] != -1);
      result.facets[i][j] = mask[k];
    }
  }

  // That's it.
  return result;
}

}

#endif
