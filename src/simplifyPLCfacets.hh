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
#include "PLC_CSG_2d.hh"

namespace polytope {

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<2, RealType>
simplifyPLCfacets(const PLC<2, RealType>& plc,
                  const std::vector<RealType>& points,
                  const RealType* xmin,
                  const RealType* xmax,
                  const RealType tol) {
  typedef Point2<RealType> PointType;
  typedef std::pair<int, int> EdgeHash;
  ReducedPLC<2, RealType> result;
  
  // Reduce to the set of unique points.
  std::vector<unsigned> old2new;
  geometry::uniquePoints<2, RealType>(points, xmin, xmax, tol, result.points, old2new);
  const unsigned nnewpoints = result.points.size()/2;
  std::cerr << "Original points: " << std::endl;
  for (unsigned i = 0; i != points.size()/2; ++i) {
    std::cerr << "    " << i << " : (" << points[2*i] << " " << points[2*i+1] << ")" << std::endl;
  }
  std::cerr << "Unique points: " << std::endl;
  for (unsigned i = 0; i != nnewpoints; ++i) {
    std::cerr << "    " << i << " : (" << result.points[2*i] << " " << result.points[2*i+1] << ")" << std::endl;
  }

  // Sort the edges around input PLC so they're in order.
  const unsigned noldfacets = plc.facets.size();
  std::vector<EdgeHash> inputEdges(noldfacets);
  for (unsigned i = 0; i != noldfacets; ++i) {
    POLY_ASSERT(plc.facets[i].size() == 2);
    inputEdges[i] = internal::hashEdge(old2new[plc.facets[i][0]], 
                                       old2new[plc.facets[i][1]]);
  }
  std::vector<int> inputEdgeOrder;
  internal::computeSortedFaceEdges(inputEdges, inputEdgeOrder);
  POLY_ASSERT(inputEdges.size() == noldfacets);
  std::vector<EdgeHash> orderedEdges;
  for (unsigned i = 0; i != noldfacets; ++i) {
    int j = inputEdgeOrder[i];
    int k = j < 0 ? ~j : j;
    orderedEdges.push_back(j < 0 ? 
                           EdgeHash(inputEdges[k].second, inputEdges[k].first) :
                           inputEdges[k]);
  }
  POLY_ASSERT(orderedEdges.size() == noldfacets);

  // {
  //   std::cerr << "input edges : ";
  //   for (unsigned i = 0; i != noldfacets; ++i) std::cerr << " (" << inputEdges[i].first << " " << inputEdges[i].second << ")";
  //   std::cerr << std::endl
  //             << "ordered edges : ";
  //   for (unsigned i = 0; i != noldfacets; ++i) std::cerr << " (" << orderedEdges[i].first << " " << orderedEdges[i].second << ")";
  //   std::cerr << std::endl;
  // }

  // Find the first edge not collinear with previous one.
  int istart = 0;
  while (istart != noldfacets and 
         geometry::collinear<2, RealType>(&result.points[2*orderedEdges[istart].first],
                                          &result.points[2*orderedEdges[istart].second],
                                          &result.points[2*orderedEdges[(istart + noldfacets - 1) % noldfacets].first],
                                          tol)) ++istart;
  POLY_ASSERT(istart != noldfacets);

  // Now walk the edges and build up the new result, removing any sequential 
  // collinear edges.
  int nwalked = 0;
  while (nwalked < noldfacets) {
    int iedge1 = (istart + nwalked) % noldfacets, ncoll = 0;
    while (geometry::collinear<2, RealType>(&result.points[2*orderedEdges[iedge1].first],
                                            &result.points[2*orderedEdges[iedge1].second],
                                            &result.points[2*orderedEdges[(iedge1 + ncoll + 1) % noldfacets].second],
                                            tol)) ++ncoll;
    int iedge2 = (iedge1 + ncoll) % noldfacets;
    result.facets.push_back(std::vector<int>());
    result.facets.back().push_back(orderedEdges[iedge1].first);
    result.facets.back().push_back(orderedEdges[iedge2].second);
    nwalked += ncoll + 1;
  }
  POLY_ASSERT(result.facets.size() > 2);

  // Make a final pass and remove any points we're not actually using.
  std::vector<int> mask(result.points.size()/2, -1);
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
      result.points[2*mask[i]]   = result.points[2*i];
      result.points[2*mask[i]+1] = result.points[2*i+1];
    }
  }
  result.points.resize(2*n);
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

//------------------------------------------------------------------------------
// 3D
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<3, RealType>
simplifyPLCfacets(const PLC<3, RealType>& plc,
                  const std::vector<RealType>& points,
                  const RealType* xmin,
                  const RealType* xmax,
                  const RealType tol) {
  typedef Point3<RealType> PointType;
  typedef std::pair<int, int> EdgeHash;
  typedef geometry::Hasher<2, RealType> HasherType2d;
  ReducedPLC<3, RealType> result;
  
  // Reduce to the set of unique points.
  std::vector<unsigned> old2new;
  geometry::uniquePoints<3, RealType>(points, xmin, xmax, tol, result.points, old2new);
  const unsigned nnewpoints = result.points.size()/3;
  std::cerr << "Original points: " << std::endl;
  for (unsigned i = 0; i != points.size()/3; ++i) {
    std::cerr << "    " << i << " : (" << points[3*i] << " " << points[3*i+1] << " " << points[3*i+2] << ")" << std::endl;
  }
  std::cerr << "Unique points: " << std::endl;
  for (unsigned i = 0; i != nnewpoints; ++i) {
    std::cerr << "    " << i << " : (" << result.points[3*i] << " " << result.points[3*i+1] << " " << result.points[3*i+2] << ")" << std::endl;
  }

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

  // Find the normals for each starting facet.  It's possible we may have degenerate (zero area)
  // facets here, so screen those out and flag them as "used".
  std::vector<PointType> facetNormals(noldfacets);
  std::vector<double> facetNormalMag(noldfacets);
  std::vector<unsigned> facetUsed(noldfacets, 0);
  for (unsigned i = 0; i != noldfacets; ++i) {
    POLY_ASSERT(plc.facets[i].size() >= 3);
    unsigned k1 = plc.facets[i][0], k2, k3;
    unsigned j = 1;
    while (j < plc.facets[i].size() and 
           (geometry::distance<3, RealType>(&points[3*k1], &points[3*plc.facets[i][j]]) < tol)) ++j;
    if (j == plc.facets[i].size()) {
      facetUsed[i] = 1;
    } else {
      k2 = plc.facets[i][j++];
      while (j < plc.facets[i].size() and 
             geometry::collinear<3, RealType>(&points[3*k1], &points[3*k2], &points[3*plc.facets[i][j]], tol)) ++j;
      if (j == plc.facets[i].size()) {
        facetUsed[i] = 1;
      } else {
        k3 = plc.facets[i][j];
        RealType a[3] = {points[3*k3  ] - points[3*k1  ],
                         points[3*k3+1] - points[3*k1+1],
                         points[3*k3+2] - points[3*k1+2]},
                 b[3] = {points[3*k2  ] - points[3*k1  ],
                         points[3*k2+1] - points[3*k1+1],
                         points[3*k2+2] - points[3*k1+2]};
          geometry::cross<3, RealType>(b, a, &facetNormals[i].x);
          facetNormalMag[i] = sqrt(double(geometry::dot<3, RealType>(&facetNormals[i].x, &facetNormals[i].x)));
      }
    }
  }

  // Walk each (valid) facet of the input PLC.
  for (unsigned i = 0; i != noldfacets; ++i) {
    if (facetUsed[i] == 0) {
      facetUsed[i] = 1;

      // Build the set of old facets that will contribute to this new one
      // by checking neighbor facets (by node connectivity) to see if they're
      // coplanar.
      std::set<unsigned> oldfacets, nodes2check;
      std::map<unsigned, unsigned> nodeFacetUsedCount;
      oldfacets.insert(i);
      std::cerr << " gluing together old facets : " << i;
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
              (std::abs(geometry::dot<3, RealType>(&facetNormals[i].x, &facetNormals[*otherFacetItr].x) - facetNormalMag[i]*facetNormalMag[*otherFacetItr]) < tol)) {
            std::cerr << " " << *otherFacetItr;
            facetUsed[*otherFacetItr] = 1;
            oldfacets.insert(*otherFacetItr);
            for (std::vector<int>::const_iterator otherNodeItr = plc.facets[*otherFacetItr].begin();
                 otherNodeItr != plc.facets[*otherFacetItr].end();
                 ++otherNodeItr) nodes2check.insert(old2new[*otherNodeItr]);
          }
        }
      }
      std::cerr << std::endl;

      // oldfacets now contains all the facets we're going to glue together into a single new one.
      // Project to a 2D plane and use the 2D method to combine the triangles.
      bool flipFacet;
      unsigned ix, iy;
      if (std::abs(facetNormals[i].x) >= std::max(std::abs(facetNormals[i].y), std::abs(facetNormals[i].z))) {
        ix = 1; iy = 2; flipFacet = (facetNormals[i].x < 0 ? true : false);
      } else if (std::abs(facetNormals[i].y) >= std::max(std::abs(facetNormals[i].x), std::abs(facetNormals[i].z))) {
        ix = 2; iy = 0; flipFacet = (facetNormals[i].y < 0 ? true : false);
      } else {
        POLY_ASSERT(std::abs(facetNormals[i].z) >= std::max(std::abs(facetNormals[i].x), std::abs(facetNormals[i].y)));
        ix = 0; iy = 1; flipFacet = (facetNormals[i].z < 0 ? true : false);
      }
      RealType xmin2d[2] = {xmin[ix], xmin[iy]}, xmax2d[2] = {xmax[ix], xmax[iy]};
      std::map<uint64_t, int> point2dindex;
      std::vector<RealType> newpoints2d(2*nnewpoints);
      for (unsigned j = 0; j != nnewpoints; ++j) {
        newpoints2d[2*j]     = result.points[3*j + ix];
        newpoints2d[2*j + 1] = result.points[3*j + iy];
      }
      ReducedPLC<2, RealType> facet2d;
      for (std::set<unsigned>::const_iterator facetItr = oldfacets.begin();
           facetItr != oldfacets.end();
           ++facetItr) {
        ReducedPLC<2, RealType> partialfacet2d;
        const std::vector<int>& facetNodes = plc.facets[*facetItr];
        const unsigned n = facetNodes.size();
        partialfacet2d.facets.resize(n);
        for (unsigned j = 0; j != n; ++j) {
          const unsigned i1 = old2new[facetNodes[j]],
                         i2 = old2new[facetNodes[(j + 1) % n]];
          const uint64_t hash1 = HasherType2d::hashPosition(&newpoints2d[2*i1], xmin2d, xmax2d, xmin2d, xmax2d, tol),
                         hash2 = HasherType2d::hashPosition(&newpoints2d[2*i2], xmin2d, xmax2d, xmin2d, xmax2d, tol);
          point2dindex[hash1] = i1;
          point2dindex[hash2] = i2;
          if (flipFacet) {
            partialfacet2d.facets[j].push_back(i2);
            partialfacet2d.facets[j].push_back(i1);
          } else {
            partialfacet2d.facets[j].push_back(i1);
            partialfacet2d.facets[j].push_back(i2);
          }
        }
        if (flipFacet) std::reverse(partialfacet2d.facets.begin(), partialfacet2d.facets.end());
        partialfacet2d.points = newpoints2d;
        std::cerr << "  partial facet : " << partialfacet2d << std::endl;
        if (facet2d.points.size() == 0) {
          facet2d = partialfacet2d;
        } else {
          facet2d = CSG::csg_union<RealType>(facet2d, partialfacet2d);
        }
        std::cerr << "Facet so far: " << facet2d << std::endl;
      }
      facet2d = simplifyPLCfacets(facet2d, facet2d.points, xmin2d, xmax2d, tol);

      // Explode the 2d reduced facet back to the full 3d.
      std::cerr << "Exploding back to 3d." << std::endl;
      const unsigned nfacetnodes = facet2d.facets.size();
      result.facets.push_back(std::vector<int>(nfacetnodes));
      std::vector<int>& newfacet = result.facets.back();
      for (unsigned j = 0; j != nfacetnodes; ++j) {
        const unsigned k = facet2d.facets[j][0];
        std::cerr << j << " " << k << " " << " (" << facet2d.points[2*k] << " " << facet2d.points[2*k+1] << ") " << HasherType2d::hashPosition(&facet2d.points[2*k], xmin2d, xmax2d, xmin2d, xmax2d, tol) << std::endl;
        std::map<uint64_t, int>::const_iterator itr = point2dindex.find(HasherType2d::hashPosition(&facet2d.points[2*k], xmin2d, xmax2d, xmin2d, xmax2d, tol));
        POLY_ASSERT(itr != point2dindex.end());
        newfacet[j] = itr->second;
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
      if (geometry::dot<3, RealType>(&oldnorm.x, &newnorm.x) < 0) std::reverse(newfacet.begin(), newfacet.end());
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
