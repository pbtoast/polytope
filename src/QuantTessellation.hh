//------------------------------------------------------------------------------
// An internal handy intermediate representation of a tessellation.
//------------------------------------------------------------------------------
#include <algorithm>
#include "polytope_internal.hh"
#include "DimensionTraits.hh"

namespace polytope {
namespace internal {

template<int Dimension, typename RealType>
class QuantTessellation {
public:
  typedef uint64_t PointHash;
  typedef std::pair<int, int> EdgeHash;
  typedef typename DimensionTraits<Dimension, RealType>::CoordHash CoordHash;
  typedef typename DimensionTraits<Dimension, RealType>::Point IntPoint;
  typedef typename DimensionTraits<Dimension, RealType>::RealPoint RealPoint;

  // The normalized generator coordinates.
  std::vector<RealType> generators;

  // The bounds for hashing positions.
  RealPoint low_labframe, high_labframe;
  RealPoint low_inner, high_inner, low_outer, high_outer;

  // The degeneracy we're using for quantizing.
  RealType degeneracy;

  //----------------------------------------------------------------------------
  // The mesh elements and connectivity.
  //----------------------------------------------------------------------------
  std::map<PointHash, int> point2id;              // PointHash -> unique point ID
  std::map<EdgeHash, int> edge2id;                // EdgeHash  -> unique edge ID
  std::vector<PointHash> points;                  // Hashed node positions.
  std::vector<EdgeHash> edges;                    // Hashed edges (node index pairs).
  std::vector<std::vector<int> > faces;           // Faces made of edges (with orientation)
  std::vector<std::vector<int> > cells;           // Cells made of faces (with orientation)
  std::vector<unsigned> infNodes;                 // Indices of nodes projected to the infSphere
  std::vector<unsigned> infEdges;                 // Indices of edges projected to the infSphere
  std::vector<unsigned> infFaces;                 // Indices of faces projected to the infSphere

  //----------------------------------------------------------------------------
  // Hash the given position.
  //----------------------------------------------------------------------------
  PointHash hashPosition(const RealPoint& p) {
    return geometry::Hasher<Dimension, RealType>::hashPosition(const_cast<RealType*>(&(p.x)), 
                                                               &low_inner.x, &high_inner.x, 
                                                               &low_outer.x, &high_outer.x);
  }
  RealPoint unhashPosition(const PointHash ip) {
    RealPoint result;
    geometry::Hasher<Dimension, RealType>::unhashPosition(&result.x, 
                                                          &low_inner.x, &high_inner.x, 
                                                          &low_outer.x, &high_outer.x, 
                                                          ip);
    return result;
  }

  //----------------------------------------------------------------------------
  // Add new elements, and return the unique index.
  //----------------------------------------------------------------------------
  int addNewNode(const RealPoint& x) {
    const PointHash ix = hashPosition(x);
    const int k = point2id.size();
    const int result = internal::addKeyToMap(ix, point2id);
    if (result == k) {
      points.push_back(ix);
    }
    POLY_ASSERT(points.size() == point2id.size());
    return result;
  }

  int addNewEdge(const EdgeHash& x) {
    const int k = edge2id.size();
    const int result = internal::addKeyToMap(x, edge2id);
    if (result == k) {
      edges.push_back(x);
    }
    POLY_ASSERT(edges.size() == edge2id.size());
    return result;
  }

  //----------------------------------------------------------------------------
  // Floating position for a point (normalized coordinates).
  //----------------------------------------------------------------------------
  RealPoint nodePosition(const unsigned i) {
    POLY_ASSERT(i < points.size());
    return unhashPosition(points[i]);
  }

  //----------------------------------------------------------------------------
  // Floating position for a point (lab frame).
  //----------------------------------------------------------------------------
  RealPoint labNodePosition(const unsigned i) {
    POLY_ASSERT(i < points.size());
    RealPoint result = nodePosition(i);
    for (unsigned j = 0; j != Dimension; ++j) {
      result[j] = result[j]*(high_labframe[j] - low_labframe[j]) + low_labframe[j];
    }
    return result;
  }

  //----------------------------------------------------------------------------
  // Floating position for an edge.
  //----------------------------------------------------------------------------
  RealPoint edgePosition(const EdgeHash& ehash) {
    RealPoint result = nodePosition(ehash.first) + nodePosition(ehash.second);
    result *= 0.5;
    return result;
  }

  //----------------------------------------------------------------------------
  // Compute the node->edge connectivity.
  //----------------------------------------------------------------------------
  std::vector<std::vector<unsigned> > nodeEdges() const {
    std::vector<std::vector<unsigned> > result(points.size());
    for (unsigned iedge = 0; iedge != edges.size(); ++iedge) {
      POLY_ASSERT(edges[iedge].first < points.size());
      POLY_ASSERT(edges[iedge].second < points.size());
      result[edges[iedge].first].push_back(iedge);
      result[edges[iedge].second].push_back(iedge);
    }
    return result;
  }

  //----------------------------------------------------------------------------
  // Compute the edge->face connectivity.
  //----------------------------------------------------------------------------
  std::vector<std::vector<int> > edgeFaces() const {
    std::vector<std::vector<int> > result(edges.size());
    for (int iface = 0; iface != faces.size(); ++iface) {
      for (unsigned k = 0; k != faces[iface].size(); ++k) {
        if (faces[iface][k] < 0) {
          POLY_ASSERT(~faces[iface][k] < edges.size());
          result[~faces[iface][k]].push_back(~iface);
        } else {
          POLY_ASSERT(faces[iface][k] < edges.size());
          result[faces[iface][k]].push_back(iface);
        }
      }
    }
    return result;
  }

  //----------------------------------------------------------------------------
  // Compute the face->cell connectivity.
  //----------------------------------------------------------------------------
  std::vector<std::vector<int> > faceCells() const {
    std::vector<std::vector<int> > result(faces.size());
    for (int icell = 0; icell != cells.size(); ++icell) {
      for (unsigned k = 0; k != cells[icell].size(); ++k) {
        if (cells[icell][k] < 0) {
          POLY_ASSERT(~cells[icell][k] < faces.size());
          result[~cells[icell][k]].push_back(~icell);
        } else {
          POLY_ASSERT(cells[icell][k] < faces.size());
          result[cells[icell][k]].push_back(icell);
        }
      }
    }
    return result;
  }

  //----------------------------------------------------------------------------
  // Convert our internal data to a standard polytope Tessellation.
  //----------------------------------------------------------------------------
  void tessellation(Tessellation<Dimension, RealType>& mesh) {
    mesh.clear();

    // Nodes.
    mesh.nodes.resize(Dimension*points.size());
    for (unsigned i = 0; i != points.size(); ++i) {
      RealPoint p = labNodePosition(i);
      std::copy(&p.x, &p.x + Dimension, &mesh.nodes[Dimension*i]);
    }

    // Faces.
    mesh.faces.reserve(faces.size());
    for (unsigned i = 0; i != faces.size(); ++i) {
      mesh.faces.push_back(std::vector<unsigned>());
      for (unsigned j = 0; j != faces[i].size(); ++j) {
        if (faces[i][j] >= 0) {
          mesh.faces[i].push_back(edges[faces[i][j]].first);
        } else {
          mesh.faces[i].push_back(edges[~faces[i][j]].second);
        }
        POLY_ASSERT(mesh.faces[i].back() < mesh.nodes.size()/Dimension);
      }
      POLY_ASSERT(mesh.faces[i].size() == faces[i].size());
    }

    // Much of our data can simply be copied over wholesale.
    mesh.cells = cells;
    mesh.infNodes = infNodes;
    mesh.infFaces = infFaces;
    mesh.faceCells = this->faceCells();
  }

};

}
}