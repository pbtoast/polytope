//------------------------------------------------------------------------------
// An internal handy intermediate representation of a tessellation.
//------------------------------------------------------------------------------
#include "polytope_internal.hh"

namespace polytope {
namespace internal {

template<int Dimension, typename RealType>
class QuantTessellation {
public:
  typedef uint64_t PointHash;
  typedef std::pair<int, int> EdgeHash;
  typedef typename DimensionTraits<Dimension, RealType>::RealPoint RealPoint;

  // The bounds for hashing positions.
  RealType low_inner[Dimension], high_inner[Dimension], 
           low_outer[Dimension], high_outer[Dimension];

  // The degeneracy we're using for quantizing.
  RealType degeneracy;

  // The mesh elements.
  std::vector<PointHash> points;             // Hashed node positions.
  std::map<PointHash, int> point2id;         // PointHash -> unique point ID
  std::map<EdgeHash, int> edge2id;           // EdgeHash  -> unique edge ID
  std::vector<std::vector<int> > faces;      // Faces made of edges (with orientation)
  std::vector<std::vector<int> > cells;      // Cells made of faces (with orientation)
  std::vector<std::vector<int> > faceCells;  // Cells of each face.
  std::vector<unsigned> infNodes;            // Indices of nodes projected to the infSphere
  std::vector<unsigned> infEdges;            // Indices of edges projected to the infSphere
  std::vector<unsigned> infFaces;            // Indices of faces projected to the infSphere

  // Hash the given position.
  PointHash hashPosition(const RealPoint& p) {
    return geometry::Hasher<Dimension, RealType>::hashPosition(const_cast<RealType*>(&(p.x)), &low_inner[0], &high_inner[0], &low_outer[0], &high_outer[0]);
  }

  // Add new elements, and return the unique index.
  int addNewNode(const RealPoint& x) {
    const PointHash ix = hashPosition(x);
    const int k = points.size();
    const int result = internal::addKeyToMap(ix, point2id);
    if (result == k) points.push_back(ix);
    POLY_ASSERT(points.size() == point2id.size());
    return result;
  }
  int addNewEdge(const EdgeHash& x) {
    return internal::addKeyToMap(x, edge2id);
  }

  // Floating position for a point.
  RealPoint nodePosition(const unsigned i) {
    POLY_ASSERT(i < points.size());
    RealPoint result;
    geometry::Hasher<Dimension, RealType>::unhashPosition(&result.x, &low_inner[0], &high_inner[0], &low_outer[0], &high_outer[0], points[i]);
    return result;
  }

  // Floating position for an edge.
  RealPoint edgePosition(const EdgeHash& ehash) {
    RealPoint result = nodePosition(ehash.first) + nodePosition(ehash.second);
    result *= 0.5;
    return result;
  }

};

}
}
