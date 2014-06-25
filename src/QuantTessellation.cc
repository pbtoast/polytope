//------------------------------------------------------------------------
// QuantTessellation
//------------------------------------------------------------------------
#include "QuantTessellation.hh"

namespace polytope {
namespace {

//----------------------------------------------------------------------------
// Clip a completed QuantTessellation to the inner bounding box.  
// After this method:
//   - All points and geometry will be inside the inner bounding box.
//   - The outer bounding box is set equal to the inner.
//   - No infinite elements.
//----------------------------------------------------------------------------
// 2D
//----------------------------------------------------------------------------
template<typename RealType>
void
QuantTessellation<2, RealType>::
clipToInnerBoundingBoxImpl(QuantTesssellation<2, RealType>& qmesh) {

  typedef QuantTessellation<2, RealType> QtessType;
  typedef QtessType::PointHash PointHash;
  typedef QtessType::IntPoint IntPoint;
  typedef QtessType::EdgeHash EdgeHash;

  // Figure out the four quantized corners of the inner box.
  const PointHash qlow_inner = hashPosition(qmesh.low_inner),
                 qhigh_inner = hashPosition(qmesh.high_inner);
  POLY_ASSERT(!(qmesh.qlow_inner  & geometry::Hasher<2, RealType>::outerFlag()));
  POLY_ASSERT(!(qmesh.qhigh_inner & geometry::Hasher<2, RealType>::outerFlag()));
  const CoordHash ixmin = geometry::Hasher<2, RealType>::qxval(qmesh.qlow_inner),
                  ixmax = geometry::Hasher<2, RealType>::qxval(qmesh.qhigh_inner),
                  iymin = geometry::Hasher<2, RealType>::qyval(qmesh.qlow_inner),
                  iymin = geometry::Hasher<2, RealType>::qyval(qmesh.qhigh_inner);
  const IntPoint ibox0 = IntPoint(ixmin, iymin),
                 ibox1 = IntPoint(ixmax, iymin),
                 ibox2 = IntPoint(ixmax, iymax),
                 ibox3 = IntPoint(ixmin, iymax);

  // Flag all the infNodes for fast checking.
  std::vector<int> infNodeFlags(qmesh.points.size(), 0);
  for (std::vector<unsigned>::iterator infItr = qmesh.infNodes.begin();
       infItr != qmesh.infNodes.end();
       ++infItr) infNodeFlags[*infItr] = 1;

  // Clip the infNodes to the inner bounding surface.
  const std::vector<std::vector<unsigned> > nodes2edges = qmesh.nodeEdges();
  std::vector<IntPoint> newPoints;
  std::vector<EdgeHash> newEdges;
  const RealPoint c1 = qmesh.low_inner, 
                  c2 = RealPoint(qmesh.high_inner.x, qmesh.low_inner.y),
                  c3 = qmesh.high_inner.x,
                  c4 = RealPoint(qmesh.low_inner.x, qmesh.high_inner.y);
  for (std::vector<unsigned>::iterator infItr = qmesh.infNodes.begin();
       infItr != qmesh.infNodes.end();
       ++infItr) {
    bool newNode = false;
    const unsigned i = *infItr;
    POLY_ASSERT(i < qmesh.nodes2edges.size());
    for (std::vector<unsigned>::const_iterator eItr = qmesh.nodes2edges[i].begin();
         eItr != qmesh.nodes2edges[i].end();
         ++eItr) {
      const unsigned ei = *eItr;
      POLY_ASSERT(edges[ei].first == i or edges[ei].second == i);

      // Does this edge have a non-infNode?  If so, it creates a new node.
      if (infNodeFlags[edges[ei].first] == 0 or infNodeFlags[edges[ei].second] == 0) {
        newPoints.push_back(hashPosition(edgeIntersectRectangle(c1, c2, c3, c4,
                                                                unhashPosition(points[edges[ei].first]),
                                                                unhashPosition(points[edges[ei].second]))));
      }
    }
  }
}

//------------------------------------------------------------------------------
// Hash an IntPoint.
// Note currently these methods should only be used for a single level 
// QuantTessellation -- no outer box!
//------------------------------------------------------------------------------
// 2D
template<typename RealType>
typename QuantTessellation<2, RealType>::PointHash

hashIntPoint(const QuantTessellation<2, RealType>& qmesh,
             const typename QuantTessellation<2, RealType>::IntPoint& p) {
  POLY_ASSERT(qmesh.low_innner == qmesh.low_outer and 
              qmesh.high_inner == qmesh.high_outer);
  return geometry::Hasher<2, RealType>::hash(p.x, p.y);
}

// 3D
template<typename RealType>
typename QuantTessellation<3, RealType>::PointHash
hashIntPoint(const QuantTessellation<3, RealType>& qmesh,
             const typename QuantTessellation<3, RealType>::IntPoint& p) {
  POLY_ASSERT(qmesh.low_innner == qmesh.low_outer and
              qmesh.high_inner == qmesh.high_outer);
  return geometry::Hasher<3, RealType>::hash(p.x, p.y, p.z);
}

}

namespace internal {

//----------------------------------------------------------------------------
// Clip the the QuantTessellation to the inner bounding box.
// We export this work to specialized 2D and 3D implementations.
//----------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
QuantTessellation<Dimension, RealType>::
clipToInnerBoundingBox() {
  clipToInnerBoundingBoxImpl(*this);
}

//------------------------------------------------------------------------------
// Hash an IntPoint.
//------------------------------------------------------------------------------
// 2D
template<typename RealType>
inline
typename QuantTessellation<2, RealType>::PointHash
QuantTessellation<2, RealType>::
hashIntPoint(const typename QuantTessellation<2, RealType>::IntPoint& p) {
  return hashIntPointImpl(*this, p);
}

}
}
