//------------------------------------------------------------------------
// QuantTessellation
//------------------------------------------------------------------------
#include "QuantTessellation.hh"

namespace polytope {
namespace {
}

namespace internal {

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
clipToInnerBoundingBox() {

  // Figure out the four quantized corners of the inner box.
  const PointHash qlow_inner = hashPosition(low_inner),
                 qhigh_inner = hashPosition(high_inner);
  POLY_ASSERT(!(qlow_inner & geometry::Hasher<2, RealType>::outerFlag()));
  POLY_ASSERT(!(qhigh_inner & geometry::Hasher<2, RealType>::outerFlag()));
  const CoordHash ixmin = geometry::Hasher<2, RealType>::qxval(qlow_inner),
                  ixmax = geometry::Hasher<2, RealType>::qxval(qhigh_inner),
                  iymin = geometry::Hasher<2, RealType>::qyval(qlow_inner),
                  iymin = geometry::Hasher<2, RealType>::qyval(qhigh_inner);
  const IntPoint ibox0 = IntPoint(ixmin, iymin),
                 ibox1 = IntPoint(ixmax, iymin),
                 ibox2 = IntPoint(ixmax, iymax),
                 ibox3 = IntPoint(ixmin, iymax);

  // Flag all the infNodes for fast checking.
  std::vector<int> infNodeFlags(points.size(), 0);
  for (std::vector<unsigned>::iterator infItr = infNodes.begin();
       infItr != infNodes.end();
       ++infItr) infNodeFlags[*infItr] = 1;

  // Clip the infNodes to the inner bounding surface.
  const std::vector<std::vector<unsigned> > nodes2edges = this->nodeEdges();
  std::vector<IntPoint> newPoints;
  std::vector<EdgeHash> newEdges;
  for (std::vector<unsigned>::iterator infItr = infNodes.begin();
       infItr != infNodes.end();
       ++infItr) {
    bool newNode = false;
    const unsigned i = *infItr;
    POLY_ASSERT(i < nodes2edges.size());
    for (std::vector<unsigned>::const_iterator eItr = nodes2edges[i].begin();
         eItr != nodes2edges[i].end();
         ++eItr) {
      const unsigned ei = *eItr;
      POLY_ASSERT(edges[ei].first == i or edges[ei].second == i);

      // Does this edge have a non-infNode?  If so, it creates a new node.
      if (infNodeFlags[edges[ei].first] == 0 or infNodeFlags[edges[ei].second] == 0) {
        newPoints.push_back(edgeIntersectBox(
      }
    }
  }
}

}
}
