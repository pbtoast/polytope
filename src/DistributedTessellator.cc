//----------------------------------------------------------------------------//
// DistributedTessellator
//----------------------------------------------------------------------------//
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>

#include "polytope.hh"
#include "convexHull_2d.hh"

namespace polytope {

using namespace std;
using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
DistributedTessellator<Dimension, RealType>::
DistributedTessellator(const Tessellator<Dimension, RealType>& tessellator):
  mSerialTessellator(tessellator) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
DistributedTessellator<Dimension, RealType>::
~DistributedTessellator() {
}

//------------------------------------------------------------------------------
// Compute an unbounded tessellation.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
tessellate(const vector<RealType>& points,
           Tessellation<Dimension, RealType>& mesh) const {
  mType = unbounded;
  this->computeDistributedTessellation(points, mesh);
}

//------------------------------------------------------------------------------
// Compute the tessellation in the box.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
tessellate(const vector<RealType>& points,
           RealType* low,
           RealType* high,
           Tessellation<Dimension, RealType>& mesh) const {
  mType = box;
  mLow = low;
  mHigh = high;
  this->computeDistributedTessellation(points, mesh);
}

//------------------------------------------------------------------------------
// Compute the tessellation in a PLC.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
tessellate(const vector<RealType>& points,
           const PLC<Dimension, RealType>& geometry,
           Tessellation<Dimension, RealType>& mesh) const {
  mType = plc;
  mPLCptr = &geometry;
  this->computeDistributedTessellation(points, mesh);
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class DistributedTessellator<2, double>;

}
