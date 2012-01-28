//------------------------------------------------------------------------
// TriangleTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>

#include "polytope.hh" // Pulls in ASSERT and TriangleTessellator.hh.
#include "triangle.h" 

namespace polytope {

using namespace std;

//------------------------------------------------------------------------------
template<typename Real>
TriangleTessellator<Real>::
TriangleTessellator():
  Tessellator<Real>()
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
TriangleTessellator<Real>::
~TriangleTessellator() {
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
void
TriangleTessellator<Real>::
tessellate(const vector<Real>& points,
           Tessellation<Real>& mesh) const 
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
void
TriangleTessellator<Real>::
tessellate(const vector<Real>& points,
           const PLC<Real>& geometry,
           Tessellation<Real>& mesh) const 
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TriangleTessellator<double>;

}
