//------------------------------------------------------------------------
// TetgenTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>

#include "polytope.hh" // Pulls in ASSERT and TetgenTessellator.hh.

#include "tetgen.h" 

namespace polytope {

using namespace std;

//------------------------------------------------------------------------------
template<typename Real>
TetgenTessellator<Real>::
TetgenTessellator():
  Tessellator<Real>()
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
TetgenTessellator<Real>::
~TetgenTessellator() {
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
void
TetgenTessellator<Real>::
tessellate(const vector<Real>& points,
           Tessellation<Real>& mesh) const 
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
void
TetgenTessellator<Real>::
tessellate(const vector<Real>& points,
           const PLC<Real>& geometry,
           Tessellation<Real>& mesh) const 
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TetgenTessellator<double>;

}
