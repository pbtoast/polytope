//------------------------------------------------------------------------
// Triangle
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>

#include "polytope.hh" // Pulls in ASSERT and Triangle.hh.
//#include "triangle.h" 

namespace polytope {

using namespace std;

//------------------------------------------------------------------------------
template<typename Real>
Triangle<Real>::
Triangle():
  Tessellator<Real>()
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
Triangle<Real>::
~Triangle() {
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
void
Triangle<Real>::
tessellate(const vector<Real>& points,
           Tessellation<Real>& mesh) const 
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
void
Triangle<Real>::
tessellate(const vector<Real>& points,
           const PLC<Real>& geometry,
           Tessellation<Real>& mesh) const 
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class Triangle<double>;

}
