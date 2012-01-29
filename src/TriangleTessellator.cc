//------------------------------------------------------------------------
// TriangleTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>

#include "polytope.hh" // Pulls in ASSERT and TriangleTessellator.hh.

extern "C"
{
#include "triangle.h" 
}

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
~TriangleTessellator() 
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
void
TriangleTessellator<Real>::
tessellate(const vector<Real>& points,
           Tessellation<Real>& mesh) const 
{
  triangulateio in, delaunay;

  // Define input points.
  in.numberofpoints = points.size()/2;
  in.pointlist = new Real[points.size()];
  copy(points.begin(), points.end(), in.pointlist);

  // No point attributes or markers.
  in.numberofpointattributes = 0;
  in.pointattributelist = 0; 
  in.pointmarkerlist = 0; // No point markers.

  // No segments or holes.
  in.numberofsegments = 0;
  in.numberofholes = 0;

  // No regions.
  in.numberofregions = 0;
  in.regionlist = 0;

  // Set up the structure for the triangulation.
  delaunay.pointlist = 0;
  delaunay.pointattributelist = 0;
  delaunay.pointmarkerlist = 0;
  delaunay.trianglelist = 0;
  delaunay.triangleattributelist = 0;
  delaunay.neighborlist = 0;
  delaunay.edgelist = 0;
  delaunay.edgemarkerlist = 0;

  // Do the triangulation.
  triangulate((char*)"QzevnB", &in, &delaunay, 0);
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
