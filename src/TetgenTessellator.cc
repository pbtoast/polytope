//------------------------------------------------------------------------
// TetgenTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>

#include "polytope.hh" // Pulls in POLY_ASSERT and TetgenTessellator.hh.

// Pull in tetgen stuff.
#define TETLIBRARY
#include "tetgen.cxx" 

namespace polytope {

using namespace std;

namespace {

//------------------------------------------------------------------------
// This function computes the determinant of the 3x3 matrix A with 
// the given 9 coefficients.
//------------------------------------------------------------------------
double 
det3(double A11, double A12, double A13,
     double A21, double A22, double A23,
     double A31, double A32, double A33)
{
  return A11*(A22*A33-A32*A23) - A12*(A21*A33-A31*A23) + A13*(A21*A31-A31*A22);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
// This function computes the determinant of the 4x4 matrix A with 
// the given 16 coefficients.
//------------------------------------------------------------------------
double 
det4(double A11, double A12, double A13, double A14,
     double A21, double A22, double A23, double A24,
     double A31, double A32, double A33, double A34,
     double A41, double A42, double A43, double A44)
{
  double det31 = det3(A22, A23, A24,
                      A32, A33, A34,
                      A42, A43, A44);
  double det32 = det3(A21, A23, A24,
                      A31, A33, A34,
                      A41, A43, A44);
  double det33 = det3(A21, A22, A24,
                      A31, A32, A34,
                      A41, A42, A44);
  double det34 = det3(A21, A22, A23,
                      A31, A32, A33,
                      A41, A42, A43);
  return A11*det31 - A12*det32 + A13*det33 - A14*det34;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
// This function computes the circumcenter of a tetrahedron with vertices
// A = (Ax, Ay), B = (Bx, By), C = (Cx, Cy), and D = (Dx, Dy), and places 
// the result in X.
//------------------------------------------------------------------------
void 
computeCircumcenter(double* A, double* B, double* C, double* D, double* X)
{
  // This solution was taken from 
  // http://mathworld.wolfram.com/Circumsphere.html.
  double r12 = A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
  double r22 = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
  double r32 = C[0]*C[0] + C[1]*C[1] + C[2]*C[2];
  double r42 = D[0]*D[0] + D[1]*D[1] + D[2]*D[2];
  double a = det4(A[0], A[1], A[2], 1.0,
                  B[0], B[1], B[2], 1.0,
                  C[0], C[1], C[2], 1.0,
                  D[0], D[1], D[2], 1.0);
  double Dx = det4(r12, A[1], A[2], 1.0,
                   r22, B[1], B[2], 1.0,
                   r32, C[1], C[2], 1.0,
                   r42, D[1], D[2], 1.0);
  double Dy = det4(r12, A[0], A[2], 1.0,
                   r22, B[0], B[2], 1.0,
                   r32, C[0], C[2], 1.0,
                   r42, D[0], D[2], 1.0);
  double Dz = det4(r12, A[0], A[1], 1.0,
                   r22, B[0], B[1], 1.0,
                   r32, C[0], C[1], 1.0,
                   r42, D[0], D[1], 1.0);
  X[0] = 0.5*Dx/a;
  X[1] = 0.5*Dy/a;
  X[2] = 0.5*Dz/a;
}
//------------------------------------------------------------------------

} // end anonymous namespace

//------------------------------------------------------------------------------
template<typename RealType>
TetgenTessellator<RealType>::
TetgenTessellator():
  Tessellator<3, RealType>()
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
TetgenTessellator<RealType>::
~TetgenTessellator() 
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TetgenTessellator<RealType>::
tessellate(const vector<RealType>& points,
           Tessellation<3, RealType>& mesh) const 
{
  // Create a PLC representing the bounding box for these points.
//  PLC<3, RealType> box = boundingBox(points);
  PLC<3, RealType> noBoundary;
  tessellate(points, noBoundary, mesh);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TetgenTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const PLC<3, RealType>& geometry,
           Tessellation<3, RealType>& mesh) const 
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TetgenTessellator<double>;

}
