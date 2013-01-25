#ifndef POLYTOPE_GEOMETRIC_UTILITIES_HH
#define POLYTOPE_GEOMETRIC_UTILITIES_HH
//------------------------------------------------------------------------------
// A semi-random collection of stuff related to geometric computations for use
// internally in polytope.
//------------------------------------------------------------------------------
#include <limits>
#include <algorithm>

#if USE_MPI
#include "mpi.h"
#include "polytope_parallel_utilities.hh"
#endif

namespace polytope {
namespace geometry {

//------------------------------------------------------------------------------
// Distance between points.
//------------------------------------------------------------------------------
// Functor's for use in partial specializations.
template<int Dimension, typename RealType> struct DistanceFunctor;

// 2D
template<typename RealType> 
struct DistanceFunctor<2, RealType> {
  static RealType impl(const RealType* a, const RealType* b) {
    const RealType dx = a[0] - b[0];
    const RealType dy = a[1] - b[1];
    return std::sqrt(dx*dx + dy*dy);
  }
};

// 3D
template<typename RealType> 
struct DistanceFunctor<3, RealType> {
  static RealType impl(const RealType* a, const RealType* b) {
    const RealType dx = a[0] - b[0];
    const RealType dy = a[1] - b[1];
    const RealType dz = a[2] - b[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }
};

//****************************************
// Functional interface.
// This is the one you call!
//****************************************
template<int Dimension, typename RealType>
RealType
distance(const RealType* a, const RealType* b) {
  return DistanceFunctor<Dimension, RealType>::impl(a, b);
}

//------------------------------------------------------------------------------
// Compute the unit vector.
//------------------------------------------------------------------------------
// Functor's for use in partial specializations.
template<int Dimension, typename RealType> struct UnitVectorFunctor;

// 2D
template<typename RealType> 
struct UnitVectorFunctor<2, RealType> {
  static void impl(RealType* a) {
    const RealType mag = std::max(1.0e-100, std::sqrt(a[0]*a[0] + a[1]*a[1]));
    a[0] /= mag;
    a[1] /= mag;
  }
};

// 3D
template<typename RealType> 
struct UnitVectorFunctor<3, RealType> {
  static void impl(RealType* a) {
    const RealType mag = std::max(1.0e-100, std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
    a[0] /= mag;
    a[1] /= mag;
    a[2] /= mag;
  }
};

//****************************************
// Functional interface.
// This is the one you call!
//****************************************
template<int Dimension, typename RealType>
void
unitVector(RealType* a) {
  UnitVectorFunctor<Dimension, RealType>::impl(a);
}

//------------------------------------------------------------------------------
// Dot product.
//------------------------------------------------------------------------------
// Functor's for use in partial specializations.
template<int Dimension, typename RealType> struct DotFunctor;

// 2D
template<typename RealType> 
struct DotFunctor<2, RealType> {
  static RealType impl(const RealType* a, const RealType* b) {
    return a[0]*b[0] + a[1]*b[1];
  }
};

// 3D
template<typename RealType> 
struct DotFunctor<3, RealType> {
  static RealType impl(const RealType* a, const RealType* b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }
};

//****************************************
// Functional interface.
// This is the one you call!
//****************************************
template<int Dimension, typename RealType>
RealType
dot(const RealType* a, const RealType* b) {
  return DotFunctor<Dimension, RealType>::impl(a, b);
}

//------------------------------------------------------------------------------
// Cross product.
//------------------------------------------------------------------------------
// Functor's for use in partial specializations.
template<int Dimension, typename RealType> struct CrossFunctor;

// 2D
template<typename RealType> 
struct CrossFunctor<2, RealType> {
  static void impl(const RealType* a, const RealType* b, RealType* c) {
    c[0] = 0.0; c[1] = 0.0; c[2] = a[0]*b[1] - a[1]*b[0];
  }
};

// 3D
template<typename RealType> 
struct CrossFunctor<3, RealType> {
  static void impl(const RealType* a, const RealType* b, RealType* c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
  }
};

//****************************************
// Functional interface.
// This is the one you call!
//****************************************
template<int Dimension, typename RealType>
void
cross(const RealType* a, const RealType* b, RealType* c) {
  CrossFunctor<Dimension, RealType>::impl(a, b, c);
}

//------------------------------------------------------------------------------
// Determine if the given points are collinear to some accuracy.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
bool
collinear(const RealType* a, const RealType* b, const RealType* c, const RealType tol) {
  RealType ab[Dimension], ac[Dimension], abmag = 0.0, acmag = 0.0;
  for (unsigned j = 0; j != Dimension; ++j) {
    ab[j] = b[j] - a[j];
    ac[j] = c[j] - a[j];
    abmag += ab[j]*ab[j];
    acmag += ac[j]*ac[j];
  }
  if (abmag < tol or acmag < tol) return true;
  abmag = std::sqrt(abmag);
  acmag = std::sqrt(acmag);
  for (unsigned j = 0; j != Dimension; ++j) {
    ab[j] /= abmag;
    ac[j] /= acmag;
  }
  return std::abs(std::abs(dot<Dimension, RealType>(ab, ac)) - 1.0) < tol;
}

//------------------------------------------------------------------------------
// Find the closest point on a line segment (2D).
//------------------------------------------------------------------------------
template<typename RealType> 
void
closestPointOnSegment2D(const RealType* point, 
                        const RealType* s1,
                        const RealType* s2,
                        RealType* result) {
  RealType shat[2] = {s2[0] - s1[0], s2[1] - s1[1]};
  const RealType seglength = std::sqrt(shat[0]*shat[0] + shat[1]*shat[1]);
  if (seglength < 1.0e-10) {
    result[0] = 0.5*(s1[0] + s2[0]);
    result[1] = 0.5*(s1[1] + s2[1]);
  } else {
    shat[0] /= seglength;
    shat[1] /= seglength;
    const RealType s1p[2] = {point[0] - s1[0], point[1] - s1[1]};
    const RealType ptest = std::max(0.0, std::min(seglength, geometry::dot<2, RealType>(s1p, shat)));
    result[0] = s1[0] + ptest*shat[0];
    result[1] = s1[1] + ptest*shat[1];
  }
}

//------------------------------------------------------------------------------
// Find the global bounding box for a set of coordinates.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
void
computeBoundingBox(const RealType* pos,
                   const unsigned pos_size,
                   const bool globalReduce,
                   RealType xmin[Dimension],
                   RealType xmax[Dimension]) {
   POLY_ASSERT(pos_size % Dimension == 0);
   for (unsigned j = 0; j != Dimension; ++j) {
      xmin[j] = std::numeric_limits<RealType>::max();
      xmax[j] = (std::numeric_limits<RealType>::is_signed ? -xmin[j] : std::numeric_limits<RealType>::min());
   }
   const unsigned n = pos_size/Dimension;
   for (unsigned i = 0; i != n; ++i) {
      for (unsigned j = 0; j != Dimension; ++j) {
         xmin[j] = std::min(xmin[j], pos[Dimension*i + j]);
         xmax[j] = std::max(xmax[j], pos[Dimension*i + j]);
      }
   }
#if USE_MPI
   if (globalReduce) {
     for (unsigned j = 0; j != Dimension; ++j) {
       xmin[j] = allReduce(xmin[j], MPI_MIN, MPI_COMM_WORLD);
       xmax[j] = allReduce(xmax[j], MPI_MAX, MPI_COMM_WORLD);
     }
   }
#endif
}

template<int Dimension, typename RealType>
inline
void
computeBoundingBox(const std::vector<RealType>& pos,
                   const bool globalReduce,
                   RealType xmin[Dimension],
                   RealType xmax[Dimension]) {
  computeBoundingBox<Dimension, RealType>(&pos.front(),
                                          pos.size(),
                                          globalReduce,
                                          xmin,
                                          xmax);
}

//------------------------------------------------------------------------------
// Find the ray-plane intersection, return true if there's a valid answer and
// false if not.
// Arguments:
//   p_ray : point origin of ray.
//   n_ray : unit normal in direction of ray.
// p_plane : point on plane.
// n_plane : unit normal of plane.
//     tol : the tolerance for zero (check for parallel lines and such)
//  result : the intersection point (if one is possible).
//------------------------------------------------------------------------------
template<typename RealType> 
bool
rayPlaneIntersection(const RealType* p_ray,
                     const RealType* n_ray,
                     const RealType* p_plane,
                     const RealType* n_plane,
                     const RealType& tol,
                     RealType* result) {
  const RealType ndots = dot<3, RealType>(n_ray, n_plane);
  const RealType delta[3] = {p_plane[0] - p_ray[0], p_plane[1] - p_ray[1], p_plane[2] - p_ray[2]};
  if (ndots < tol) {
    // Ray is parallel to the plane.
    if (dot<3, RealType>(p_ray, delta) < tol) {
      // Ray is in the plane, just choose the origin point.
      result[0] = p_ray[0]; result[1] = p_ray[1]; result[2] = p_ray[2];
      return true;
    } else {
      // Ray is out of the plane.
      return false;
    }
  }

  // The ray is not parallel to the plane, so check for intersection.
  const RealType d = dot<3, RealType>(delta, n_plane)/ndots;
  if (d >= 0.0) {
    // Yep, they intersect.
    result[0] = p_ray[0] + d*p_ray[0];
    result[1] = p_ray[1] + d*p_ray[1];
    result[2] = p_ray[2] + d*p_ray[2];
    return true;
  } else {
    // The intersection is the wrong direction from the origin of the ray,
    // so they don't intersect.
    return false;
  }
}

//------------------------------------------------------------------------------
// Find the ray-sphere intersection, assuming the ray's origin is inside the 
// sphere.
// Arguments:
//   p_ray : point origin of ray.
//   n_ray : unit normal in direction of ray.
// p_sphere: center of the sphere.
// r_sphere: radius of the sphere.
//     tol : the tolerance for zero (check for parallel lines and such)
//  result : the intersection point
//------------------------------------------------------------------------------
template<typename RealType> 
bool
raySphereIntersection(const RealType* p_ray,
                      const RealType* n_ray,
                      const RealType* p_sphere,
                      const RealType r_sphere,
                      const RealType& tol,
                      RealType* result) {
  POLY_ASSERT2((distance<3, RealType>(p_ray, p_sphere) <= r_sphere), 
               "(" << p_ray[0] << " " << p_ray[1] << " " << p_ray[2] 
               << ") (" << p_sphere[0] << " " << p_sphere[1] << " " << p_sphere[2] << ")" << " : " 
               << (distance<3, RealType>(p_ray, p_sphere)) << " " << r_sphere);
  POLY_ASSERT(std::abs(n_ray[0]*n_ray[0] + n_ray[1]*n_ray[1] + n_ray[2]*n_ray[2] - 1.0) < 1.0e-10);
  const RealType rs0[3] = {p_ray[0] - p_sphere[0],
                           p_ray[1] - p_sphere[1],
                           p_ray[2] - p_sphere[2]};
  const RealType b = 2.0*dot<3, RealType>(n_ray, rs0);
  const RealType c = dot<3, RealType>(rs0, rs0) - r_sphere*r_sphere;
  RealType d = b*b - 4.0*c;
  POLY_ASSERT(d >= -tol);
  if (d < tol) {
    // Glancing intersection at the edge of the sphere -- it must be the origin of the ray!
    result[0] = p_ray[0]; result[1] = p_ray[1]; result[2] = p_ray[2];
    return true;
  } else {
    d = std::sqrt(d);
    const RealType t = 0.5*std::max(-b - d, -b + d);
    POLY_ASSERT(t >= 0.0);
    result[0] = p_ray[0] + t*n_ray[0];
    result[1] = p_ray[1] + t*n_ray[1];
    result[2] = p_ray[2] + t*n_ray[2];
  }
  POLY_ASSERT(std::abs(distance<3, RealType>(result, p_sphere) - r_sphere) < tol);
  return true;
}

//------------------------------------------------------------------------------
// Find the closest point of intersection for ray and a box.
// Return true if there was an intersection, false if not.
// Arguments:
//   p_ray : point origin of ray.
//   n_ray : unit normal in direction of ray.
// box_min : minimum coordinate of box.
// box_max : maximum coordinate of box.
//     tol : the tolerance for zero (check for parallel lines and such)
//  result : the intersection point (if one is possible).
//------------------------------------------------------------------------------
template<typename RealType> 
bool
rayBoxIntersection(RealType* p_ray,
                   RealType* n_ray,
                   RealType* box_min,
                   RealType* box_max,
                   const RealType& tol,
                   RealType* result) {

  // Prepare.
  bool found1 = false;
  RealType dist, minDist = std::numeric_limits<RealType>::max(), candidate[3], n_plane[3];

  // x=xmin plane.
  n_plane[0] = 1.0; n_plane[1] = 0.0; n_plane[2] = 0.0;
  if (rayPlaneIntersection(p_ray, n_ray, box_min, n_plane, tol, candidate)) {
    found1 = true;
    dist = distance<3, RealType>(p_ray, candidate);
    if (dist < minDist) {
      minDist = dist;
      std::copy(candidate, candidate + 3, result);
    }
  }

  // x=xmax plane.
  if (rayPlaneIntersection(p_ray, n_ray, box_max, n_plane, tol, candidate)) {
    found1 = true;
    dist = distance<3, RealType>(p_ray, candidate);
    if (dist < minDist) {
      minDist = dist;
      std::copy(candidate, candidate + 3, result);
    }
  }

  // y=ymin plane.
  n_plane[0] = 0.0; n_plane[1] = 1.0; n_plane[2] = 0.0;
  if (rayPlaneIntersection(p_ray, n_ray, box_min, n_plane, tol, candidate)) {
    found1 = true;
    dist = distance<3, RealType>(p_ray, candidate);
    if (dist < minDist) {
      minDist = dist;
      std::copy(candidate, candidate + 3, result);
    }
  }

  // y=ymax plane.
  if (rayPlaneIntersection(p_ray, n_ray, box_max, n_plane, tol, candidate)) {
    found1 = true;
    dist = distance<3, RealType>(p_ray, candidate);
    if (dist < minDist) {
      minDist = dist;
      std::copy(candidate, candidate + 3, result);
    }
  }

  // z=zmin plane.
  n_plane[0] = 0.0; n_plane[1] = 0.0; n_plane[2] = 1.0;
  if (rayPlaneIntersection(p_ray, n_ray, box_min, n_plane, tol, candidate)) {
    found1 = true;
    dist = distance<3, RealType>(p_ray, candidate);
    if (dist < minDist) {
      minDist = dist;
      std::copy(candidate, candidate + 3, result);
    }
  }

  // z=zmax plane.
  if (rayPlaneIntersection(p_ray, n_ray, box_max, n_plane, tol, candidate)) {
    found1 = true;
    dist = distance<3, RealType>(p_ray, candidate);
    if (dist < minDist) {
      minDist = dist;
      std::copy(candidate, candidate + 3, result);
    }
  }

  return found1;
}

//------------------------------------------------------------------------------
// Compute the centroid of the mesh cell.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
computeCellCentroid(const Tessellation<Dimension, RealType>& mesh,
                    const unsigned ci,
                    RealType* ccent) {
  POLY_ASSERT(ci < mesh.cells.size());
  const unsigned nf = mesh.cells[ci].size();
  unsigned i, j, k;
  for (j = 0; j != Dimension; ++j) ccent[j] = 0.0;
  std::vector<unsigned> uniqueNodes;
  for (k = 0; k != nf; ++k) {
    const std::vector<unsigned>& faceNodes = mesh.faces[internal::positiveID(mesh.cells[ci][k])];
    std::copy(faceNodes.begin(), faceNodes.end(), std::back_inserter(uniqueNodes));
  }

  // Now reduce to the unique vertices.
  std::sort(uniqueNodes.begin(), uniqueNodes.end());
  std::vector<unsigned>::iterator endItr = std::unique(uniqueNodes.begin(), uniqueNodes.end());
  for (std::vector<unsigned>::iterator itr = uniqueNodes.begin(); itr != endItr; ++itr) {
    i = *itr;
    POLY_ASSERT2(i < mesh.nodes.size()/Dimension, i << " " << mesh.nodes.size()/Dimension);
    for (j = 0; j != Dimension; ++j) ccent[j] += mesh.nodes[Dimension*i + j];
  }
  const unsigned n = std::distance(uniqueNodes.begin(), endItr);
  POLY_ASSERT(n > 0);
  for (j = 0; j != Dimension; ++j) ccent[j] /= n;
}

//------------------------------------------------------------------------------
// Compute the centroid and unit normal of a Tessellation face.
//------------------------------------------------------------------------------
template<typename RealType>
void
computeFaceCentroidAndNormal(const Tessellation<3, RealType>& mesh,
                             const unsigned fi,
                             RealType* fcent,
                             RealType* fhat) {
  POLY_ASSERT(fi < mesh.faces.size());
  const unsigned n = mesh.faces[fi].size();
  POLY_ASSERT(n >= 3);
  unsigned i, j, ni, nj;
  std::vector<unsigned> verts;
  const double degeneracy = 1.0e-10;

  // Compute the centroid, and look for three vertices in the face that are 
  // not collinear.
  fcent[0] = 0.0; fcent[1] = 0.0; fcent[2] = 0.0;
  for (i = 0; i != n; ++i) {
    ni = mesh.faces[fi][i];
    POLY_ASSERT(ni < mesh.nodes.size()/3);
    fcent[0] += mesh.nodes[3*ni];
    fcent[1] += mesh.nodes[3*ni+1];
    fcent[2] += mesh.nodes[3*ni+2];
    if (verts.size() < 2 or
        (verts.size() == 2 and not collinear<3, RealType>(&mesh.nodes[3*verts[0]],
                                                          &mesh.nodes[3*verts[1]],
                                                          &mesh.nodes[3*ni],
                                                          degeneracy))) verts.push_back(ni);
  }
  POLY_ASSERT(n > 0);
  fcent[0] /= n; fcent[1] /= n; fcent[2] /= n;
  
  // Now we can compute the unit normal.
  POLY_ASSERT(verts.size() == 3);
  RealType ab[3], ac[3];
  ab[0] = mesh.nodes[3*verts[1]  ] - mesh.nodes[3*verts[0]  ];
  ab[1] = mesh.nodes[3*verts[1]+1] - mesh.nodes[3*verts[0]+1];
  ab[2] = mesh.nodes[3*verts[1]+2] - mesh.nodes[3*verts[0]+2];
  ac[0] = mesh.nodes[3*verts[2]  ] - mesh.nodes[3*verts[0]  ];
  ac[1] = mesh.nodes[3*verts[2]+1] - mesh.nodes[3*verts[0]+1];
  ac[2] = mesh.nodes[3*verts[2]+2] - mesh.nodes[3*verts[0]+2];
  cross<3, RealType>(ab, ac, fhat);
  unitVector<3, RealType>(fhat);
}

//------------------------------------------------------------------------------
// This function computes the circumcenter of a triangle with vertices
// A = (Ax, Ay), B = (Bx, By), and C = (Cx, Cy), and places the result 
// in X.
//------------------------------------------------------------------------------
inline
void 
computeCircumcenter2d(double* A, double* B, double* C, double* X) {
  // This solution was taken from Wikipedia's entry:
  // http://en.wikipedia.org/wiki/Circumscribed_circle
  double D = 2.0*(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1]));
  X[0] = ((A[0]*A[0] + A[1]*A[1])*(B[1]-C[1]) + (B[0]*B[0] + B[1]*B[1])*(C[1]-A[1]) + 
          (C[0]*C[0] + C[1]*C[1])*(A[1]-B[1]))/D;
  X[1] = ((A[0]*A[0] + A[1]*A[1])*(C[0]-B[0]) + (B[0]*B[0] + B[1]*B[1])*(A[0]-C[0]) + 
          (C[0]*C[0] + C[1]*C[1])*(B[0]-A[0]))/D;
}

//------------------------------------------------------------------------
// This function computes the determinant of the 3x3 matrix A with 
// the given 9 coefficients.
//------------------------------------------------------------------------
inline
double 
det3(double A11, double A12, double A13,
     double A21, double A22, double A23,
     double A31, double A32, double A33) {
  return A11*(A22*A33-A32*A23) - A12*(A21*A33-A31*A23) + A13*(A21*A31-A31*A22);
}

//------------------------------------------------------------------------
// This function computes the determinant of the 4x4 matrix A with 
// the given 16 coefficients.
//------------------------------------------------------------------------
inline
double 
det4(double A11, double A12, double A13, double A14,
     double A21, double A22, double A23, double A24,
     double A31, double A32, double A33, double A34,
     double A41, double A42, double A43, double A44) {
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
// This function computes the circumcenter of a tetrahedron with vertices
// A = (Ax, Ay, Az), B = (Bx, By, Bz), C = (Cx, Cy, Cz), and D = (Dx, Dy, Dz).
// The result is returned in X.
//------------------------------------------------------------------------
inline
void
computeCircumcenter3d(double* A, double* B, double* C, double* D, double* X) {
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
  POLY_ASSERT(a != 0.0);
  X[0] = 0.5*Dx/a;
  X[1] = 0.5*Dy/a;
  X[2] = 0.5*Dz/a;
}

//------------------------------------------------------------------------------
// Compute the circumcenter of a triangle in 3D.
// Same argument conventions as above.
// Based on the form presented at
// http://en.wikipedia.org/wiki/Circumscribed_circle
//------------------------------------------------------------------------------
inline
bool
computeTriangleCircumcenter3d(double* A, double* B, double* C, double* X) {
  const double degeneracy = 1.0e-10;
  if (collinear<3, double>(A, B, C, degeneracy)) return false;
  
  const double a[3] = {A[0] - C[0], 
                       A[1] - C[1], 
                       A[2] - C[2]};
  const double b[3] = {B[0] - C[0],
                       B[1] - C[1],
                       B[2] - C[2]};
  const double
    a2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2],
    b2 = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
  const double c[3] = {a2*b[0] - b2*a[0],
                       a2*b[1] - b2*a[1],
                       a2*b[2] - b2*a[2]};
  double d[3], e[3];
  cross<3, double>(a, b, d);
  cross<3, double>(c, d, e);
  const double d2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
  POLY_ASSERT(d2 > 0.0);
  X[0] = 0.5*e[0]/d2 + C[0];
  X[1] = 0.5*e[1]/d2 + C[1];
  X[2] = 0.5*e[2]/d2 + C[2];
  return true;
}

}
}

#endif
