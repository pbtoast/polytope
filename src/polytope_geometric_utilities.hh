#ifndef POLYTOPE_GEOMETRIC_UTILITIES_HH
#define POLYTOPE_GEOMETRIC_UTILITIES_HH
//------------------------------------------------------------------------------
// A semi-random collection of stuff related to geometric computations for use
// internally in polytope.
//------------------------------------------------------------------------------
#include <limits>

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
computeBoundingBox(const std::vector<RealType>& pos,
                   const bool globalReduce,
                   RealType* xmin,
                   RealType* xmax) {
   POLY_ASSERT(pos.size() % Dimension == 0);
   for (unsigned j = 0; j != Dimension; ++j) {
      xmin[j] = std::numeric_limits<RealType>::max();
      xmax[j] = (std::numeric_limits<RealType>::is_signed ? -xmin[j] : std::numeric_limits<RealType>::min());
   }
   const unsigned n = pos.size()/Dimension;
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
       xmax[j] = allReduce(xmax[j], MPI_MIN, MPI_COMM_WORLD);
     }
   }
#endif
}

}
}

#endif
