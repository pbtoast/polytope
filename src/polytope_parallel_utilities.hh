#ifndef POLYTOPE_PARALLEL_UTILITIES_HH
#define POLYTOPE_PARALLEL_UTILITIES_HH
//------------------------------------------------------------------------------
// A semi-random collection of stuff that is helpful for polytope in parallel
// runs.  We just assume MPI is available here, so don't include this header
// unless that is so!
//------------------------------------------------------------------------------
#include <vector>
#include <limits>

#include "mpi.h"

#include "KeyTraits.hh"
#include "convexHull_2d.hh"
#include "convexHull_3d.hh"

namespace polytope {

//-----------------------------------------------------------------------------
// Traits to handle mapping RealType -> MPI data type.
//-----------------------------------------------------------------------------
template<typename RealType> struct DataTypeTraits;

template<> struct DataTypeTraits<unsigned> { 
  static MPI_Datatype MpiDataType() { return MPI_UNSIGNED; }
};

template<> struct DataTypeTraits<float> { 
  static MPI_Datatype MpiDataType() { return MPI_FLOAT; }
};

template<> struct DataTypeTraits<double> { 
  static MPI_Datatype MpiDataType() { return MPI_DOUBLE; }
};

//------------------------------------------------------------------------------
// A convenient interface to MPI_Allreduce.
//------------------------------------------------------------------------------
template<typename Value>
inline
Value
allReduce(const Value& value, 
          const MPI_Op op, 
          const MPI_Comm comm) {
  Value tmp = value;
  Value result;
  MPI_Allreduce(&tmp, &result, 1, DataTypeTraits<Value>::MpiDataType(), op, comm);
  return result;
}

//------------------------------------------------------------------------------
// Hide dimensional specializations.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType> struct DimensionTraits {};

// 2D
template<typename RealType>
struct DimensionTraits<2, RealType> {
  typedef typename polytope::ReducedPLC<2, RealType> ConvexHull;
  typedef polytope::KeyTraits::Key CoordHash;
  typedef polytope::Point2<CoordHash> Point;
  typedef polytope::Point2<RealType> RealPoint;

  static ConvexHull convexHull(const std::vector<RealType>& points, 
                               const RealType* low,
                               const RealType& dx) { 
    return ConvexHull(polytope::convexHull_2d(points, low, dx), points);
  }
  static Point constructPoint(const RealType* ri,
                              const RealType* rlow,
                              const RealType& dx,
                              const size_t i) {
    return Point(ri[0], ri[1], 
                 rlow[0], rlow[1], 
                 dx, i);
  }
  static RealPoint constructPoint(const RealType* ri) {
     return RealPoint(ri[0], ri[1]);
  }
  static Point faceCentroid(const polytope::Tessellation<2, RealType>& mesh,
                            const unsigned iface,
                            const RealType* rlow,
                            const RealType& dx) {
    POLY_ASSERT(iface < mesh.faces.size());
    POLY_ASSERT(mesh.faces[iface].size() == 2);
    const unsigned n1 = mesh.faces[iface][0], n2 = mesh.faces[iface][1];
    POLY_ASSERT(n1 < mesh.nodes.size()/2);
    POLY_ASSERT(n2 < mesh.nodes.size()/2);
    const Point pn1 = constructPoint(&mesh.nodes[2*n1], rlow, dx, n1);
    const Point pn2 = constructPoint(&mesh.nodes[2*n2], rlow, dx, n2);
    const Point result((pn1.x >> 1) + (pn2.x >> 1),
                       (pn1.y >> 1) + (pn2.y >> 1),
                       iface);
    return result;
  }
  static RealType maxLength(const RealType* low, const RealType* high) {
    return std::max(high[0] - low[0], high[1] - low[1]);
  }
  static std::vector<RealType> extractCoords(const std::vector<RealType>& allCoords,
                                             const std::vector<unsigned>& indices) {
    std::vector<RealType> result;
    result.reserve(2*indices.size());
    for (std::vector<unsigned>::const_iterator itr = indices.begin();
         itr != indices.end();
         ++itr) {
      const unsigned i = *itr;
      POLY_ASSERT(i < allCoords.size()/2);
      result.push_back(allCoords[2*i]);
      result.push_back(allCoords[2*i + 1]);
    }
    return result;
  }
};

// 3D
template<typename RealType>
struct DimensionTraits<3, RealType> {
  typedef typename polytope::ReducedPLC<3, RealType> ConvexHull;
  typedef polytope::KeyTraits::Key CoordHash;
  typedef polytope::Point3<CoordHash> Point;
  typedef polytope::Point3<RealType> RealPoint;

  static ConvexHull convexHull(const std::vector<RealType>& points, 
                               const RealType* low,
                               const RealType& dx) { 
    return ConvexHull(polytope::convexHull_3d(points, low, dx), points);
  }
  static Point constructPoint(const RealType* ri,
                              const RealType* rlow,
                              const RealType& dx,
                              const size_t i) {
    return Point(ri[0], ri[1], ri[2], 
                 rlow[0], rlow[1], rlow[2],
                 dx, i);
  }
  static RealPoint constructPoint(const RealType* ri) {
    return RealPoint(ri[0], ri[1], ri[2]);
  }
  static Point faceCentroid(const polytope::Tessellation<3, RealType>& mesh,
                            const unsigned iface,
                            const RealType* rlow,
                            const RealType& dx) {
    POLY_ASSERT(iface < mesh.faces.size());
    const unsigned nnodes = mesh.faces[iface].size();
    POLY_ASSERT(nnodes >= 3);
    RealType pface[3] = {0.0, 0.0, 0.0};
    unsigned ni;
    for (typename std::vector<unsigned>::const_iterator itr = mesh.faces[iface].begin();
         itr != mesh.faces[iface].end();
         ++itr) {
      ni = mesh.faces[iface][*itr];
      pface[0] += mesh.nodes[3*ni];
      pface[1] += mesh.nodes[3*ni + 1];
      pface[2] += mesh.nodes[3*ni + 2];
    }
    pface[0] /= nnodes;
    pface[1] /= nnodes;
    pface[2] /= nnodes;
    return constructPoint(pface, rlow, dx, iface);
  }
  static RealType maxLength(const RealType* low, const RealType* high) {
    return std::max(std::max(high[0] - low[0], high[1] - low[1]), high[2] - low[2]);
  }
  static std::vector<RealType> extractCoords(const std::vector<RealType>& allCoords,
                                             const std::vector<unsigned>& indices) {
    std::vector<RealType> result;
    result.reserve(3*indices.size());
    for (std::vector<unsigned>::const_iterator itr = indices.begin();
         itr != indices.end();
         ++itr) {
      const unsigned i = *itr;
      POLY_ASSERT(i < allCoords.size()/3);
      result.push_back(allCoords[3*i]);
      result.push_back(allCoords[3*i + 1]);
      result.push_back(allCoords[3*i + 2]);
    }
    return result;
  }
};

}

#endif
