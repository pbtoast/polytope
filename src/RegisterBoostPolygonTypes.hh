//------------------------------------------------------------------------------
// This file is used to register our types with the Boost.Polygon library.
//------------------------------------------------------------------------------
#ifndef __polytope_RegisterBoostPolygonType__
#define __polytope_RegisterBoostPolygonType__

// The Voronoi tools in Boost.Polygon
#include <boost/polygon/voronoi.hpp>

// Include Boost multiprecision types
#include <boost/multiprecision/cpp_int.hpp>

#include "Point.hh"
#include "Segment.hh"

//------------------------------------------------------------------------
// Map Polytope's point class to Boost.Polygon
//------------------------------------------------------------------------
namespace boost {
namespace polygon {

typedef double fptType;
typedef int CoordHash;
typedef polytope::Point2<CoordHash> IntPoint;
typedef polytope::Segment<2> IntSegment;

template <>
struct geometry_concept<IntPoint> { typedef point_concept type; };
  
template <>
struct point_traits<IntPoint> {
  typedef IntPoint point_type;
  typedef CoordHash coordinate_type;
   
  static inline coordinate_type get(const point_type& point, orientation_2d orient) {
    return (orient == HORIZONTAL) ? point.x : point.y;
  }
};

template <>
struct point_mutable_traits<IntPoint> {
  typedef IntPoint point_type;
  typedef CoordHash coordinate_type;
   
  static inline void set(point_type& point, orientation_2d orient, coordinate_type value) {
    if (orient == HORIZONTAL)
      point.x = value;
    else
      point.y = value;
  }
  static inline point_type construct(coordinate_type x, coordinate_type y) {
    return point_type(x,y);
  }
};

//------------------------------------------------------------------------------
// We also need segments.
//------------------------------------------------------------------------------
template <>
struct geometry_concept<IntSegment> { typedef segment_concept type; };

template <>
struct segment_traits<IntSegment> {
  typedef CoordHash coordinate_type;
  typedef IntPoint point_type;

  static inline point_type get(const IntSegment& segment, direction_1d dir) {
    return dir.to_int() ? segment.b : segment.a;
  }
};

// //------------------------------------------------------------------------
// // Custom comparison operator
// //------------------------------------------------------------------------
// struct polytope_ulp_comparison {
//   enum Result {
//     LESS  = -1,
//     EQUAL =  0,
//     MORE  =  1
//   };
//   Result operator()(fptType a, fptType b, unsigned int maxUlps) const {
//     if (a > b) {
//       return a - b <= maxUlps ? EQUAL : LESS;
//     }
//     return b - a <= maxUlps ? EQUAL : MORE;
//   }
// };

// //------------------------------------------------------------------------
// // Custom floating point converter
// //------------------------------------------------------------------------
// struct polytope_fpt_converter {
//   template <typename T>
//   fptType operator()(const T& that) const {
//     return static_cast<fptType>(that);
//   }

//   template <size_t N>
//   fptType operator()(const typename detail::extended_int<N>& that) const {
//     fptType result = 0.0;
//     for (size_t i = 1; i <= (std::min)((size_t)3, that.size()); ++i) {
//       if (i != 1)
//         result *= static_cast<fptType>(0x100000000ULL);
//       result += that.chunks()[that.size() - i];
//     }
//     return (that.count() < 0) ? -result : result;
//   }
// };

// //------------------------------------------------------------------------
// // Custom voronoi diagram traits
// //------------------------------------------------------------------------
// struct polytope_voronoi_diagram_traits {
//   typedef fptType coordinate_type;
//   typedef voronoi_cell<coordinate_type> cell_type;
//   typedef voronoi_vertex<coordinate_type> vertex_type;
//   typedef voronoi_edge<coordinate_type> edge_type;
//   typedef struct {
//   public:
//     enum {ULPS = 128};
//     bool operator()(const vertex_type &v1, const vertex_type &v2) const {
//       return (ulp_cmp(v1.x(), v2.x(), ULPS) == polytope_ulp_comparison::EQUAL and
//               ulp_cmp(v1.y(), v2.y(), ULPS) == polytope_ulp_comparison::EQUAL);
//     }
//   private:
//      polytope_ulp_comparison ulp_cmp;
//   } vertex_equality_predicate_type;
// };

// //------------------------------------------------------------------------
// // Voronoi coordinate type traits
// //------------------------------------------------------------------------
// struct polytope_voronoi_ctype_traits {
//   typedef CoordHash int_type;
//   typedef detail::extended_int<3> int_x2_type;
//   typedef detail::extended_int<3> uint_x2_type;
//   typedef detail::extended_int<128> big_int_type;
//   typedef fptType fpt_type;
//   typedef fptType efpt_type;
//   typedef polytope_ulp_comparison ulp_cmp_type;
//   typedef polytope_fpt_converter to_fpt_converter_type;
//   typedef polytope_fpt_converter to_efpt_converter_type;
// };

} //end boost namespace
} //end polygon namespace


// // The Boost.Polygon Voronoi diagram
// typedef boost::polygon::voronoi_builder<polytope::DimensionTraits<2, RealType>::CoordHash, 
// 					boost::polygon::polytope_voronoi_ctype_traits> VB;
// typedef boost::polygon::voronoi_diagram<double, 
// 					boost::polygon::polytope_voronoi_diagram_traits> VD;

#endif
