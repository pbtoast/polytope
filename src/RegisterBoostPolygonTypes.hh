//------------------------------------------------------------------------------
// This file is used to register our types with the Boost.Polygon library.
//------------------------------------------------------------------------------
#ifndef __polytope_RegisterBoostPolygonType__
#define __polytope_RegisterBoostPolygonType__

// The Voronoi tools in Boost.Polygon
#include <boost/polygon/voronoi.hpp>
#include <boost/polygon/polygon.hpp>

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
typedef std::vector<IntPoint> IntPolygon;
typedef std::list<IntPolygon> IntPolygonSet;

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

//------------------------------------------------------------------------------
// Map a vector of our points to be a polygon in Boost.Polygon.
//------------------------------------------------------------------------------
template <>
struct geometry_concept<IntPolygon>{ typedef polygon_concept type; };

template <>
struct polygon_traits<IntPolygon> {
  typedef CoordHash coordinate_type;
  typedef IntPolygon::const_iterator iterator_type;
  typedef IntPoint point_type;

  // Get the begin iterator
  static inline iterator_type begin_points(const IntPolygon& t) {
    return t.begin();
  }

  // Get the end iterator
  static inline iterator_type end_points(const IntPolygon& t) {
    return t.end();
  }

  // Get the number of sides of the polygon
  static inline std::size_t size(const IntPolygon& t) {
    return t.size();
  }

  // Get the winding direction of the polygon
  static inline winding_direction winding(const IntPolygon& t) {
    return unknown_winding;
  }
};

template <>
struct polygon_mutable_traits<IntPolygon> {
  //expects stl style iterators
  template <typename iT>
  static inline IntPolygon& set_points(IntPolygon& t, 
                                       iT input_begin, iT input_end) {
    t.clear();
    t.insert(t.end(), input_begin, input_end);
    return t;
  }

};

//------------------------------------------------------------------------------
// Register our IntPolygonSet.
//------------------------------------------------------------------------------
//vector isn't automatically a polygon set in the library
//because it is a standard container there is a shortcut
//for mapping it to polygon set concept, but I'll do it
//the long way that you would use in the general case.

//first we register IntPolygonSet as a polygon set
template <>
struct geometry_concept<IntPolygonSet> { typedef polygon_set_concept type; };

//next we map to the concept through traits
template <>
struct polygon_set_traits<IntPolygonSet> {
  typedef int coordinate_type;
  typedef IntPolygonSet::const_iterator iterator_type;
  typedef IntPolygonSet operator_arg_type;

  static inline iterator_type begin(const IntPolygonSet& polygon_set) {
    return polygon_set.begin();
  }

  static inline iterator_type end(const IntPolygonSet& polygon_set) {
    return polygon_set.end();
  }

  //don't worry about these, just return false from them
  static inline bool clean(const IntPolygonSet& polygon_set) { return false; }
  static inline bool sorted(const IntPolygonSet& polygon_set) { return false; }
};

template <>
struct polygon_set_mutable_traits<IntPolygonSet> {
  template <typename input_iterator_type>
  static inline void set(IntPolygonSet& polygon_set, input_iterator_type input_begin, input_iterator_type input_end) {
    polygon_set.clear();
    //this is kind of cheesy. I am copying the unknown input geometry
    //into my own polygon set and then calling get to populate the
    //deque
    polygon_set_data<int> ps;
    ps.insert(input_begin, input_end);
    ps.get(polygon_set);
    //if you had your own odd-ball polygon set you would probably have
    //to iterate through each polygon at this point and do something
    //extra
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
