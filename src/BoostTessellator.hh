//------------------------------------------------------------------------
// BoostTessellator
// 
// Polytope wrapper for the native 2D Voronoi tessellator in Boost.Polygon
// v1.52 or greater
//------------------------------------------------------------------------
#ifndef __Polytope_BoostTessellator__
#define __Polytope_BoostTessellator__

#ifdef HAVE_BOOST_VORONOI

#include <vector>
#include <cmath>

#include "Tessellator.hh"
#include "Clipper2d.hh"
#include "BoostOrphanage.hh"
#include "Point.hh"
#include "polytope_tessellator_utilities.hh"
#include "BoostTessellatorTraits.hh"

// The Voronoi tools in Boost.Polygon
#include <boost/polygon/voronoi.hpp>

// We use the Boost.Geometry library to handle polygon intersections and such.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

// Include Boost multiprecision types
#include <boost/multiprecision/cpp_int.hpp>

//------------------------------------------------------------------------
// Map Polytope's point class to Boost.Polygon
//------------------------------------------------------------------------
namespace boost{
namespace polygon{

typedef double fptType;
typedef polytope::DimensionTraits<2, double>::CoordHash CoordHash;
typedef polytope::Point2<CoordHash> IntPoint;

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

//------------------------------------------------------------------------
// Custom comparison operator
//------------------------------------------------------------------------
struct polytope_ulp_comparison {
  enum Result {
    LESS  = -1,
    EQUAL =  0,
    MORE  =  1
  };
  Result operator()(fptType a, fptType b, unsigned int maxUlps) const {
    if (a > b) {
      return a - b <= maxUlps ? EQUAL : LESS;
    }
    return b - a <= maxUlps ? EQUAL : MORE;
  }
};

//------------------------------------------------------------------------
// Custom floating point converter
//------------------------------------------------------------------------
struct polytope_fpt_converter {
  template <typename T>
  fptType operator()(const T& that) const {
    return static_cast<fptType>(that);
  }

  template <size_t N>
  fptType operator()(const typename detail::extended_int<N>& that) const {
    fptType result = 0.0;
    for (size_t i = 1; i <= (std::min)((size_t)3, that.size()); ++i) {
      if (i != 1)
        result *= static_cast<fptType>(0x100000000ULL);
      result += that.chunks()[that.size() - i];
    }
    return (that.count() < 0) ? -result : result;
  }
};

//------------------------------------------------------------------------
// Custom voronoi diagram traits
//------------------------------------------------------------------------
struct polytope_voronoi_diagram_traits {
  typedef fptType coordinate_type;
  typedef voronoi_cell<coordinate_type> cell_type;
  typedef voronoi_vertex<coordinate_type> vertex_type;
  typedef voronoi_edge<coordinate_type> edge_type;
  typedef struct {
  public:
    enum {ULPS = 128};
    bool operator()(const vertex_type &v1, const vertex_type &v2) const {
      return (ulp_cmp(v1.x(), v2.x(), ULPS) == polytope_ulp_comparison::EQUAL and
              ulp_cmp(v1.y(), v2.y(), ULPS) == polytope_ulp_comparison::EQUAL);
    }
  private:
     polytope_ulp_comparison ulp_cmp;
  } vertex_equality_predicate_type;
};

//------------------------------------------------------------------------
// Voronoi coordinate type traits
//------------------------------------------------------------------------
struct polytope_voronoi_ctype_traits {
  typedef CoordHash int_type;
  typedef detail::extended_int<3> int_x2_type;
  typedef detail::extended_int<3> uint_x2_type;
  typedef detail::extended_int<128> big_int_type;
  typedef fptType fpt_type;
  typedef fptType efpt_type;
  typedef polytope_ulp_comparison ulp_cmp_type;
  typedef polytope_fpt_converter to_fpt_converter_type;
  typedef polytope_fpt_converter to_efpt_converter_type;
};

} //end boost namespace
} //end polygon namespace


// // The Boost.Polygon Voronoi diagram
// typedef boost::polygon::voronoi_builder<polytope::DimensionTraits<2, RealType>::CoordHash, 
// 					boost::polygon::polytope_voronoi_ctype_traits> VB;
// typedef boost::polygon::voronoi_diagram<double, 
// 					boost::polygon::polytope_voronoi_diagram_traits> VD;

namespace polytope {

template<typename RealType>
class BoostTessellator: public Tessellator<2, RealType> {
public:

  // The Boost.Polygon Voronoi diagram
  typedef boost::polygon::voronoi_diagram<RealType> VD;

  // Some useful typedefs
  typedef std::pair<int, int> EdgeHash;
  typedef typename polytope::DimensionTraits<2, RealType>::CoordHash CoordHash;

  // The Big Switch: geometric operations in integers or doubles?
  typedef BoostTessellatorTraits<RealType, RealType>  BTT;     // uncomment to use floating point
  // typedef BoostTessellatorTraits<RealType, CoordHash> BTT;  // uncomment to use integers
   
  // The typedefs that follow from this choice
  typedef typename BTT::RealPoint RealPoint;
  typedef typename BTT::IntPoint  IntPoint;
  typedef typename BTT::PointType PointType;
  typedef typename BTT::CoordType CoordType;

  // Constructor, destructor.
  BoostTessellator();
  ~BoostTessellator();

  // Tessellate the given generators. A bounding box is constructed about
  // the generators, and the corners of the bounding box are added as 
  // additional generators if they are not present in the list.
  void tessellate(const std::vector<RealType>& points,
                  Tessellation<2, RealType>& mesh) const;

  // Tessellate with a bounding box representing the boundaries.
  void tessellate(const std::vector<RealType>& points,
                  RealType* low,
                  RealType* high,
                  Tessellation<2, RealType>& mesh) const;

  // Tessellate obeying the given PLC + points boundary.
  void tessellate(const std::vector<RealType>& points,
                  const std::vector<RealType>& PLCpoints,
                  const PLC<2, RealType>& geometry,
                  Tessellation<2, RealType>& mesh) const;

  // Tessellate obeying the given ReducedPLC boundary.
  void tessellate(const std::vector<RealType>& points,
                  const ReducedPLC<2, RealType>& geometry,
                  Tessellation<2, RealType>& mesh) const;

  // This Tessellator handles PLCs!
  bool handlesPLCs() const { return true; }

  // The name of the tessellator
  std::string name() const { return "BoostTessellator"; }

  //! Returns the accuracy to which this tessellator can distinguish coordinates.
  //! Should be returned appropriately for normalized coordinates, i.e., if all
  //! coordinates are in the range xi \in [0,1], what is the minimum allowed 
  //! delta in x.
  virtual RealType degeneracy() const { return mDegeneracy; }
  void degeneracy(const RealType val) const { mDegeneracy = val; }

private:
  //-------------------- Private interface ---------------------- //

  // ------------------------------------------------- //
  // Specialized tessellations based on the point set  //
  // ------------------------------------------------- //

  // Compute the nodes around a collection of generators
  void computeCellNodes(const std::vector<RealType>& points,
                        std::vector<std::vector<unsigned> >& cellNodes,
                        std::map<int, PointType>& id2node,
                        std::vector<unsigned>& infNodes) const;

  // Compute the nodes around a linear, 1d collection of generators
  void computeCellNodesCollinear(const std::vector<RealType>& points,
                                 std::vector<std::vector<unsigned> >& cellNodes,
                                 std::map<int, PointType>& id2node,
                                 std::vector<unsigned>& infNodes) const;

  void constructBoundedTopology(const std::vector<RealType>& points,
                                const ReducedPLC<2, RealType>& geometry,
                                const std::vector<ReducedPLC<2, CoordType> >& cellRings,
                                Tessellation<2, RealType>& mesh) const;

  // ----------------------------------------------------- //
  // Private tessellate calls used by internal algorithms  //
  // ----------------------------------------------------- //

  // Bounded tessellation with prescribed bounding box
  void tessellate(const std::vector<RealType>& points,
                  const ReducedPLC<2, CoordHash>& intGeometry,
                  const QuantizedCoordinates<2,RealType>& coords,
                  std::vector<ReducedPLC<2, CoordHash> >& intCells) const;

  // -------------------------- //
  // Private member variables   //
  // -------------------------- //

  // The quantized coordinates for this tessellator (inner and outer)
  static CoordHash coordMax;
  static RealType mDegeneracy; 
  mutable QuantizedCoordinates<2,RealType> mCoords;

  friend class BoostOrphanage<RealType>;
};



//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<typename RealType> 
typename BoostTessellator<RealType>::CoordHash
BoostTessellator<RealType>::coordMax = (1LL << 26);

template<typename RealType> 
RealType  
BoostTessellator<RealType>::mDegeneracy = 1.0/BoostTessellator<RealType>::coordMax;

} //end polytope namespace

#endif
#endif
