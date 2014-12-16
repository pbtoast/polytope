#include <iostream>
#include <vector>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>

#include "polytope_test_utilities.hh"

#define BOOST_GEOMETRY_NO_ROBUSTNESS true

using namespace std;

//------------------------------------------------------------------------
// main
//------------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif


  typedef int32_t CoordType;
  typedef boost::geometry::model::d2::point_xy<CoordType> PointType;
  typedef boost::geometry::model::polygon<PointType> PolygonType;
  typedef boost::geometry::model::ring<PointType> RingType;

  { // Test 1
    PolygonType a, b;
    boost::geometry::read_wkt("POLYGON((42817136 -3774506, 43029074 -3929862, 31446819 18947953, 30772384 19615678, 30101303 19612322, 30114725 16928001, 33520458 6878575, 35332375 2413654, 35725796 2024148))", a);
    boost::geometry::read_wkt("POLYGON((-33386239 -33721784, 33721785 -33386239, 33386240 33721785, -33721784 33386240))", b);
    boost::geometry::correct(a);
    boost::geometry::correct(b);
    POLY_CHECK(not boost::geometry::intersects(a));
    POLY_CHECK(not boost::geometry::intersects(b));
    std::vector<PolygonType> result;
    boost::geometry::intersection(a, b, result);
    POLY_CHECK2(result.size() == 1, result.size());
    boost::geometry::correct(result[0]);
    for (RingType::const_iterator itr = result[0].outer().begin();
         itr != result[0].outer().end();
         ++itr)   cout << itr->get<0>() << " " << itr->get<1>() << endl;
    POLY_CHECK(not boost::geometry::intersects(result[0]));
  }

#if HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
//------------------------------------------------------------------------
