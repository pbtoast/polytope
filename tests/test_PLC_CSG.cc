// Unit tests for Computational Solid Geometry (CSG) operations on PLC's.

#include <iostream>
#include <vector>
#include <cassert>
#include <ctime>

#include "polytope.hh"
#include "PLC_CSG_3d.hh"
#include "simplifyPLCfacets.hh"
#include "polytope_internal.hh"
#include "polytope_geometric_utilities.hh"
#include "polytope_test_utilities.hh"
#include "polytope_write_OOGL.hh"
#include "polytope_plc_canned_geometries.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

namespace {

//------------------------------------------------------------------------------
// Rotate a set of points.
//------------------------------------------------------------------------------
template<typename RealType>
void
rotatePoints(std::vector<RealType>& points,
             const RealType theta,
             const RealType phi,
             const RealType psi) {
  POLY_ASSERT(points.size() % 3 == 0);
  const RealType R[3][3] = {{ cos(theta)*cos(psi), cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi), sin(phi)*sin(psi) - cos(phi)*sin(theta)*cos(psi)},
                            {-cos(theta)*sin(psi), cos(phi)*cos(psi) - sin(phi)*sin(theta)*sin(psi), sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi)},
                            { sin(theta),         -sin(phi)*cos(theta),                              cos(phi)*cos(theta)}};
  const unsigned n = points.size()/3;
  for (unsigned i = 0; i != n; ++i) {
    const RealType x = points[3*i], y = points[3*i+1], z = points[3*i+2];
    points[3*i  ] = R[0][0]*x + R[0][1]*y + R[0][2]*z;
    points[3*i+1] = R[1][0]*x + R[1][1]*y + R[1][2]*z;
    points[3*i+2] = R[2][0]*x + R[2][1]*y + R[2][2]*z;
  }
}

//------------------------------------------------------------------------------
// Emergency dump.
//------------------------------------------------------------------------------
template<typename RealType>
std::string
escapePod(const std::string nameEnd,
          const ReducedPLC<3, RealType>& plc) {
    std::stringstream os;
    os << "test_PLC_CSG_" << nameEnd;
    writePLCtoOFF(plc, plc.points, os.str());
    return " : attempted to write to file " + os.str();
}

} // anonymous namespace

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv) {
  typedef Point3<double> PointType;
  typedef Point3<int64_t> IntPointType;

#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // //----------------------------------------------------------------------
  // // Convert from PLC->CSG_internal_3d::Polygons->PLC
  // //----------------------------------------------------------------------
  // {
  //   double low[3] = {0.0, 0.0, 0.0}, high[3] = {1.0, 1.0, 1.0};
  //   const ReducedPLC<3, double> box1 = plc_box<3, double>(low, high);
  //   const std::vector<CSG::CSG_internal_3d::Polygon<double> > polys = CSG::CSG_internal_3d::ReducedPLCtoPolygons(box1);
  //   POLY_CHECK2(polys.size() == 12, "Num polys = " << polys.size() << ", expected 12.");
  //   const ReducedPLC<3, double> box2 = CSG::CSG_internal_3d::ReducedPLCfromPolygons(polys);
  //   const ReducedPLC<3, double> box3 = polytope::simplifyPLCfacets(box2, box2.points, low, high, 1.0e-10);
  //   POLY_CHECK(box3.facets.size() == 6);
  //   POLY_CHECK(box3.points.size() == 3*8);
  // }

  //----------------------------------------------------------------------
  // Operate on two boxes offset in x.
  //----------------------------------------------------------------------
  {
    double low1[3] = {0.0, 0.0, 0.0}, high1[3] = {1.0, 1.0, 1.0},
           low2[3] = {0.5, 0.0, 0.0}, high2[3] = {1.5, 1.0, 1.0};
    const ReducedPLC<3, double> box1 = plc_box<3, double>(low1, high1),
                                box2 = plc_box<3, double>(low2, high2);
    const ReducedPLC<3, double> box_union = CSG::csg_union(box1, box2);
    POLY_CHECK(box_union.points.size() % 3 == 0);
    // escapePod("box1", box1);
    // escapePod("box2", box2);
    // cerr << escapePod("box_union_test", box_union) << endl;
    const ReducedPLC<3, double> box_union_simplify = polytope::simplifyPLCfacets(box_union, box_union.points, low1, high2, 1.0e-10);
    cerr << escapePod("box_union_test_simplify", box_union_simplify) << endl;
    const ReducedPLC<3, double> box_intersect = CSG::csg_intersect(box1, box2);
    // cerr << escapePod("box_intersect_test", box_intersect) << endl;
    const ReducedPLC<3, double> box_intersect_simplify = polytope::simplifyPLCfacets(box_intersect, box_intersect.points, low1, high2, 1.0e-10);
    // cerr << escapePod("box_intersect_test_simplify", box_intersect_simplify) << endl;
    POLY_CHECK(box_intersect_simplify.facets.size() == 6);
    POLY_CHECK(box_intersect_simplify.points.size() == 3*8);
  }

  // //----------------------------------------------------------------------
  // // Generate a complicated sphere with holes cut out of it.
  // //----------------------------------------------------------------------
  // {
  //   const PointType origin(0.0, 0.0, 0.0);
  //   const unsigned nphi = 20;
  //   clock_t t0 = clock();
  //   double low[3] = {-1.0, -1.0, -1.0}, high[3] = {1.0, 1.0, 1.0};
  //   const ReducedPLC<3, double> a = plc_box<3, double>(low, high),
  //                               b = plc_sphere<double>(origin, 1.35, nphi),
  //                               c = plc_cylinder(origin, 0.7, 2.0, nphi);
  //   ReducedPLC<3, double> d(c), e(c);
  //   rotatePoints(d.points, 0.0, M_PI/2.0, 0.0);
  //   rotatePoints(e.points, 0.0, M_PI/2.0, M_PI/2.0);
  //   clock_t t1 = clock();
  //   cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds to generate input PLCs." << endl;
  //   t0 = clock();
  //   const ReducedPLC<3, double> holey_sphere = CSG::csg_subtract(CSG::csg_intersect(a, b), 
  //                                                                CSG::csg_union(CSG::csg_union(c, d), e));
  //   t1 = clock();
  //   cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds to perform CSG operations." << endl;
  //   escapePod("holey_sphere", holey_sphere);
  // }

  // ********************************************************************************
  // Repeat above tests for int types.
  //----------------------------------------------------------------------
  // Convert from PLC->CSG_internal_3d::Polygons->PLC
  //----------------------------------------------------------------------
  {
    int64_t low[3] = {0, 0, 0}, high[3] = {1 << 10, 1 << 10, 1 << 10};
    const ReducedPLC<3, int64_t> box1 = plc_box<3, int64_t>(low, high);
    const std::vector<CSG::CSG_internal_3d::Polygon<int64_t> > polys = CSG::CSG_internal_3d::ReducedPLCtoPolygons(box1);
    POLY_CHECK2(polys.size() == 12, "Num polys = " << polys.size() << ", expected 12.");
    const ReducedPLC<3, int64_t> box2 = CSG::CSG_internal_3d::ReducedPLCfromPolygons(polys);
    const ReducedPLC<3, int64_t> box3 = polytope::simplifyPLCfacets(box2, box2.points, low, high, int64_t(1));
    POLY_CHECK(box3.facets.size() == 6);
    POLY_CHECK(box3.points.size() == 3*8);
  }

  //----------------------------------------------------------------------
  // Operate on two boxes offset in x.
  //----------------------------------------------------------------------
  {
    int64_t low1[3] = {0, 0, 0},          high1[3] = {1 << 10, 1 << 10, 1 << 10},
            low2[3] = {high1[0]/2, 0, 0}, high2[3] = {low2[0] + high1[0], high1[1], high1[2]};
    const ReducedPLC<3, int64_t> box1 = plc_box<3, int64_t>(low1, high1),
                                 box2 = plc_box<3, int64_t>(low2, high2);
    const ReducedPLC<3, int64_t> box_union = CSG::csg_union(box1, box2);
    POLY_CHECK(box_union.points.size() % 3 == 0);
    escapePod("box1_int", box1);
    escapePod("box2_int", box2);
    cerr << escapePod("box_union_test_int", box_union) << endl;
    const ReducedPLC<3, int64_t> box_union_simplify = polytope::simplifyPLCfacets(box_union, box_union.points, low1, high2, int64_t(1));
    cerr << escapePod("box_union_test_simplify_int", box_union_simplify) << endl;
    const ReducedPLC<3, int64_t> box_intersect = CSG::csg_intersect(box1, box2);
    cerr << escapePod("box_intersect_test_int", box_intersect) << endl;
    const ReducedPLC<3, int64_t> box_intersect_simplify = polytope::simplifyPLCfacets(box_intersect, box_intersect.points, low1, high2, int64_t(1));
    cerr << escapePod("box_intersect_test_simplify_int", box_intersect_simplify) << endl;
    POLY_CHECK(box_intersect_simplify.facets.size() == 6);
    POLY_CHECK(box_intersect_simplify.points.size() == 3*8);
  }

  cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
  return 0;
}
