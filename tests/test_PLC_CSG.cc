// Unit tests for Computational Solid Geometry (CSG) operations on PLC's.

#include <iostream>
#include <vector>
#include <cassert>
#include <ctime>

#include "polytope.hh"
#include "PLC_CSG.hh"
#include "simplifyPLCfacets.hh"
#include "polytope_internal.hh"
#include "polytope_geometric_utilities.hh"
#include "polytope_test_utilities.hh"
#include "polytope_write_OOGL.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

namespace {

//------------------------------------------------------------------------------
// Create a PLC representation of a sphere.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<3, RealType>
plc_sphere(const Point3<RealType>& center,
           const RealType radius,
           const unsigned nphi) {
  POLY_ASSERT(nphi > 0);

  // Prepare the result.
  ReducedPLC<3, RealType> result;

  // Put a bunch of points on the sphere.
  const double dphi = M_PI/nphi;
  for (unsigned iphi = 0; iphi != nphi; ++iphi) {
    const double phi = (iphi + 0.5)*dphi;
    const double dl = radius*dphi;
    const double rxy = radius*sin(phi);
    const double circ = 2.0*M_PI*rxy;
    const unsigned ntheta = std::max(1U, unsigned(circ/dl + 0.5));
    const double dtheta = 2.0*M_PI/ntheta;
    for (unsigned itheta = 0; itheta != ntheta; ++itheta) {
      const double theta = (itheta + 0.5)*dtheta;
      result.points.push_back(center.x + radius*cos(theta)*sin(phi));
      result.points.push_back(center.y + radius*sin(theta)*sin(phi));
      result.points.push_back(center.z + radius*cos(phi));
    }
  }

  // Build the convex hull of our points.
  const RealType low[3] = {center.x - 1.1*radius,
                           center.y - 1.1*radius,
                           center.z - 1.1*radius};
  const RealType dx = radius/(1 << 21);
  const PLC<3, RealType> hull = polytope::convexHull_3d(result.points, low, dx);

  // Put that topology in our result and we're done.
  result.facets = hull.facets;
  return result;
}

//------------------------------------------------------------------------------
// Create a PLC representation of a cylinder.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<3, RealType>
plc_cylinder(const Point3<RealType>& center,
             const RealType radius,
             const RealType length,
             const unsigned ncirc) {
  POLY_ASSERT(ncirc >= 3);

  // Prepare the result.
  ReducedPLC<3, RealType> result;

  // Put points on the two end circles.
  const double dtheta = 2.0*M_PI/ncirc;
  for (unsigned itheta = 0; itheta != ncirc; ++itheta) {
    const double theta = (itheta + 0.5)*dtheta;
    result.points.push_back(center.x + radius*cos(theta));
    result.points.push_back(center.y + radius*sin(theta));
    result.points.push_back(center.z + 0.5*length);
    result.points.push_back(center.x + radius*cos(theta));
    result.points.push_back(center.y + radius*sin(theta));
    result.points.push_back(center.z - 0.5*length);
  }

  // Build the convex hull of our points.
  const RealType low[3] = {center.x - 1.1*radius,
                           center.y - 1.1*radius,
                           center.z - 0.6*length};
  const RealType dx = length/(1 << 21);
  const PLC<3, RealType> hull = polytope::convexHull_3d(result.points, low, dx);

  // Put that topology in our result and we're done.
  result.facets = hull.facets;
  return result;
}

//------------------------------------------------------------------------------
// Return a 3D PLC box.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<3, RealType>
plc_box(const RealType x1, const RealType x2,
        const RealType y1, const RealType y2,
        const RealType z1, const RealType z2) {

  // Create the piecewise linear complex representing the box. Note that 
  // the box consists of facets that are defined by their connections to 
  // generating points.
  // Should look like the following:
  //
  //        6--------7            y
  //       /        /|            |
  //      /        / |            |
  //     2--------3  |             ------x
  //     |  .     |  |           /
  //     |  4.....|..5          z
  //     | .      | / 
  //     |.       |/
  //     0--------1             
  //
  // Create the vertices for our bounding surface.
  ReducedPLC<3, RealType> box;
  box.points.resize(3*8);
  box.points[3*0+0] = x1; box.points[3*0+1] = y1; box.points[3*0+2] = z2;
  box.points[3*1+0] = x2; box.points[3*1+1] = y1; box.points[3*1+2] = z2;
  box.points[3*2+0] = x1; box.points[3*2+1] = y2; box.points[3*2+2] = z2;
  box.points[3*3+0] = x2; box.points[3*3+1] = y2; box.points[3*3+2] = z2;
  box.points[3*4+0] = x1; box.points[3*4+1] = y1; box.points[3*4+2] = z1;
  box.points[3*5+0] = x2; box.points[3*5+1] = y1; box.points[3*5+2] = z1;
  box.points[3*6+0] = x1; box.points[3*6+1] = y2; box.points[3*6+2] = z1;
  box.points[3*7+0] = x2; box.points[3*7+1] = y2; box.points[3*7+2] = z1;

  // 6 facets
  box.facets.resize(6);

  // facet 0 -- bottom face.
  box.facets[0].resize(4);
  box.facets[0][0] = 0;
  box.facets[0][1] = 4;
  box.facets[0][2] = 5;
  box.facets[0][3] = 1;

  // facet 1 -- top face.
  box.facets[1].resize(4);
  box.facets[1][0] = 2;
  box.facets[1][1] = 3;
  box.facets[1][2] = 7;
  box.facets[1][3] = 6;

  // facet 2 -- left face.
  box.facets[2].resize(4);
  box.facets[2][0] = 0;
  box.facets[2][1] = 2;
  box.facets[2][2] = 6;
  box.facets[2][3] = 4;

  // facet 3 -- right face.
  box.facets[3].resize(4);
  box.facets[3][0] = 1;
  box.facets[3][1] = 5;
  box.facets[3][2] = 7;
  box.facets[3][3] = 3;

  // facet 4 -- front face.
  box.facets[4].resize(4);
  box.facets[4][0] = 0;
  box.facets[4][1] = 1;
  box.facets[4][2] = 3;
  box.facets[4][3] = 2;

  // facet 5 -- back face.
  box.facets[5].resize(4);
  box.facets[5][0] = 5;
  box.facets[5][1] = 4;
  box.facets[5][2] = 6;
  box.facets[5][3] = 7;

  // That's it.
  return box;
}

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
std::string
escapePod(const std::string nameEnd,
          const ReducedPLC<3, double>& plc) {
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

#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  //----------------------------------------------------------------------
  // Convert from PLC->CSG_internal::Polygons->PLC
  //----------------------------------------------------------------------
  {
    const ReducedPLC<3, double> box1 = plc_box<double>(0.0, 1.0,
                                                       0.0, 1.0,
                                                       0.0, 1.0);
    const std::vector<CSG::CSG_internal::Polygon<double> > polys = CSG::CSG_internal::ReducedPLCtoPolygons(box1);
    POLY_CHECK2(polys.size() == 12, "Num polys = " << polys.size() << ", expected 12.");
    const ReducedPLC<3, double> box2 = CSG::CSG_internal::ReducedPLCfromPolygons(polys);
    const ReducedPLC<3, double> box3 = polytope::simplifyPLCfacets(box2, box2.points, 1.0e-10);
    POLY_CHECK(box3.facets.size() == 6);
    POLY_CHECK(box3.points.size() == 3*8);
  }

  //----------------------------------------------------------------------
  // Operate on two boxes offset in x.
  //----------------------------------------------------------------------
  {
    const ReducedPLC<3, double> box1 = plc_box<double>(0.0, 1.0,
                                                       0.0, 1.0,
                                                       0.0, 1.0),
                                box2 = plc_box<double>(0.5, 1.5,
                                                       0.0, 1.0,
                                                       0.0, 1.0);
    const ReducedPLC<3, double> box_union = CSG::csg_union(box1, box2);
    POLY_CHECK(box_union.points.size() % 3 == 0);
    // escapePod("box1", box1);
    // escapePod("box2", box2);
    // cerr << escapePod("box_union_test", box_union) << endl;
    const ReducedPLC<3, double> box_intersect = CSG::csg_intersect(box1, box2);
    // cerr << escapePod("box_intersect_test", box_intersect) << endl;
    const ReducedPLC<3, double> box_intersect_simplify = polytope::simplifyPLCfacets(box_intersect, box_intersect.points, 1.0e-10);
    // cerr << escapePod("box_intersect_test_simplify", box_intersect_simplify) << endl;
    POLY_CHECK(box_intersect_simplify.facets.size() == 6);
    POLY_CHECK(box_intersect_simplify.points.size() == 3*8);
  }

  //----------------------------------------------------------------------
  // Generate a complicated sphere with holes cut out of it.
  //----------------------------------------------------------------------
  {
    const PointType origin(0.0, 0.0, 0.0);
    const unsigned nphi = 20;
    clock_t t0 = clock();
    const ReducedPLC<3, double> a = plc_box<double>(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0),
                                b = plc_sphere<double>(origin, 1.35, nphi),
                                c = plc_cylinder(origin, 0.7, 2.0, nphi);
    ReducedPLC<3, double> d(c), e(c);
    rotatePoints(d.points, 0.0, M_PI/2.0, 0.0);
    rotatePoints(e.points, 0.0, M_PI/2.0, M_PI/2.0);
    clock_t t1 = clock();
    cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds to generate input PLCs." << endl;
    t0 = clock();
    const ReducedPLC<3, double> holey_sphere = CSG::csg_subtract(CSG::csg_intersect(a, b), 
                                                                 CSG::csg_union(CSG::csg_union(c, d), e));
    t1 = clock();
    cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds to perform CSG operations." << endl;
    escapePod("holey_sphere", holey_sphere);
  }

  cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
  return 0;
}
