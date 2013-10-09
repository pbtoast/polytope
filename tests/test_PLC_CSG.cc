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

#if HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

namespace {

//------------------------------------------------------------------------------
// Create a single zone pseudo-tessellation from a PLC.
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
// Create a single zone pseudo-tessellation from a PLC.
//------------------------------------------------------------------------------
template<typename RealType>
void
tessellationFromPLC(const PLC<3, RealType>& plc,
                    const std::vector<RealType>& points,
                    Tessellation<3, RealType>& mesh) {
  // Nodes.
  mesh.nodes = points;

  // Faces.
  const unsigned nfaces = plc.facets.size();
  mesh.faces = vector<vector<unsigned> >(nfaces);
  for (unsigned i = 0; i != nfaces; ++i) {
    const unsigned n = plc.facets[i].size();
    for (unsigned j = 0; j != n; ++j) mesh.faces[i].push_back(plc.facets[i][j]);
  }

  // Face cells.
  mesh.faceCells.resize(nfaces, vector<int>(size_t(1), int(0)));

  // The one cell.
  mesh.cells.resize(1);
  for (int i = 0; i != nfaces; ++i) mesh.cells[0].push_back(i);
}

//------------------------------------------------------------------------------
// Emergency dump.
//------------------------------------------------------------------------------
std::string
escapePod(const std::string nameEnd,
          const ReducedPLC<3, double>& plc) {
    std::stringstream os;
    os << "test_PLC_CSG_" << nameEnd;
    Tessellation<3, double> mesh;
    tessellationFromPLC(plc, plc.points, mesh);
    outputMesh(mesh, os.str());
    return " : attempted to write to file " + os.str();
}
}

//------------------------------------------------------------------------------
// Return a 3D PLC box.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<3, RealType>
boxPLC(const RealType x1, const RealType x2,
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
    const ReducedPLC<3, double> box1 = boxPLC<double>(0.0, 1.0,
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
    const ReducedPLC<3, double> box1 = boxPLC<double>(0.0, 1.0,
                                                      0.0, 1.0,
                                                      0.0, 1.0),
                                box2 = boxPLC<double>(0.5, 1.5,
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
    // cerr << "Simplified PLC: " << box_intersect_simplify << endl
    //      << "  points: " << endl;
    // for (unsigned i = 0; i != box_intersect_simplify.points.size()/3; ++i) cerr << "    " << i << " (" << box_intersect_simplify.points[3*i] << " " << box_intersect_simplify.points[3*i+1] << " " << box_intersect_simplify.points[3*i+2] << ")" << endl;
    POLY_CHECK(box_intersect_simplify.facets.size() == 6);
    POLY_CHECK(box_intersect_simplify.points.size() == 3*8);
  }

  //----------------------------------------------------------------------
  // Generate a complicated sphere with holes cut out of it.
  //----------------------------------------------------------------------
  {
    const PointType origin(0.0, 0.0, 0.0);
    const double router = 2.0;
    const unsigned nphi = 2;
    const ReducedPLC<3, double> outer_sphere = plc_sphere<double>(origin, router, nphi);
    cerr << "sphere : " << outer_sphere << endl;
    cerr << "points : " << endl;
    for (unsigned i = 0; i != outer_sphere.points.size()/3; ++i) cerr << "    " << i << " (" << outer_sphere.points[3*i] << " " << outer_sphere.points[3*i+1] << " " << outer_sphere.points[3*i+2] << ")" << endl;
    escapePod("outer_sphere", outer_sphere);
  }

  cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
  return 0;
}
