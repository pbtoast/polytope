// Unit tests for Computational Solid Geometry (CSG) operations on PLC's.

#include <iostream>
#include <vector>
#include <cassert>
#include <ctime>

#include "polytope.hh"
#include "PLC_CSG.hh"
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
// Simplify the facets of a PLC -- combine coplanar facets and remove collinear 
// points around the facets.
//------------------------------------------------------------------------------
template<typename RealType>
PLC<3, RealType>
simplifyPLCfacets(const PLC<3, RealType>& plc,
                  const std::vector<RealType>& points) {
  // PLC<3, RealType> result;
  
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
    const ReducedPLC<3, double> box2 = CSG::CSG_internal::ReducedPLCfromPolygons(polys);
    // escapePod("box1", box1);
    // escapePod("box1_from_polys", box2);
    POLY_CHECK2(polys.size() == 12, "Num polys = " << polys.size() << ", expected 12.");
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
    escapePod("box1", box1);
    escapePod("box2", box2);
    cerr << escapePod("box_union_test", box_union) << endl;
    const ReducedPLC<3, double> box_intersect = CSG::csg_intersect(box1, box2);
    cerr << escapePod("box_intersect_test", box_intersect) << endl;
  }

  cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
  return 0;
}
