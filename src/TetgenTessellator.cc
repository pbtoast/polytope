//------------------------------------------------------------------------
// TetgenTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>

#include "polytope.hh" // Pulls in POLY_ASSERT and TetgenTessellator.hh.
#include "Point.hh"

// Pull in tetgen stuff.
#define TETLIBRARY
#include "tetgen.h"

namespace polytope {

using namespace std;

namespace {

//------------------------------------------------------------------------
// This function computes the determinant of the 3x3 matrix A with 
// the given 9 coefficients.
//------------------------------------------------------------------------
double 
det3(double A11, double A12, double A13,
     double A21, double A22, double A23,
     double A31, double A32, double A33) {
  return A11*(A22*A33-A32*A23) - A12*(A21*A33-A31*A23) + A13*(A21*A31-A31*A22);
}

//------------------------------------------------------------------------
// This function computes the determinant of the 4x4 matrix A with 
// the given 16 coefficients.
//------------------------------------------------------------------------
double 
det4(double A11, double A12, double A13, double A14,
     double A21, double A22, double A23, double A24,
     double A31, double A32, double A33, double A34,
     double A41, double A42, double A43, double A44) {
  double det31 = det3(A22, A23, A24,
                      A32, A33, A34,
                      A42, A43, A44);
  double det32 = det3(A21, A23, A24,
                      A31, A33, A34,
                      A41, A43, A44);
  double det33 = det3(A21, A22, A24,
                      A31, A32, A34,
                      A41, A42, A44);
  double det34 = det3(A21, A22, A23,
                      A31, A32, A33,
                      A41, A42, A43);
  return A11*det31 - A12*det32 + A13*det33 - A14*det34;
}

//------------------------------------------------------------------------
// This function computes the circumcenter of a tetrahedron with vertices
// A = (Ax, Ay, Az), B = (Bx, By, Bz), C = (Cx, Cy, Cz), and D = (Dx, Dy, Dz),
// and places the result in X.
//------------------------------------------------------------------------
void 
computeCircumcenter(double* A, double* B, double* C, double* D, double* X) {
  // This solution was taken from 
  // http://mathworld.wolfram.com/Circumsphere.html.
  double r12 = A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
  double r22 = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
  double r32 = C[0]*C[0] + C[1]*C[1] + C[2]*C[2];
  double r42 = D[0]*D[0] + D[1]*D[1] + D[2]*D[2];
  double a = det4(A[0], A[1], A[2], 1.0,
                  B[0], B[1], B[2], 1.0,
                  C[0], C[1], C[2], 1.0,
                  D[0], D[1], D[2], 1.0);
  double Dx = det4(r12, A[1], A[2], 1.0,
                   r22, B[1], B[2], 1.0,
                   r32, C[1], C[2], 1.0,
                   r42, D[1], D[2], 1.0);
  double Dy = det4(r12, A[0], A[2], 1.0,
                   r22, B[0], B[2], 1.0,
                   r32, C[0], C[2], 1.0,
                   r42, D[0], D[2], 1.0);
  double Dz = det4(r12, A[0], A[1], 1.0,
                   r22, B[0], B[1], 1.0,
                   r32, C[0], C[1], 1.0,
                   r42, D[0], D[1], 1.0);
  X[0] = 0.5*Dx/a;
  X[1] = 0.5*Dy/a;
  X[2] = 0.5*Dz/a;
}

//------------------------------------------------------------------------------
// An implementation of the map specialized to help constructing counters.
// This thing just overloads the index operator to start the count at zero
// for new key values.
//------------------------------------------------------------------------------
template<typename Key, 
         typename Comparator = std::less<Key> >
class CounterMap: public std::map<Key, unsigned> {
public:
  CounterMap(): std::map<Key, unsigned>() {}
  virtual ~CounterMap() {}
  unsigned& operator[](const Key& key) {
    typename std::map<Key, unsigned>::iterator itr = this->find(key);
    if (itr == this->end()) {
      std::map<Key, unsigned>::operator[](key) = 0U;
      itr = this->find(key);
    }
    POLY_ASSERT(itr != this->end());
    return itr->second;
  }
};

//------------------------------------------------------------------------------
// Increment the given position by 
//------------------------------------------------------------------------------
void
incrementPosition(Point3<double>& val,
                  const vector<double>& points,
                  const vector<int>& indices) {
  Point3<double> xi(0.0, 0.0, 0.0);
  for (unsigned k = 0; k != indices.size(); ++k) {
    const unsigned i = indices[k];
    POLY_ASSERT(i < points.size()/3);
    xi.x += points[3*i];
    xi.y += points[3*i + 1];
    xi.z += points[3*i + 2];
  }
  xi /= indices.size();
  val += xi;
}

} // end anonymous namespace

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
TetgenTessellator::
TetgenTessellator(const bool directComputation):
  Tessellator<3, double>(),
  mDirectComputation(directComputation) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
TetgenTessellator::
~TetgenTessellator() {
}

//------------------------------------------------------------------------------
// Unbounded tessellation.
//------------------------------------------------------------------------------
void
TetgenTessellator::
tessellate(const vector<double>& points,
           Tessellation<3, double>& mesh) const {
  PLC<3, double> geometry;
  this->tessellate(points, points, geometry, mesh);
}

//------------------------------------------------------------------------------
// Tessellate within a bounding box.
//------------------------------------------------------------------------------
void
TetgenTessellator::
tessellate(const vector<double>& points,
           double* low,
           double* high,
           Tessellation<3, double>& mesh) const {
  ReducedPLC<3, double> box = this->boundingBox(low, high);
  this->tessellate(points, box.points, box, mesh);
}

//------------------------------------------------------------------------------
// Tessellate within a PLC boundary.
//------------------------------------------------------------------------------
void
TetgenTessellator::
tessellate(const vector<double>& points,
           const vector<double>& PLCpoints,
           const PLC<3, double>& geometry,
           Tessellation<3, double>& mesh) const {
  if (mDirectComputation) {
    this->computeDirectVoronoi(points, PLCpoints, geometry, mesh);
  } else {
    this->computeDualOfTetrahedralization(points, PLCpoints, geometry, mesh);
  }
}

//------------------------------------------------------------------------------
// Internal method to build the tessellation directly using Tetgen's native 
// Voronoi capability.
//------------------------------------------------------------------------------
void
TetgenTessellator::
computeDirectVoronoi(const vector<double>& points,
                     const vector<double>& PLCpoints,
                     const PLC<3, double>& geometry,
                     Tessellation<3, double>& mesh) const {

  // Pre-conditions.
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 3 == 0);

  typedef int64_t CoordHash;
  typedef set<int> FaceHash;
  typedef Point3<CoordHash> IntPoint;
  typedef Point3<RealType> RealPoint;

  // Compute the normalized generators.
  RealType low[3], high[3];
  const unsigned numGenerators = points.size() / 3;
  vector<double> generators = this->computeNormalizedPoints(points, PLCpoints, low, high);

  // Build the input to tetgen.
  tetgenio in;
  in.firstnumber = 0;
  in.mesh_dim = 3;
  in.pointlist = &generators.front();
  in.pointattributelist = 0;
  in.pointmtrlist = 0;
  in.pointmarkerlist = 0;
  in.numberofpoints = generators.size() / 3;
  in.numberofpointattributes = 0;
  in.numberofpointmtrs = 0;

  // Copy the PLC boundary info to the tetgen input.
  vector<tetgenio::polygon> plcFacetPolygons(geometry.facets.size());
  for (unsigned ifacet = 0; ifacet != geometry.facets.size(); ++ifacet) {
    POLY_ASSERT(geometry.facets[ifacet].size() >= 3);
    plcFacetPolygons[ifacet].numberofvertices = geometry.facets[ifacet].size();
    plcFacetPolygons[ifacet].vertexlist = const_cast<int*>(&geometry.facets[ifacet].front());
  }
  unsigned nfacets = plcFacetPolygons.size();
  vector<RealPoint> holeCentroids(geometry.holes.size(), RealPoint(0.0, 0.0, 0.0));
  for (unsigned ihole = 0; ihole != geometry.holes.size(); ++ihole) {
    plcFacetPolygons.resize(nfacets + geometry.holes[ihole].size());
    for (unsigned ifacet = 0; ifacet != geometry.holes[ihole].size(); ++ifacet) {
      POLY_ASSERT(geometry.holes[ihole][ifacet].size() >= 3);
      plcFacetPolygons[nfacets + ifacet].numberofvertices = geometry.holes[ihole][ifacet].size();
      plcFacetPolygons[nfacets + ifacet].vertexlist = const_cast<int*>(&geometry.holes[ihole][ifacet].front());
      incrementPosition(holeCentroids[ihole], PLCpoints, geometry.holes[ihole][ifacet]);
    }
    holeCentroids[ihole] /= RealType(geometry.holes[ihole].size());
    nfacets = plcFacetPolygons.size();
  }
  POLY_ASSERT(nfacets == plcFacetPolygons.size());
  vector<tetgenio::facet> plcFacets(nfacets);
  for (unsigned ifacet = 0; ifacet != nfacets; ++ifacet) {
    in.init(&plcFacets[ifacet]);
    plcFacets[ifacet].polygonlist = &plcFacetPolygons[ifacet];
    plcFacets[ifacet].numberofpolygons = 1;
  }
  in.numberoffacets = nfacets;
  in.facetlist = &plcFacets.front();

  // Do the tetrahedralization.
  tetgenio out;
  tetrahedralize((char*)"vV", &in, &out);

  // Make sure we got something.
  if (out.numberoftetrahedra == 0)
    error("TetgenTessellator: Delauney tetrahedralization produced 0 tetrahedra!");
  if (out.numberofpoints != generators.size()/3) {
    char err[1024];
    snprintf(err, 1024, "TetgenTessellator: Delauney tetrahedralization produced %d tetrahedra\n(%d generating points given)", 
             out.numberofpoints, (int)numGenerators);
    error(err);
  }

  // // Find the circumcenters of each tetrahedron.
  // RealType  clow[3] = { numeric_limits<RealType>::max(),  
  //                       numeric_limits<RealType>::max(),
  //                       numeric_limits<RealType>::max()};
  // RealType chigh[3] = {-numeric_limits<RealType>::max(), 
  //                      -numeric_limits<RealType>::max(),
  //                      -numeric_limits<RealType>::max()};
  // CounterMap<FaceHash> faceCounter;
  // vector<RealPoint> circumcenters(out.numberoftetrahedra);
  // map<int, set<int> > gen2tet;
  // map<int, set<int> > neighbors;
  // unsigned p1, p2, p3, p4;
  // for (unsigned i = 0; i != out.numberoftetrahedra; ++i) {
  // }
}

//------------------------------------------------------------------------------
// Internal method to build the tessellation by first doing the 
// tetrahedralization, and then computing the dual.
//------------------------------------------------------------------------------
void
TetgenTessellator::
computeDualOfTetrahedralization(const vector<double>& points,
                                const vector<double>& PLCpoints,
                                const PLC<3, double>& geometry,
                                Tessellation<3, double>& mesh) const {

  // Pre-conditions.
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 3 == 0);

  typedef int64_t CoordHash;
  typedef set<int> FaceHash;
  typedef Point3<CoordHash> IntPoint;
  typedef Point3<RealType> RealPoint;

  // // Normalize coordinates in the input box.
  // vector<double> generators;
  // const unsigned numGenerators = points.size() / 3;
  // generators.reserve(3*(numGenerators + 8));
  // double box[3] = {high[0] - low[0],
  //                  high[1] - low[1],
  //                  high[2] - low[2]};
  // POLY_ASSERT(box[0] > 0.0 and box[1] > 0.0 and box[2] > 0.0);
  // for (unsigned i = 0; i != points.size(); ++i) {
  //   const unsigned j = i % 3;
  //   generators.push_back((points[i] - low[j])/box[j]);
  //   POLY_ASSERT(generators.back() >= 0.0 and generators.back() <= 1.0);
  // }
  // POLY_ASSERT(generators.size() == points.size());

  // // Add the boundary points to the generator list.
  // generators.push_back(0.0); generators.push_back(0.0); generators.push_back(0.0);
  // generators.push_back(1.0); generators.push_back(0.0); generators.push_back(0.0);
  // generators.push_back(1.0); generators.push_back(1.0); generators.push_back(0.0);
  // generators.push_back(0.0); generators.push_back(1.0); generators.push_back(0.0);
  // generators.push_back(0.0); generators.push_back(0.0); generators.push_back(1.0);
  // generators.push_back(1.0); generators.push_back(0.0); generators.push_back(1.0);
  // generators.push_back(1.0); generators.push_back(1.0); generators.push_back(1.0);
  // generators.push_back(0.0); generators.push_back(1.0); generators.push_back(1.0);
  // POLY_ASSERT(generators.size() == 3*(numGenerators + 8));

  // // Build the input to tetgen.
  // tetgenio in;
  // in.firstnumber = 0;
  // in.mesh_dim = 3;
  // in.pointlist = &generators.front();
  // in.pointattributelist = 0;
  // in.pointmtrlist = 0;
  // in.pointmarkerlist = 0;
  // in.numberofpoints = generators.size() / 3;
  // in.numberofpointattributes = 0;
  // in.numberofpointmtrs = 0;

  // // Do the tetrahedralization.
  // tetgenio out;
  // tetrahedralize((char*)"qV", &in, &out);

  // // Make sure we got something.
  // if (out.numberoftetrahedra == 0)
  //   error("TetgenTessellator: Delauney tetrahedralization produced 0 tetrahedra!");
  // if (out.numberofpoints != generators.size()/3) {
  //   char err[1024];
  //   snprintf(err, 1024, "TetgenTessellator: Delauney tetrahedralization produced %d tetrahedra\n(%d generating points given)", 
  //            out.numberofpoints, (int)numGenerators);
  //   error(err);
  // }

  // // Find the circumcenters of each tetrahedron.
  // RealType  clow[3] = { numeric_limits<RealType>::max(),  
  //                       numeric_limits<RealType>::max(),
  //                       numeric_limits<RealType>::max()};
  // RealType chigh[3] = {-numeric_limits<RealType>::max(), 
  //                      -numeric_limits<RealType>::max(),
  //                      -numeric_limits<RealType>::max()};
  // CounterMap<FaceHash> faceCounter;
  // vector<RealPoint> circumcenters(out.numberoftetrahedra);
  // map<int, set<int> > gen2tet;
  // map<int, set<int> > neighbors;
  // unsigned p1, p2, p3, p4;
  // for (unsigned i = 0; i != out.numberoftetrahedra; ++i) {
  // }
}

//------------------------------------------------------------------------------
// Access the internal attribute to compute directly or not.
//------------------------------------------------------------------------------
bool
TetgenTessellator::
directComputation() const {
  return mDirectComputation;
}

void
TetgenTessellator::
directComputation(const bool x) {
  mDirectComputation = x;
}

}
