//---------------------------------Spheral++----------------------------------//
// VoroPP_3d
//----------------------------------------------------------------------------//
#include <stdint.h>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>
#include <set>

#include "polytope.hh" // Pulls in ASSERT and VoroPP_3d.hh.
#include "container.hh"

namespace polytope {

using namespace std;
using namespace voro;
using std::min;
using std::max;
using std::abs;

namespace { // We hide internal functions in an anonymous namespace.

//------------------------------------------------------------------------------
// A integer version of the simple 3D point.
//------------------------------------------------------------------------------
template<typename Uint>
struct Point3 {
  Uint x, y, z;
  Point3(): x(0), y(0), z(0) {}
  Point3(const Uint& xi, const Uint& yi, const Uint& zi): x(xi), y(yi), z(zi) {}
  Point3& operator=(const Point3& rhs) { x = rhs.x; y = rhs.y; z = rhs.z; return *this; }
  bool operator==(const Point3& rhs) const { return (x == rhs.x and y == rhs.y and z == rhs.z); }
  bool operator<(const Point3& rhs) const {
    return (x < rhs.x                               ? true :
            x == rhs.x and y < rhs.y                ? true :
            x == rhs.x and y == rhs.y and z < rhs.z ? true :
            false);
  }
  template<typename Real>
  Point3(const Real& xi, const Real& yi, const Real& zi, const Real& dx): 
    x(static_cast<Uint>(xi/dx + 0.5)),
    y(static_cast<Uint>(yi/dx + 0.5)),
    z(static_cast<Uint>(zi/dx + 0.5)) {}
  template<typename Real> Real realx(const Real& xmin, const Real& dx) { return static_cast<Real>(x*dx) + xmin; }
  template<typename Real> Real realy(const Real& ymin, const Real& dy) { return static_cast<Real>(y*dy) + ymin; }
  template<typename Real> Real realz(const Real& zmin, const Real& dz) { return static_cast<Real>(z*dz) + zmin; }
};

// It's nice being able to print these things.
template<typename Uint>
std::ostream&
operator<<(std::ostream& os, const Point3<Uint>& p) {
  os << "(" << p.x << " " << p.y << " " << p.z <<  ")";
  return os;
}

//------------------------------------------------------------------------------
// Take the given set of vertex positions, and either find them in the known
// mesh nodes or add them to the known set.
//------------------------------------------------------------------------------
template<typename Real, typename Uint>
map<unsigned, unsigned>
updateMeshVertices(vector<Point3<Uint> >& vertices,
                   map<Point3<Uint>, unsigned>& vertexHash2ID,
                   Tessellation<3, Real>& mesh,
                   const Real& xmin,
                   const Real& ymin,
                   const Real& zmin,
                   const Real& fconv) {
  const unsigned n = vertices.size();
  bool newVertex;
  unsigned i, j;
  Uint ix, iy, iz, ix0, iy0, iz0, ix1, iy1, iz1;
  Point3<Uint> ipt, ipt1;
  map<unsigned, unsigned> result;
  for (i = 0; i !=n; ++i) {
    ipt = vertices[i];
    ix0 = ipt.x > 0 ? ipt.x - 1 : ipt.x;
    iy0 = ipt.y > 0 ? ipt.y - 1 : ipt.y;
    iz0 = ipt.z > 0 ? ipt.z - 1 : ipt.z;
    ix1 = ipt.x + 1;
    iy1 = ipt.y + 1;
    iz1 = ipt.z + 1;
    newVertex = true;
    iz = iz0;
    while (newVertex and iz != iz1) {
      iy = iy0;
      while (newVertex and iy != iy1) {
        ix = ix0;
        while (newVertex and ix != ix1) {
          ipt1 = Point3<Uint>(ix, iy, iz);
          newVertex = (vertexHash2ID.find(ipt1) == vertexHash2ID.end());
          ++ix;
        }
        ++iy;
      }
      ++iz;
    }
    if (newVertex) {
      j = vertexHash2ID.size();
      vertexHash2ID[ipt] = j;
      mesh.nodes.push_back(vertices[i].realx(xmin, fconv));
      mesh.nodes.push_back(vertices[i].realy(ymin, fconv));
      mesh.nodes.push_back(vertices[i].realz(zmin, fconv));
      ASSERT(mesh.nodes.size()/3 == j + 1);
      result[i] = j;
    } else {
      result[i] = vertexHash2ID[ipt1];
    }
  }
  ASSERT(result.size() == vertices.size());
  return result;
}

//------------------------------------------------------------------------------
// Helper method to update our face info.
//------------------------------------------------------------------------------
template<typename Real>
inline
void
insertFaceInfo(const set<unsigned>& fhashi,
               const unsigned icell,
               const vector<unsigned>& faceNodeIDs,
               map<set<unsigned>, unsigned>& faceHash2ID,
               Tessellation<3, Real>& mesh) {
  typedef set<unsigned> FaceHash;

  // Is this a new face?
  map<FaceHash, unsigned>::const_iterator faceItr = faceHash2ID.find(fhashi);
  if (faceItr == faceHash2ID.end()) {

    // Yep, it's a new face.
    const unsigned iface = mesh.faces.size();
    faceHash2ID[fhashi] = iface;
    mesh.cells[icell].push_back(iface);
    mesh.faces.push_back(faceNodeIDs);
    mesh.faceCells.push_back(vector<unsigned>());
    mesh.faceCells.back().push_back(icell);
    ASSERT(count(mesh.cells[icell].begin(), mesh.cells[icell].end(), iface) == 1);
    ASSERT(mesh.faces.size() == iface + 1);
    ASSERT(mesh.faceCells.size() == iface + 1);
    ASSERT(mesh.faceCells[iface].size() == 1 and mesh.faceCells[iface][0] == icell);

  } else {

    // Nope, this is an existing face, so we record it in the 
    // cell list as the 1's complement.
    const unsigned iface = faceItr->second;
    ASSERT(iface < mesh.faces.size());
    mesh.cells[icell].push_back(~iface);
    mesh.faceCells[iface].push_back(icell);
    ASSERT(count(mesh.cells[icell].begin(), mesh.cells[icell].end(), ~iface) == 1);
    ASSERT(mesh.faceCells[iface].size() == 2 and mesh.faceCells[iface][1] == icell);
  }
}

} // end anonymous namespace

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Real>
VoroPP_3d<Real>::
VoroPP_3d(const unsigned nx,
          const unsigned ny,
          const unsigned nz,
          const Real degeneracy):
  mNx(nx),
  mNy(ny),
  mNz(nz),
  mDegeneracy2(degeneracy*degeneracy) {
  ASSERT(mDegeneracy2 > 0.0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Real>
VoroPP_3d<Real>::
~VoroPP_3d() {
}

//------------------------------------------------------------------------------
// Compute the tessellation in the box.
//------------------------------------------------------------------------------
template<typename Real>
void
VoroPP_3d<Real>::
tessellate(vector<Real>& points,
           Real* low,
           Real* high,
           Tessellation<3, Real>& mesh) const {

  typedef set<unsigned> FaceHash;
  typedef Point3<uint64_t> VertexHash;

  const unsigned ncells = points.size()/3;
  const Real xmin = low[0], ymin = low[1], zmin = low[2];
  const Real xmax = high[0], ymax = high[1], zmax = high[2];
  const Real scale = max(xmax - xmin, max(ymax - ymin, zmax - zmin));
  const Real dx = this->degeneracy();
  const Real fconv = dx*scale;

  // Pre-conditions.
  ASSERT(xmin < xmax);
  ASSERT(ymin < ymax);
  ASSERT(zmin < zmax);
  ASSERT(scale > 0.0);
  ASSERT(points.size() % 3 == 0);
  for (int i = 0; i != ncells; ++i) {
    ASSERT(points[3*i]     >= xmin and points[3*i]     <= xmax);
    ASSERT(points[3*i + 1] >= ymin and points[3*i + 1] <= ymax);
    ASSERT(points[3*i + 2] >= zmin and points[3*i + 2] <= zmax);
  }
  ASSERT(mesh.nodes.size() == 0);
  ASSERT(mesh.cells.size() == 0);
  ASSERT(mesh.faces.size() == 0);
  ASSERT(mesh.faceCells.size() == 0);

  unsigned i, j, k, iv, iface, nf, nvf, icell;
  double xc, yc, zc;

  // Size the output arrays.
  mesh.cells.resize(ncells);

  // Scale the input coordinates to a unit box, which seems to be more robust for
  // Voro++.
  vector<Real> generators;
  generators.reserve(3*ncells);
  for (i = 0; i != ncells; ++i) {
    generators.push_back((points[3*i]     - xmin)/scale);
    generators.push_back((points[3*i + 1] - ymin)/scale);
    generators.push_back((points[3*i + 2] - zmin)/scale);
    ASSERT(generators[3*i]     >= 0.0 and generators[3*i]     <= 1.0);
    ASSERT(generators[3*i + 1] >= 0.0 and generators[3*i + 1] <= 1.0);
    ASSERT(generators[3*i + 2] >= 0.0 and generators[3*i + 2] <= 1.0);
  }
  ASSERT(generators.size() == 3*ncells);

  // Build the Voro++ container, and add the generators.
  container con(0.0, 1.0,
                0.0, 1.0,
                0.0, 1.0,
                mNx, mNy, mNz,
                false, false, false, 8);
  for (i = 0; i != ncells; ++i) con.put(i, generators[3*i], generators[3*i + 1], generators[3*i + 2]);
  ASSERT(con.total_particles() == ncells);

  // Build the tessellation cell by cell.
  voronoicell_neighbor cell;                       // Use cells with neighbor tracking.
  map<FaceHash, unsigned> faceHash2ID;             // map from face hash to mesh ID.
  map<VertexHash, unsigned> vertexHash2ID;         // map from vertex hash to mesh ID.
  c_loop_all loop(con); // Loop over all cells.
  if (loop.start()) {
    do {
      if (con.compute_cell(cell, loop)) {
        icell = loop.pid();                   // The cell index.

        // Get the cell centroid.
        double *pp = con.p[loop.ijk] + con.ps*loop.q;   // Man, this is obvious!
        cell.centroid(xc, yc, zc);
        xc += pp[0];
        yc += pp[1];
        zc += pp[2];

        // Read out the vertices into a temporary array.
        vector<Point3<uint64_t> > vertices;
        for (k = 0; k != cell.p; ++k) vertices.push_back(Point3<uint64_t>(Real(xc + 0.5*cell.pts[3*k]),
                                                                          Real(yc + 0.5*cell.pts[3*k + 1]),
                                                                          Real(zc + 0.5*cell.pts[3*k + 2]),
                                                                          dx));
        ASSERT(vertices.size() >= 4);

        // Add any new vertices from this cell to the global set, and update the vertexMap
        // to point to the global (mesh) node IDs.
        map<unsigned, unsigned> vertexMap = updateMeshVertices(vertices, vertexHash2ID, mesh, xmin, ymin, zmin, fconv);

        // // Blago!
        // std::cout << "Mesh vertices for cell " << icell << " : ";
        // std::copy(vertices.begin(), vertices.end(), ostream_iterator<Point3<uint64_t> >(std::cout, " "));
        // cout << endl;
        // std::cout << "                           ";
        // for (k = 0; k != vertices.size(); ++k) std::cout << "(" 
        //                                                  << vertices[k].realx(xmin, fconv) << " "
        //                                                  << vertices[k].realy(ymin, fconv) << " "
        //                                                  << vertices[k].realz(zmin, fconv) << ") ";

        // cout << endl;
        // std::cout << "                           ";
        // for (k = 0; k != vertices.size(); ++k) std::cout << vertexMap[k] << " ";
        // cout << endl;
        // ASSERT(vertices.size() == 8);
        // // Blago!

        // Read the face vertex indices to a temporary array as well.
        vector<int> voroFaceVertexIndices;
        cell.face_vertices(voroFaceVertexIndices);

        // Walk the faces.
        nf = cell.number_of_faces();
        ASSERT(nf >= 4);
        k = 0;
        for (iface = 0; iface != nf; ++iface) {
          ASSERT(k < voroFaceVertexIndices.size());
          FaceHash fhashi;
          vector<unsigned> faceNodeIDs;

          // Read the vertices for this face.  We assume they are listed in the
          // proper counter-clockwise order (viewed from outside)!
          nvf = voroFaceVertexIndices[k++];
          ASSERT(nvf >= 3);
          for (iv = 0; iv != nvf; ++iv) {
            ASSERT(k < voroFaceVertexIndices.size());
            ASSERT(vertexMap.find(voroFaceVertexIndices[k]) != vertexMap.end());
            j = vertexMap[voroFaceVertexIndices[k++]];

            // Is this vertex new to the face?
            if (fhashi.find(j) == fhashi.end()) {
              fhashi.insert(j);
              faceNodeIDs.push_back(j);
            }
          }
          ASSERT(faceNodeIDs.size() >= 3);

          // Add this face to the cell.
          insertFaceInfo(fhashi, icell, faceNodeIDs, faceHash2ID, mesh);
        }
      }
    } while (loop.inc());
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class VoroPP_3d<double>;
template class VoroPP_3d<float>;

}
