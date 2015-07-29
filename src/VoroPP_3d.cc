//----------------------------------------------------------------------------//
// VoroPP_3d
//----------------------------------------------------------------------------//
#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>
#include <set>

#include "polytope.hh"
#include "polytope_internal.hh" // Pulls in POLY_ASSERT and VoroPP_3d.hh.
#include "Point.hh"
#include "container.hh"

namespace polytope {

using namespace std;
using namespace voro;
using std::min;
using std::max;
using std::abs;

namespace { // We hide internal functions in an anonymous namespace.

//------------------------------------------------------------------------------
// Take the given set of vertex positions, and either find them in the known
// mesh nodes or add them to the known set.
//------------------------------------------------------------------------------
template<typename RealType, typename UintType>
map<unsigned, unsigned>
updateMeshVertices(vector<Point3<UintType> >& vertices,
                   map<Point3<UintType>, unsigned>& vertexHash2ID,
                   Tessellation<3, RealType>& mesh,
                   const RealType& xmin,
                   const RealType& ymin,
                   const RealType& zmin,
                   const RealType& fconv) {
  const unsigned n = vertices.size();
  bool newVertex;
  unsigned i, j;
  UintType ix, iy, iz, ix0, iy0, iz0, ix1, iy1, iz1;
  Point3<UintType> ipt, ipt1;
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
          ipt1 = Point3<UintType>(ix, iy, iz);
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
      POLY_ASSERT(mesh.nodes.size()/3 == j + 1);
      result[i] = j;
    } else {
      result[i] = vertexHash2ID[ipt1];
    }
  }
  POLY_ASSERT(result.size() == vertices.size());
  return result;
}

//------------------------------------------------------------------------------
// Helper method to update our face info.
//------------------------------------------------------------------------------
template<typename RealType>
inline
void
insertFaceInfo(const set<unsigned>& fhashi,
               const unsigned icell,
               const vector<unsigned>& faceNodeIDs,
               map<set<unsigned>, unsigned>& faceHash2ID,
               Tessellation<3, RealType>& mesh) {
  typedef set<unsigned> FaceHash;

  // Is this a new face?
  map<FaceHash, unsigned>::const_iterator faceItr = faceHash2ID.find(fhashi);
  if (faceItr == faceHash2ID.end()) {

    // Yep, it's a new face.
    const unsigned iface = mesh.faces.size();
    faceHash2ID[fhashi] = iface;
    mesh.cells[icell].push_back(iface);
    mesh.faces.push_back(faceNodeIDs);
    mesh.faceCells.push_back(vector<int>());
    mesh.faceCells.back().push_back(icell);
    POLY_ASSERT(count(mesh.cells[icell].begin(), mesh.cells[icell].end(), iface) == 1);
    POLY_ASSERT(mesh.faces.size() == iface + 1);
    POLY_ASSERT(mesh.faceCells.size() == iface + 1);
    POLY_ASSERT(mesh.faceCells[iface].size() == 1 and mesh.faceCells[iface][0] == icell);

  } else {

    // Nope, this is an existing face, so we record it in the 
    // cell list as the 1's complement.
    const int iface = faceItr->second;
    POLY_ASSERT(iface < mesh.faces.size());
    mesh.cells[icell].push_back(~iface);
    mesh.faceCells[iface].push_back(~int(icell));
    POLY_ASSERT(count(mesh.cells[icell].begin(), mesh.cells[icell].end(), ~iface) == 1);
    POLY_ASSERT(mesh.faceCells[iface].size() == 2 and mesh.faceCells[iface][1] == ~int(icell));
  }
}

//------------------------------------------------------------------------------
// Define our own local random number generator wrapping the standard srand &
// rand methods.
//------------------------------------------------------------------------------
double random01() {
  return double(rand())/RAND_MAX;
}

} // end anonymous namespace

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename RealType>
VoroPP_3d<RealType>::
VoroPP_3d(const unsigned nx,
          const unsigned ny,
          const unsigned nz,
          const RealType degeneracy):
  mNx(nx),
  mNy(ny),
  mNz(nz),
  mDegeneracy2(degeneracy*degeneracy) {
  POLY_ASSERT(mDegeneracy2 > 0.0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename RealType>
VoroPP_3d<RealType>::
~VoroPP_3d() {
}

//------------------------------------------------------------------------------
// Compute the tessellation in the box.
//------------------------------------------------------------------------------
template<typename RealType>
void
VoroPP_3d<RealType>::
tessellate(const vector<RealType>& points,
           RealType* low,
           RealType* high,
           Tessellation<3, RealType>& mesh) const {

  typedef set<unsigned> FaceHash;
  typedef typename polytope::DimensionTraits<3, RealType>::IntPoint VertexHash;

  const unsigned ncells = points.size()/3;
  const RealType xmin = low[0], ymin = low[1], zmin = low[2];
  const RealType xmax = high[0], ymax = high[1], zmax = high[2];
  const RealType scale = max(xmax - xmin, max(ymax - ymin, zmax - zmin));
  const RealType dx = this->degeneracy();
  const RealType fconv = dx*scale;
  const RealType randomfuzz = 0.0; // 2.0*dx;

  // Pre-conditions.
  POLY_ASSERT(xmin < xmax);
  POLY_ASSERT(ymin < ymax);
  POLY_ASSERT(zmin < zmax);
  POLY_ASSERT(scale > 0.0);
  POLY_ASSERT(points.size() % 3 == 0);
  for (int i = 0; i != ncells; ++i) {
    POLY_ASSERT(points[3*i]     >= xmin and points[3*i]     <= xmax);
    POLY_ASSERT(points[3*i + 1] >= ymin and points[3*i + 1] <= ymax);
    POLY_ASSERT(points[3*i + 2] >= zmin and points[3*i + 2] <= zmax);
  }
  POLY_ASSERT(mesh.nodes.size() == 0);
  POLY_ASSERT(mesh.cells.size() == 0);
  POLY_ASSERT(mesh.faces.size() == 0);
  POLY_ASSERT(mesh.faceCells.size() == 0);

  unsigned i, j, k, iv, iface, nf, nvf, icell;
  double xc, yc, zc;

  // Size the output arrays.
  mesh.cells.resize(ncells);

  // Scale the input coordinates to a unit box, which seems to be more robust for
  // Voro++.
  vector<RealType> generators;
  generators.reserve(3*ncells);
  for (i = 0; i != ncells; ++i) {
    generators.push_back((points[3*i]     - xmin)/scale);
    generators.push_back((points[3*i + 1] - ymin)/scale);
    generators.push_back((points[3*i + 2] - zmin)/scale);
    POLY_ASSERT(generators[3*i]     >= 0.0 and generators[3*i]     <= 1.0);
    POLY_ASSERT(generators[3*i + 1] >= 0.0 and generators[3*i + 1] <= 1.0);
    POLY_ASSERT(generators[3*i + 2] >= 0.0 and generators[3*i + 2] <= 1.0);
  }
  POLY_ASSERT(generators.size() == 3*ncells);

  // Build the Voro++ container, and add the generators.
  container con(0.0, 1.0,
                0.0, 1.0,
                0.0, 1.0,
                mNx, mNy, mNz,
                false, false, false, 8);
  for (i = 0; i != ncells; ++i) con.put(i, 
                                        max(0.0, min(1.0, generators[3*i]     + (random01() - 0.5)*randomfuzz)),
                                        max(0.0, min(1.0, generators[3*i + 1] + (random01() - 0.5)*randomfuzz)),
                                        max(0.0, min(1.0, generators[3*i + 2] + (random01() - 0.5)*randomfuzz)));
  POLY_ASSERT(con.total_particles() == ncells);

  // Build the tessellation cell by cell.
  voronoicell cell;                                // The Voro++ cell (without neighbor tracking).
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
        vector<VertexHash > vertices;
        for (k = 0; k != cell.p; ++k) vertices.push_back(VertexHash(RealType(xc + 0.5*cell.pts[3*k]),
								    RealType(yc + 0.5*cell.pts[3*k + 1]),
								    RealType(zc + 0.5*cell.pts[3*k + 2]),
								    dx));
        POLY_ASSERT(vertices.size() >= 4);

        // Add any new vertices from this cell to the global set, and update the vertexMap
        // to point to the global (mesh) node IDs.
        map<unsigned, unsigned> vertexMap = updateMeshVertices(vertices, vertexHash2ID, mesh, xmin, ymin, zmin, fconv);

        // // Blago!
        // std::cout << "Mesh vertices for cell " << icell << " : ";
        // std::copy(vertices.begin(), vertices.end(), ostream_iterator<VertexHash >(std::cout, " "));
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
        // POLY_ASSERT(vertices.size() == 8);
        // // Blago!

        // Read the face vertex indices to a temporary array as well.
        vector<int> voroFaceVertexIndices;
        cell.face_vertices(voroFaceVertexIndices);

        // Walk the faces.
        nf = cell.number_of_faces();
        POLY_ASSERT(nf >= 4);
        k = 0;
        for (iface = 0; iface != nf; ++iface) {
          POLY_ASSERT(k < voroFaceVertexIndices.size());
          FaceHash fhashi;
          vector<unsigned> faceNodeIDs;

          // Read the vertices for this face.  We assume they are listed in the
          // proper counter-clockwise order (viewed from outside)!
          nvf = voroFaceVertexIndices[k++];
          POLY_ASSERT(nvf >= 3);
          for (iv = 0; iv != nvf; ++iv) {
            POLY_ASSERT(k < voroFaceVertexIndices.size());
            POLY_ASSERT(vertexMap.find(voroFaceVertexIndices[k]) != vertexMap.end());
            j = vertexMap[voroFaceVertexIndices[k++]];

            // Is this vertex new to the face?
            if (fhashi.find(j) == fhashi.end()) {
              fhashi.insert(j);
              faceNodeIDs.push_back(j);
            }
          }
          POLY_ASSERT(faceNodeIDs.size() >= 3);

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
// template class VoroPP_3d<float>;

}
