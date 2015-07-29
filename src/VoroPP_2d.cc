//----------------------------------------------------------------------------//
// VoroPP_2d
//----------------------------------------------------------------------------//
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>
#include <limits>
#include <stdint.h>

#include "polytope.hh"
#include "polytope_internal.hh" // Pulls in POLY_ASSERT and VoroPP_2d.hh.
#include "convexHull_2d.hh"
#include "container_2d.hh"

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
updateMeshVertices(vector<Point2<UintType> >& vertices,
                   map<Point2<UintType>, unsigned>& vertexHash2ID,
                   Tessellation<2, RealType>& mesh,
                   const RealType& xmin,
                   const RealType& ymin,
                   const RealType& fconv) {
  const unsigned n = vertices.size();
  bool newVertex;
  unsigned i, j;
  UintType ix, iy, ix0, iy0, ix1, iy1;
  Point2<UintType> ipt, ipt1;
  map<unsigned, unsigned> result;
  for (i = 0; i != n; ++i) {
    ipt = vertices[i];
    ix0 = ipt.x > 0 ? ipt.x - 1 : ipt.x;
    iy0 = ipt.y > 0 ? ipt.y - 1 : ipt.y;
    ix1 = ipt.x + 1;
    iy1 = ipt.y + 1;
    newVertex = true;
    iy = iy0;
    while (newVertex and iy != iy1) {
      ix = ix0;
      while (newVertex and ix != ix1) {
        ipt1 = Point2<UintType>(ix, iy);
        newVertex = (vertexHash2ID.find(ipt1) == vertexHash2ID.end());
        ++ix;
      }
      ++iy;
    }
    if (newVertex) {
      j = vertexHash2ID.size();
      vertexHash2ID[ipt] = j;
      mesh.nodes.push_back(vertices[i].realx(xmin, fconv));
      mesh.nodes.push_back(vertices[i].realy(ymin, fconv));
      POLY_ASSERT(mesh.nodes.size()/2 == j + 1);
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
insertFaceInfo(const pair<unsigned, unsigned>& fhashi,
               const unsigned icell,
               const unsigned inode,
               const unsigned jnode,
               map<pair<unsigned, unsigned>, unsigned>& faceHash2ID,
               Tessellation<2, RealType>& mesh) {
  typedef pair<unsigned, unsigned> FaceHash;

  // cerr << "Looking for face hash (" << fhashi.first << " " << fhashi.second << ") in {";
  // for (map<FaceHash, unsigned>::const_iterator itr = faceHash2ID.begin();
  //      itr != faceHash2ID.end();
  //      ++itr) cerr << "((" << itr->first.first << " " << itr->first.second << ") " << itr->second << ") ";
  // cerr << "}" << endl;

  // Is this a new face?
  map<FaceHash, unsigned>::const_iterator faceItr = faceHash2ID.find(fhashi);
  if (faceItr == faceHash2ID.end()) {

    // Yep, it's a new face.
    const unsigned iface = mesh.faces.size();
    faceHash2ID[fhashi] = iface;
    mesh.cells[icell].push_back(iface);
    mesh.faces.push_back(vector<unsigned>());
    mesh.faces.back().push_back(inode);
    mesh.faces.back().push_back(jnode);
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
//     cerr << iface << " " << mesh.faceCells[iface].size() << " " << mesh.faceCells[iface][0] << " " << mesh.faceCells[iface][1] << endl;
    POLY_ASSERT(mesh.faceCells[iface].size() == 2 and mesh.faceCells[iface][1] == ~int(icell));
  }
}

//------------------------------------------------------------------------------
// A unique hash for a face as an ordered collection of node indices.
//------------------------------------------------------------------------------
pair<unsigned, unsigned> 
hashFace(const unsigned i, const unsigned j) {
  POLY_ASSERT(i != j);
  return (i < j ? make_pair(i, j) : make_pair(j, i));
}

} // end anonymous namespace

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename RealType>
VoroPP_2d<RealType>::
VoroPP_2d(const unsigned nx,
          const unsigned ny,
          const RealType degeneracy):
  mNx(nx),
  mNy(ny),
  mDegeneracy2(degeneracy*degeneracy) {
  POLY_ASSERT(mDegeneracy2 > 0.0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename RealType>
VoroPP_2d<RealType>::
~VoroPP_2d() {
}

//------------------------------------------------------------------------------
// Compute the tessellation in the box.
//------------------------------------------------------------------------------
template<typename RealType>
void
VoroPP_2d<RealType>::
tessellate(const vector<RealType>& points,
           RealType* low,
           RealType* high,
           Tessellation<2, RealType>& mesh) const {

  typedef pair<unsigned, unsigned> FaceHash;
  typedef typename polytope::DimensionTraits<2, RealType>::CoordHash CoordHash;
  typedef typename polytope::DimensionTraits<2, RealType>::IntPoint VertexHash;

  const unsigned ncells = points.size()/2;
  const RealType xmin = low[0], ymin = low[1];
  const RealType xmax = high[0], ymax = high[1];
  const RealType scale = max(xmax - xmin, ymax - ymin);
  const RealType dx = this->degeneracy();
  const RealType fconv = dx*scale;

  // Pre-conditions.
  POLY_ASSERT(points.size() % 2 == 0);
  for (int i = 0; i != ncells; ++i) {
    if (!(points[2*i]     >= xmin and points[2*i]     <= xmax)) cerr << "Blago : " << points[2*i] << " " << xmin << " " << xmax << endl;
    if (!(points[2*i + 1] >= ymin and points[2*i + 1] <= ymax)) cerr << "Blago : " << points[2*i+1] << " " << ymin << " " << ymax << endl;
    POLY_ASSERT(points[2*i]     >= xmin and points[2*i]     <= xmax);
    POLY_ASSERT(points[2*i + 1] >= ymin and points[2*i + 1] <= ymax);
  }
  POLY_ASSERT(mesh.nodes.size() == 0);
  POLY_ASSERT(mesh.cells.size() == 0);
  POLY_ASSERT(mesh.faces.size() == 0);
  POLY_ASSERT(mesh.faceCells.size() == 0);
  POLY_ASSERT(xmin < xmax);
  POLY_ASSERT(ymin < ymax);
  POLY_ASSERT(scale > 0.0);

  unsigned i, j, k, nv, icell;
  double xc, yc;
  vector<CoordHash> vertices;
  vector<VertexHash > hullVertices;
  PLC<2, CoordHash> hull;

  // Size the output arrays.
  mesh.cells.resize(ncells);

  // Scale the input coordinates to a unit box, which seems to be more robust for
  // Voro++.
  vector<RealType> generators;
  generators.reserve(2*ncells);
  for (i = 0; i != ncells; ++i) {
    generators.push_back((points[2*i]     - xmin)/scale);
    generators.push_back((points[2*i + 1] - ymin)/scale);
    POLY_ASSERT(generators[2*i]     >= 0.0 and generators[2*i]     <= 1.0);
    POLY_ASSERT(generators[2*i + 1] >= 0.0 and generators[2*i + 1] <= 1.0);
  }
  POLY_ASSERT(generators.size() == 2*ncells);

  // Build the Voro++ container, and add the generators.
  container_2d con(0.0, 1.0,
                   0.0, 1.0,
                   mNx, mNy,
                   false, false, 8);
  for (i = 0; i != ncells; ++i) con.put(i, generators[2*i], generators[2*i + 1]);
  POLY_ASSERT(con.total_particles() == ncells);

  // Dummy low coordinates in integer space for the convex hull.
  CoordHash ilow[2] = {0U, 0U};

  // Build the tessellation cell by cell.
  voronoicell_neighbor_2d cell;                    // Use cells with neighbor tracking.
  map<FaceHash, unsigned> faceHash2ID;             // map from face hash to ID.
  map<VertexHash, unsigned> vertexHash2ID;         // map from vertex hash to ID.
  c_loop_all_2d loop(con); // Loop over all cells.
  if (loop.start()) {
    do {
      if (con.compute_cell(cell, loop)) {
        icell = loop.pid();                   // The cell index.

        // Get the cell centroid.
        double *pp = con.p[loop.ij] + con.ps*loop.q;  // Man, this is obvious!
        cell.centroid(xc, yc);
        xc += pp[0];
        yc += pp[1];

        // Read out the vertices into a temporary array.
        nv = cell.p;
        POLY_ASSERT(nv >= 3);
        vertices = vector<CoordHash>();
        vertices.reserve(2*nv);
        for (k = 0; k != nv; ++k) {
          vertices.push_back(CoordHash(max(0.0, min(numeric_limits<CoordHash>::max() - 1.0, (xc + 0.5*cell.pts[2*k    ])/dx + 0.5))));
          vertices.push_back(CoordHash(max(0.0, min(numeric_limits<CoordHash>::max() - 1.0, (yc + 0.5*cell.pts[2*k + 1])/dx + 0.5))));
        }

        // Sort the vertices counter-clockwise.
        hull = convexHull_2d(vertices, ilow, CoordHash(1));
//         if (!(hull.facets.size() == nv)) {
//           cerr << "Input vertices: " << endl;
//           for (unsigned k = 0; k != nv; ++k) cerr << "    " << k << " (" << vertices[2*k] << " " << vertices[2*k + 1] << ")" << endl;
//           cerr << "Hull ordering: ";
//           for (unsigned k = 0; k != hull.facets.size(); ++k) cerr << " " << hull.facets[k][0];
//           cerr << endl;
//         }
        POLY_ASSERT(hull.facets.size() <= nv);            // There may be some degenerate vertices, but the hull method reduces to the unique set.
        nv = hull.facets.size();
        hullVertices = vector<VertexHash >();
        hullVertices.reserve(nv);
        for (k = 0; k != nv; ++k) {
          i = hull.facets[k][0];
          hullVertices.push_back(VertexHash(vertices[2*i], vertices[2*i + 1]));
        }

        // Add any new vertices from this cell to the global set, and update the vertexMap
        // to point to the global (mesh) node IDs.
        map<unsigned, unsigned> vertexMap = updateMeshVertices(hullVertices, vertexHash2ID, mesh, xmin, ymin, fconv);
//         cerr << "Cell vertex mapping: " << endl;
//         for (unsigned k = 0; k != nv; ++k) {
//           cerr << "  --> " << hullVertices[k] << " " << vertexMap[k] << endl;
//         }

        // Build the faces by walking the cell vertices counter-clockwise.
        for (k = 0; k != nv; ++k) {
          i = vertexMap[k];
          j = vertexMap[(k + 1) % nv];
          POLY_ASSERT(i < mesh.nodes.size()/2);
          POLY_ASSERT(j < mesh.nodes.size()/2);
          POLY_ASSERT(i != j);
          insertFaceInfo(hashFace(i, j), icell, i, j, faceHash2ID, mesh);
        }
      }
    } while (loop.inc());
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class VoroPP_2d<double>;
//template class VoroPP_2d<float>;

}
