//----------------------------------------------------------------------------//
// ParallelTessellator
//----------------------------------------------------------------------------//
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>

#include "polytope.hh" // Pulls in ASSERT and VoroPP_2d.hh.
#include "convexHull_2d.hh"

namespace polytope {

using namespace std;
using std::min;
using std::max;
using std::abs;

namespace { // We hide internal functions in an anonymous namespace.


} // end anonymous namespace

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
ParallelTessellator<Dimension, RealType>::
ParallelTessellator(const Tessellator<Dimension, RealType>& tessellator):
  mSerialTessellator(tessellator) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
ParallelTessellator<Dimension, RealType>::
~ParallelTessellator() {
}

//------------------------------------------------------------------------------
// Compute the tessellation in the box.
//------------------------------------------------------------------------------
template<typename int, typename RealType>
void
ParallelTessellator<Dimension, RealType>::
tessellate(vector<RealType>& points,
           RealType* low,
           RealType* high,
           Tessellation<Dimension, RealType>& mesh) const {

  const unsigned ncells = points.size()/Dimension;


  const RealType xmin = low[0], ymin = low[1];
  const RealType xmax = high[0], ymax = high[1];
  const RealType scale = max(xmax - xmin, ymax - ymin);
  const RealType dx = this->degeneracy();
  const RealType fconv = dx*scale;

  // Pre-conditions.
  ASSERT(points.size() % 2 == 0);
  for (int i = 0; i != ncells; ++i) {
    if (!(points[2*i]     >= xmin and points[2*i]     <= xmax)) cerr << "Blago : " << points[2*i] << " " << xmin << " " << xmax << endl;
    if (!(points[2*i + 1] >= ymin and points[2*i + 1] <= ymax)) cerr << "Blago : " << points[2*i+1] << " " << ymin << " " << ymax << endl;
    ASSERT(points[2*i]     >= xmin and points[2*i]     <= xmax);
    ASSERT(points[2*i + 1] >= ymin and points[2*i + 1] <= ymax);
  }
  ASSERT(mesh.nodes.size() == 0);
  ASSERT(mesh.cells.size() == 0);
  ASSERT(mesh.faces.size() == 0);
  ASSERT(mesh.faceCells.size() == 0);
  ASSERT(xmin < xmax);
  ASSERT(ymin < ymax);
  ASSERT(scale > 0.0);

  unsigned i, j, k, nv, icell;
  double xc, yc;

  // Size the output arrays.
  mesh.cells.resize(ncells);

  // Scale the input coordinates to a unit box, which seems to be more robust for
  // Voro++.
  vector<RealType> generators;
  generators.reserve(2*ncells);
  for (i = 0; i != ncells; ++i) {
    generators.push_back((points[2*i]     - xmin)/scale);
    generators.push_back((points[2*i + 1] - ymin)/scale);
    ASSERT(generators[2*i]     >= 0.0 and generators[2*i]     <= 1.0);
    ASSERT(generators[2*i + 1] >= 0.0 and generators[2*i + 1] <= 1.0);
  }
  ASSERT(generators.size() == 2*ncells);

  // Build the Voro++ container, and add the generators.
  container_2d con(0.0, 1.0,
                   0.0, 1.0,
                   mNx, mNy,
                   false, false, 8);
  for (i = 0; i != ncells; ++i) con.put(i, generators[2*i], generators[2*i + 1]);
  ASSERT(con.total_particles() == ncells);

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
        vector<Point2<uint64_t> > vertices;
        for (k = 0; k != cell.p; ++k) vertices.push_back(Point2<uint64_t>(RealType(xc + 0.5*cell.pts[2*k]),
                                                                          RealType(yc + 0.5*cell.pts[2*k + 1]),
                                                                          dx));
        ASSERT(vertices.size() >= 3);

        // Sort the vertices counter-clockwise.
        vertices = convexHull_2d(vertices);

        // Add any new vertices from this cell to the global set, and update the vertexMap
        // to point to the global (mesh) node IDs.
        map<unsigned, unsigned> vertexMap = updateMeshVertices(vertices, vertexHash2ID, mesh, xmin, ymin, fconv);

        // Build the faces by walking the cell vertices counter-clockwise.
        nv = vertices.size();
        for (k = 0; k != nv; ++k) {
          i = vertexMap[k];
          j = vertexMap[(k + 1) % nv];
          ASSERT(i < mesh.nodes.size()/2);
          ASSERT(j < mesh.nodes.size()/2);

          // If these vertices are distinct, add the face.
          if (i != j) insertFaceInfo(hashFace(i, j), icell, i, j, faceHash2ID, mesh);
        }
      }
    } while (loop.inc());
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class ParallelTessellator<double>;

}
