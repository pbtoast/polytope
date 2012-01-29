//---------------------------------Spheral++----------------------------------//
// VoroPP
//----------------------------------------------------------------------------//
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>

#include "polytope.hh" // Pulls in ASSERT and VoroPP_2d.hh.
#include "container_2d.hh"

namespace polytope {

using namespace std;
using namespace voro;
using std::min;
using std::max;
using std::abs;

namespace { // We hide internal functions in an anonymous namespace.

//------------------------------------------------------------------------------
// A simple 2D point.
//------------------------------------------------------------------------------
template<typename Real>
struct Point2 {
  Real x, y;
  Point2(): x(0.0), y(0.0) {}
  Point2(const Real& xi, const Real& yi): x(xi), y(yi) {}
  Point2& operator=(const Point2& rhs) { x = rhs.x; y = rhs.y; return *this; }
  bool operator<(const Point2<Real>& rhs) const {
    return ((x < rhs.x) or (x == rhs.x and y < rhs.y));
  }
};

// It's nice being able to print these things.
template<typename Real>
std::ostream&
operator<<(std::ostream& os, const Point2<Real>& p) {
  os << "(" << p.x << " " << p.y << ")";
  return os;
}

//------------------------------------------------------------------------------
// Z coordinate of cross product : (p2 - p1)x(p3 - p1).
//------------------------------------------------------------------------------
template<typename Real>
double zcross(const Point2<Real>& p1, const Point2<Real>& p2, const Point2<Real>& p3) {
  return (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);
}

//------------------------------------------------------------------------------
// Sort a set of Point2 points in counter-clockwise order using Andrews 
// monotone chain algorithm.
// Based on an example at http://www.algorithmist.com/index.php/Monotone_Chain_Convex_Hull.cpp
//------------------------------------------------------------------------------
template<typename Real>
vector<Point2<Real> >
sortCounterClockwise(vector<Point2<Real> >& points) {
  const unsigned n = points.size();
  int i, k, t;
  
  // Sort the input points by x coordinate.
  sort(points.begin(), points.end());

  // Prepare the result.
  vector<Point2<Real> > result(2*n);

  // Build the lower hull.
  for (i = 0, k = 0; i < n; i++) {
    while (k >= 2 and zcross(result[k - 2], result[k - 1], points[i]) <= 0.0) k--;
    result[k++] = points[i];
  }

  // Build the upper hull.
  for (i = n - 2, t = k + 1; i >= 0; i--) {
    while (k >= t and zcross(result[k - 2], result[k - 1], points[i]) <= 0.0) k--;
    result[k++] = points[i];
  }

  // Size the result and we're done.
  result.resize(n);
  return result;
}

//------------------------------------------------------------------------------
// Helper method to update our face info.
//------------------------------------------------------------------------------
template<typename Real>
inline
void
insertFaceInfo(const pair<unsigned, unsigned>& fhashi,
               const unsigned icell,
               const unsigned inode,
               const unsigned jnode,
               map<pair<unsigned, unsigned>, unsigned>& faceHash2ID,
               Tessellation<2, Real>& mesh) {
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
    // cerr << iface << " " << mesh.faceCells[iface].size() << " " << mesh.faceCells[iface][0] << " " << mesh.faceCells[iface][1] << endl;
    ASSERT(mesh.faceCells[iface].size() == 2 and mesh.faceCells[iface][1] == icell);
  }
}

//------------------------------------------------------------------------------
// Compute the distance^2 between points.
//------------------------------------------------------------------------------
template<typename Real>
inline
Real
distance2(const Real& x1, const Real& y1,
          const Real& x2, const Real& y2) {
  return (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 -y1);
}

//------------------------------------------------------------------------------
// A unique hash for a face as an ordered collection of node indices.
//------------------------------------------------------------------------------
pair<unsigned, unsigned> 
hashFace(const unsigned i, const unsigned j)
{
  ASSERT(i != j);
  return (i < j ? make_pair(i, j) : make_pair(j, i));
}

} // end anonymous namespace

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Real>
VoroPP_2d<Real>::
VoroPP_2d(const Real xmin, const Real ymin,
          const Real xmax, const Real ymax,
          const unsigned nx,
          const unsigned ny,
          const Real degeneracy):
  mNx(nx),
  mNy(ny),
  mxmin(xmin),
  mymin(ymin),
  mxmax(xmax),
  mymax(ymax),
  mDegeneracy2(degeneracy*degeneracy),
  mScale(max(xmax - xmin, ymax - ymin)) {
  ASSERT(mxmin < mxmax);
  ASSERT(mymin < mymax);
  ASSERT(mDegeneracy2 > 0.0);
  ASSERT(mScale > 0.0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Real>
VoroPP_2d<Real>::
~VoroPP_2d() {
}

//------------------------------------------------------------------------------
// Compute the tessellation in the box.
//------------------------------------------------------------------------------
template<typename Real>
void
VoroPP_2d<Real>::
tessellate(const vector<Real>& points,
           Tessellation<2, Real>& mesh) const {

  typedef pair<unsigned, unsigned> FaceHash;

  // Pre-conditions.
  ASSERT(points.size() % 2 == 0);
  for (int i = 0; i != points.size(); ++i) {
    if (!(points[2*i]     >= mxmin and points[2*i]     <= mxmax)) cerr << "Blago : " << points[2*i] << " " << mxmin << " " << mxmax << endl;
    if (!(points[2*i + 1] >= mymin and points[2*i + 1] <= mymax)) cerr << "Blago : " << points[2*i+1] << " " << mymin << " " << mymax << endl;
    ASSERT(points[2*i]     >= mxmin and points[2*i]     <= mxmax);
    ASSERT(points[2*i + 1] >= mymin and points[2*i + 1] <= mymax);
  }
  ASSERT(mesh.nodes.size() == 0);
  ASSERT(mesh.cells.size() == 0);
  ASSERT(mesh.faces.size() == 0);
  ASSERT(mesh.faceCells.size() == 0);

  bool newNode;
  const unsigned ncells = points.size()/2;
  unsigned i, n, icell, jcell;
  double xc, yc;
  Real xv, yv, xv_last, yv_last;

  // Size the output arrays.
  mesh.cells.resize(ncells);

  // Scale the input coordinates to a unit box, which seems to be more robust for
  // Voro++.
  vector<Real> generators;
  generators.reserve(2*ncells);
  for (i = 0; i != ncells; ++i) {
    generators.push_back((points[2*i]     - mxmin)/mScale);
    generators.push_back((points[2*i + 1] - mymin)/mScale);
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
  vector<vector<unsigned> > cellNeighbors(ncells); // Keep track of neighbor cells.
  vector<vector<unsigned> > cellNodes(ncells);     // Keep track of the cell nodes.
  map<FaceHash, unsigned> faceHash2ID;             // map from face hash to ID.
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

        // Read the neighbor cell IDs.  Any negative IDs indicate a boundary
        // surface, so just throw them away.
        vector<int> tmpNeighbors;
        cell.neighbors(tmpNeighbors);
        remove_copy_if(tmpNeighbors.begin(), tmpNeighbors.end(), 
                       back_inserter(cellNeighbors[icell]),
                       bind2nd(less<int>(), 0));

        // Read out the vertices into a temporary array.
        vector<Point2<Real> > vertices;
        for (unsigned k = 0; k != cell.p; ++k) vertices.push_back(Point2<Real>(xc + 0.5*cell.pts[2*k],
                                                                               yc + 0.5*cell.pts[2*k + 1]));
        ASSERT(vertices.size() >= 3);

        // Sort the vertices counter-clockwise.
        vertices = sortCounterClockwise(vertices);

        // Assign the global nodes based on the cell vertices.
        xv_last = 10.0;
        yv_last = 10.0;
        for (unsigned k = 0; k != vertices.size(); ++k) {

          // Vertex position.
          xv = vertices[k].x;
          yv = vertices[k].y;
          // xv = xc + 0.5*cell.pts[2*k];
          // yv = yc + 0.5*cell.pts[2*k + 1];

          // Is this node distinct from the last one we visited?
          if (distance2(xv, yv, xv_last, yv_last) > mDegeneracy2) {

            // Has this vertex already been created by one of our neighbors?
            newNode = true;
            vector<unsigned>::const_iterator neighborItr = cellNeighbors[icell].begin();
            while (newNode and neighborItr != cellNeighbors[icell].end()) {
              jcell = *neighborItr++;
              ASSERT(jcell < ncells);
              vector<unsigned>::const_iterator nodeItr = cellNodes[jcell].begin();
              while (newNode and nodeItr != cellNodes[jcell].end()) {
                i = *nodeItr++;
                ASSERT(i < mesh.nodes.size());
                if (distance2(xv, yv, mesh.nodes[2*i], mesh.nodes[2*i + 1]) < mDegeneracy2) {
                  // Found it!
                  newNode = false;
                  cellNodes[icell].push_back(i);
                }
              }
            }

            // This is a new vertex position, so create a new node.
            if (newNode) {
              i = mesh.nodes.size()/2;
              cellNodes[icell].push_back(i);
              mesh.nodes.push_back(xv);
              mesh.nodes.push_back(yv);
            }

            // Figure out what "face" this node and previous one in the cell represent.
            n = cellNodes[icell].size();
            if (n > 1) {
              int i = cellNodes[icell][n - 2];
              int j = cellNodes[icell][n - 1];
              insertFaceInfo(hashFace(i, j), icell, i, j, faceHash2ID, mesh);
            }
          }
          xv_last = xv;
          yv_last = yv;
        }
        ASSERT(cellNodes[icell].size() >= 3);

        // We have to tie together the first and last cell vertices in a final face.
        int i = cellNodes[icell].back();
        int j = cellNodes[icell].front();
        insertFaceInfo(hashFace(i, j), icell, i, j, faceHash2ID, mesh);
      }
    } while (loop.inc());
  }
        
  // De-normalize the vertex coordinates back to the input frame.
  for (int i = 0; i != mesh.nodes.size(); ++i) mesh.nodes[i] = mxmin + mScale*mesh.nodes[i];
}

//------------------------------------------------------------------------------
// Compute the tessellation in the box.
//------------------------------------------------------------------------------
template<typename Real>
void
VoroPP_2d<Real>::
tessellate(const vector<Real>& points,
           const PLC<2, Real>& geometry,
           Tessellation<2, Real>& mesh) const {
  ASSERT(false); // Implemenet me!
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class VoroPP_2d<double>;
template class VoroPP_2d<float>;

}
