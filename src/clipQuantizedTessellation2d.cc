//------------------------------------------------------------------------------
// 2D implementation of clipQuantizedTessellation.
//
// This method relies upon the Boost.Polygon library for geometric operations.
//------------------------------------------------------------------------------

#include <vector>
#include <map>

#include "clipQuantizedTessellation.hh"
#include "RegisterBoostPolygonTypes.hh"

using namespace std;
namespace bp = boost::polygon;
using namespace boost::polygon::operators;

namespace polytope {

template<typename IntType, typename RealType>
void clipQuantizedTessellation(QuantizedTessellation2d<IntType, RealType>& qmesh,
                               const std::vector<RealType>& PLCpoints,
                               const PLC<2, RealType>& geometry) {

  using namespace bp;
  typedef polygon_with_holes_data<IntType> Polygon;
  typedef typename polygon_traits<Polygon>::point_type Point;
  typedef std::vector<polygon_data<IntType> > PolygonSet;

  // Copy the input PLC geometry to a polygon Boost.Polygon will recognize.
  // We have registered a std::vector<IntPoint> as a polygon.
  const unsigned numBoundaryPoints = geometry.facets.size();
  std::vector<Point> boundaryPoints(numBoundaryPoints);
  IntPoint p;
  for (unsigned i = 0; i != geometry.facets.size(); ++i) {
    POLY_ASSERT(geometry.facets[i].size() == 2);
    const int k = geometry.facets[i][0];
    qmesh.quantize(&PLCpoints[2*k], &p.x);
    boundaryPoints[i] = bp::construct<Point>(p.x, p.y);
  }
  Polygon boundary;
  bp::set_points(boundary, boundaryPoints.begin(), boundaryPoints.end());

  // Add the holes.
  std::vector<std::vector<Point> > holePoints(geometry.holes.size());
  for (unsigned ihole = 0; ihole != geometry.holes.size(); ++ihole) {
    const unsigned n = geometry.holes[ihole].size();
    holePoints[ihole].resize(n);
    for (unsigned i = 0; i != n; ++i) {
      POLY_ASSERT(geometry.holes[ihole][i].size() == 2);
      const int k = geometry.holes[ihole][i][0];
      qmesh.quantize(&PLCpoints[2*k], &p.x);
      holePoints[ihole][i] = bp::construct<Point>(p.x, p.y);
    }
  }
  boundary.set_holes(holePoints.begin(), holePoints.end());
  
  // Make the polygon_set boundary.
  PolygonSet boundarySet;
  boundarySet += boundary;

  // Prepare new mesh data.
  const unsigned ncells = qmesh.cellEdges.size();
  std::vector<bp::IntPoint> newNodes;
  std::vector<std::pair<int, int> > newEdges;
  std::vector<std::vector<int> > newCellEdges(ncells);
  std::map<bp::IntPoint, int, PointComparator<IntType> > node2id(PointComparator<IntType>(2));
  std::map<std::pair<int, int>, int> edge2id;

  // We also prepare to keep track of the orphans for each cell.
  std::map<unsigned, PolygonSet> orphans;

  // Now walk each cell and clip it with the boundary.
  for (unsigned i = 0; i != ncells; ++i) {
    unsigned nverts = qmesh.cellEdges[i].size();

    // Copy the cell geometry to the quantized polygon type.
    Polygon cell;
    std::vector<Point> cellPoints(nverts);
    for (unsigned j = 0; j != nverts; ++j) {
      int k = qmesh.cellEdges[i][j];
      if (k < 0) {
        k = ~k;
        p = qmesh.nodes[qmesh.edges[k].second];
      } else {
        p = qmesh.nodes[qmesh.edges[k].first];
      }
      cellPoints[j] = bp::construct<Point>(p.x, p.y);
    }
    bp::set_points(cell, cellPoints.begin(), cellPoints.end());

    // Clip the cell against the boundary.
    PolygonSet cellSet;
    cellSet += cell;
    cellSet &= boundarySet;
    unsigned polygonIndex = 0;
    if (cellSet.size() > 1) {
      // This intersection generated orphans.  Find the polygon that contains the generator -- that becomes this new cell.
      // The other fragments we hold onto for adoption.
      Point gen = bp::construct<Point>(qmesh.generators[i].x, qmesh.generators[i].y);
      while (polygonIndex < cellSet.size() and not bp::contains(cellSet[polygonIndex], gen)) ++polygonIndex;
      POLY_ASSERT(polygonIndex < cellSet.size());
      orphans[i] = cellSet;
      orphans[i].erase(orphans[i].begin() + polygonIndex);
      std::cerr << "polytope clipping WARNING: detected " << (cellSet.size() - 1) << " orphan cell piece(s): identified generator in fragment:" << polygonIndex << std::endl;
    }

    // Read out the final cell geometry to the new QuantizedTessellation.
    // It appears Boost.Polygon gives us back the same node for the beginning and ending nodes of the
    // polygon.
    nverts = cellSet[polygonIndex].size();
    POLY_ASSERT(nverts > 3);
    POLY_ASSERT(*(cellSet[0].begin()) == *(cellSet[0].begin() + nverts - 1));
    for (unsigned j = 0; j != nverts-1; ++j) {
      const Point& v0 = *(cellSet[0].begin() + j);
      const Point& v1 = *(cellSet[0].begin() + j + 1);
   
      // {
      //   const double mag2 = (double(v0.x()) - double(v1.x()))*(double(v0.x()) - double(v1.x())) + (double(v0.y()) - double(v1.y()))*(double(v0.y()) - double(v1.y()));
      //   if (mag2 < 5) {
      //     std::cerr << " Strange point: " << i << " " << nverts 
      //               << " (" << v0.x() << " " << v0.y() << ")" 
      //               << " (" << v1.x() << " " << v1.y() << ")"
      //               << std::endl;
      //   }
      // }

      // Insert vertex 0.
      p.x = v0.x();
      p.y = v0.y();
      int old_size = node2id.size();
      const int j0 = internal::addKeyToMap(p, node2id);
      if (j0 == old_size) {
        POLY_ASSERT(j0 == newNodes.size());
        newNodes.push_back(p);
      }

      // Insert vertex 1.
      p.x = v1.x();
      p.y = v1.y();
      old_size = node2id.size();
      const int j1 = internal::addKeyToMap(p, node2id);
      if (j1 == old_size) {
        POLY_ASSERT(j1 == newNodes.size());
        newNodes.push_back(p);
      }
    
      // Now insert the edge if non-degenerate.
      if (j0 != j1) {
        const std::pair<int, int> edge = internal::hashEdge(j0, j1);
        POLY_ASSERT((edge.first == j0 and edge.second == j1) or
                    (edge.first == j1 and edge.second == j0));
        old_size = edge2id.size();
        const int e1 = internal::addKeyToMap(edge, edge2id);
        if (e1 == old_size) {
          POLY_ASSERT(e1 == newEdges.size());
          newEdges.push_back(edge);
        }
        if (edge.first == j0) {
          newCellEdges[i].push_back(e1);
        } else {
          newCellEdges[i].push_back(~e1);
        }
      }
    }
  }

  // Copy the topology and geometry back to the input QuantizedTessellation.
  qmesh.nodes = newNodes;
  qmesh.edges = newEdges;
  qmesh.cellEdges = newCellEdges;
}

//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
template
void
clipQuantizedTessellation<int, double>(QuantizedTessellation2d<int, double>& qmesh,
                                       const std::vector<double>& PLCpoints,
                                       const PLC<2, double>& geometry);

}
