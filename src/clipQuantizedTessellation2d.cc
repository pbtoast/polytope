//------------------------------------------------------------------------------
// 2D implementation of clipQuantizedTessellation.
//
// This method relies upon the Boost.Polygon library for geometric operations.
//------------------------------------------------------------------------------

#include <vector>
#include <map>

#include "clipQuantizedTessellation.hh"
#include "RegisterBoostPolygonTypes.hh"

namespace bp = boost::polygon;
using namespace boost::polygon::operators;

namespace polytope {

template<typename IntType, typename RealType>
void clipQuantizedTessellation(QuantizedTessellation2d<IntType, RealType>& qmesh,
                               const std::vector<RealType>& PLCpoints,
                               const PLC<2, RealType>& geometry) {

  using namespace bp;
  typedef polygon_data<CoordHash> Polygon;
  typedef polygon_traits<Polygon>::point_type Point;
  typedef std::vector<polygon_data<CoordHash> > PolygonSet;

  // Copy the input PLC geometry to a polygon Boost.Polygon will recognize.
  // We have registered a std::vector<IntPoint> as a polygon.
  const unsigned numBoundaryPoints = geometry.facets.size();
  std::vector<Point> boundaryPoints(numBoundaryPoints);
  IntPoint p;
  for (unsigned i = 0; i != geometry.facets.size(); ++i) {
    POLY_ASSERT(geometry.facets[i].size() == 2);
    qmesh.quantize(&PLCpoints[2*i], &p.x);
    boundaryPoints[i] = bp::construct<Point>(p.x, p.y);
  }
  Polygon boundary;
  bp::set_points(boundary, boundaryPoints.begin(), boundaryPoints.end());
  PolygonSet boundarySet;
  boundarySet += boundary;

  // TODO -- handle holes.
  
  // Prepare new mesh data.
  const unsigned ncells = qmesh.cellEdges.size();
  std::vector<bp::IntPoint> newNodes;
  std::vector<std::pair<int, int> > newEdges;
  std::vector<std::vector<int> > newCellEdges(ncells);
  std::map<bp::IntPoint, int> node2id;
  std::map<std::pair<int, int>, int> edge2id;

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
    if (cellSet.size() > 1) {
      std::cerr << "polytope clipping WARNING: detected " << (cellSet.size() - 1) << " orphan cell piece(s), ignoring." << std::endl;
    }

    // Read out the final cell geometry to the new QuantizedTessellation.
    // It appears Boost.Polygon gives us back the same node for the beginning and ending nodes of the
    // polygon.
    nverts = cellSet[0].size();
    POLY_ASSERT(nverts > 3);
    for (unsigned j = 0; j != nverts-1; ++j) {
      const Point& v0 = *(cellSet[0].begin() + j);
      const Point& v1 = *(cellSet[0].begin() + j + 1);
   
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
    
      // Now insert the edge.
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
