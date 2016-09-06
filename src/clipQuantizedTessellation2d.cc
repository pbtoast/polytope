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

  // Copy the input PLC geometry to a polygon Boost.Polygon will recognize.
  // We have registered a std::vector<IntPoint> as a polygon.
  bp::IntPolygon boundary;
  boundary.resize(geometry.facets.size());
  for (unsigned i = 0; i != geometry.facets.size(); ++i) {
    POLY_ASSERT(geometry.facets[i].size() == 2);
    qmesh.quantize(&PLCpoints[2*i], &boundary[i].x);
  }
  bp::IntPolygonSet boundarySet;
  boundarySet += boundary;

  // TODO -- handle holes.
  
  // Prepare new mesh data.
  std::vector<bp::IntPoint> newNodes;
  std::vector<std::pair<int, int> > newEdges;
  std::vector<std::vector<int> > newCellEdges;

  // Now walk each cell and clip it with the boundary.
  const unsigned ncells = qmesh.cellEdges.size();
  for (unsigned i = 0; i != ncells; ++i) {
    unsigned nverts = qmesh.cellEdges[i].size();

    // Copy the cell geometry to the quantized polygon type.
    bp::IntPolygon cell(nverts);
    for (unsigned j = 0; j != nverts; ++j) {
      int k = qmesh.cellEdges[i][j];
      if (k < 0) {
        k = ~k;
        cell[j] = qmesh.nodes[qmesh.edges[k].second];
      } else {
        cell[j] = qmesh.nodes[qmesh.edges[k].first];
      }
    }

    // Clip the cell against the boundary.
    bp::IntPolygonSet cellSet;
    cellSet += cell;
    cellSet &= boundarySet;
    std::cerr << "Final cell set: " << cellSet.size() << std::endl;

    // Read out the final cell geometry to the new QuantizedTessellation.
    nverts = cellSet[0].size();
    std::map<bp::IntPoint, int> node2id;
    std::map<std::pair<int, int>, int> edge2id;
    for (unsigned j = 0; j != nverts; ++j) {
      
      // Insert vertex 0.
      int old_size = node2id.size();
      const int j0 = internal::addKeyToMap(cellSet[0][j], node2id);
      if (j0 == old_size) {
        POLY_ASSERT(j0 == newNodes.size());
        newNodes.push_back(cellSet[0][j]);
      }

      // Insert vertex 1.
      old_size = node2id.size();
      const int j1 = internal::addKeyToMap(cellSet[0][(j+1)%nverts], node2id);
      if (j1 == old_size) {
        POLY_ASSERT(j1 == newNodes.size());
        newNodes.push_back(cellSet[0][(j+1)%nverts]);
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
