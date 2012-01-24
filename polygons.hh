#ifndef POLYTOPE_POLYGONS_HH
#define POLYTOPE_POLYGONS_HH

#include <vector>
#include "Mesh.hh"

namespace polytope
{

namespace polygons
{

//! Traverse the nodes of the given polygonal cell (within the given mesh), 
//! appending the indices of the nodes traversed to \a nodes.
void traverseNodes(const Mesh<2>& mesh,
                   const Cell<2>& cell, 
                   std::vector<int>& nodes);

//! Traverse the convex hull of the polygon with the given points in the 
//! plane (in counter-clockwise order), storing the traversal order in 
//! \a indices.
void traverseConvexHull(const std::vector<Point<2> >& points,
                        std::vector<int>& indices);


}

}
#endif
