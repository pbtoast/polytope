#ifndef TESSELLATE_HH
#define TESSELLATE_HH

#include <vector>
#include "Point.hh"
#include "Mesh.hh"
#include "MeshDiff.hh"

namespace Charybdis
{

//! Computes the Voronoi tessellation of the given generator points (including ghost generators), 
//! storing the tessellation in the given mesh and storing the changes from the initial state of 
//! the mesh to \a diff.
template <int Dimension>
void tessellate(const std::vector<Point<Dimension> >& points,
                const std::vector<Point<Dimension> >& ghostPoints,
                Mesh<Dimension>& mesh,
                MeshDiff& diff)
{
  // No general recipe.
}

// Template specializations.
template <>
void tessellate(const std::vector<Point<2> >& points,
                const std::vector<Point<2> >& ghostPoints,
                Mesh<2>& mesh,
                MeshDiff& diff);

template <>
void tessellate(const std::vector<Point<3> >& points,
                const std::vector<Point<3> >& ghostPoints,
                Mesh<3>& mesh,
                MeshDiff& diff);

//! Computes the Voronoi tessellation of the given generator points, 
//! storing the tessellation in the given mesh and storing the changes 
//! from the initial state of the mesh to \a diff.
template <int Dimension>
void tessellate(const std::vector<Point<Dimension> >& points,
                Mesh<Dimension>& mesh,
                MeshDiff& diff)
{
  std::vector<Point<Dimension> > ghosts;
  tessellate(points, ghosts, mesh, diff);
}

//! Computes the Voronoi tessellation of the given generator points (including ghost generators), 
//! storing the tessellation in the given mesh. This version does not keep track of changes and 
//! is appropriate for initial tessellations.
template <int Dimension>
void tessellate(const std::vector<Point<Dimension> >& points,
                const std::vector<Point<Dimension> >& ghostPoints,
                Mesh<Dimension>& mesh)
{
  MeshDiff diff;
  tessellate(points, ghostPoints, mesh, diff);
}

//! Computes the Voronoi tessellation of the given generator points,
//! storing the tessellation in the given mesh. This version does not 
//! keep track of changes and is appropriate for initial tessellations.
template <int Dimension>
void tessellate(const std::vector<Point<Dimension> >& points,
                Mesh<Dimension>& mesh)
{
  MeshDiff diff;
  std::vector<Point<Dimension> > ghosts;
  tessellate(points, ghosts, mesh, diff);
}

} // end namespace

#endif
