#ifndef MESH_HH
#define MESH_HH

#include <vector>
#include "Point.hh"
#include "Cell.hh"
#include "Face.hh"
#include "MeshDiff.hh"

namespace Charybdis
{

//! \class Mesh
//! This class defines a Voronoi mesh consisting of a set of polyhedral or polygonal cells determined 
//! a Voronoi graph consisting of a set of generators.
template <int Dimension>
class Mesh
{
  public:

  //! Default constructor - creates an empty Voronoi mesh.
  Mesh(): m_cells(), m_faces(), m_nodes(), m_numGhostCells(0), 
          m_numGhostFaces(0), m_numGhostNodes(0), m_firstGlobalIndex(0) {}

  //! Destructor.
  ~Mesh() 
  {
    // Clean up!
    for (int i = 0; i < m_cells.size(); ++i)
      delete m_cells[i];
    for (int i = 0; i < m_faces.size(); ++i)
      delete m_faces[i];
  }

  //! Returns the number of cells in the mesh, excluding ghost cells.
  int numCells() const 
  { 
    return m_cells.size() - m_numGhostCells;
  }

  //! Returns the number of ghost cells in the mesh.
  int numGhostCells() const { return m_numGhostCells; }

  //! Returns a pointer to the cell with the given index.
  const Cell<Dimension>* const cell(int index) const { return m_cells[index]; }

  //! Returns the number of faces in the mesh, excluding ghost faces.
  int numFaces() const 
  { 
    return m_faces.size() - m_numGhostFaces;
  }

  //! Returns the number of ghost faces in the mesh.
  int numGhostFaces() const { return m_numGhostFaces; }

  //! Returns a pointer to the face with the given index.
  const Face<Dimension>* const face(int index) const { return m_faces[index]; }

  //! Returns the number of nodes in the mesh, excluding ghost nodes.
  int numNodes() const { return m_nodes.size() - m_numGhostNodes; }

  //! Returns the number of ghost nodes in the mesh.
  int numGhostNodes() const { return m_numGhostNodes; }

  //! Returns the position of the given node.
  const Point<Dimension>& node(int index) const { return m_nodes[index]; }

  private:

  template <int D>
  friend void tessellate(const std::vector<Point<D> >& points,
                         const std::vector<Point<D> >& ghostPoints,
                         Mesh<D>& mesh,
                         MeshDiff& meshDiff);

  //---------------------------------------------------------------
  // Data members -- accessible by friends.
  //---------------------------------------------------------------

  //! Cells.
  std::vector<Cell<Dimension>*> m_cells;

  //! Faces.
  std::vector<Face<Dimension>*> m_faces;

  //! Node positions.
  std::vector<Point<Dimension> > m_nodes;

  //! Number of ghost cells, faces, nodes.
  int m_numGhostCells, m_numGhostFaces, m_numGhostNodes;

  //! (Global) index of first locally-maintained cell.
  unsigned long m_firstGlobalIndex; 

  // This object may not be copied.
  Mesh(const Mesh&);
  Mesh& operator=(const Mesh&);

};

// Useful type definitions.
typedef Mesh<2> PolygonalMesh;
typedef Mesh<3> PolyhedralMesh;

} // end namespace

#endif
