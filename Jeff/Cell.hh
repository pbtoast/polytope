#ifndef CELL_HH
#define CELL_HH

#include "MeshDiff.hh"
#include "charybdisconfig.hh"
#include <vector>

namespace Charybdis
{

template <int Dimension> class Mesh;
template <int Dimension> class Point;

//! \class Cell
//! This class defines a Voronoi cell.
template <int D>
class Cell
{
  public:

  //! Default constructor - creates an unassociated Voronoi cell.
  Cell(): faces(0), numFaces(0), centroid(), volume() {}

  //! Destructor.
  ~Cell() 
  {
    if (faces != 0)
      delete [] faces;
  }

  //! Assign faces to the cell. Note that data is 
  //! copied from \a faces.
  void assignFaces(const int* faces, int numFaces)
  {
    ASSERT(faces != 0);
    ASSERT(numFaces > 0);
    if (this->faces != 0)
      delete [] this->faces;
    this->faces = new int[numFaces];
    this->numFaces = numFaces;
    std::copy(faces, faces + numFaces, this->faces);
  }

  //! Assign faces to the cell. Note that data is 
  //! copied from \a faces.
  void assignFaces(const std::vector<int>& faces)
  {
    assignFaces(&faces[0], faces.size());
  }

  friend std::ostream& operator<<(std::ostream& s, const Cell& cell)
  {
    s << "Cell(faces = ";
    for (int f = 0; f < cell.numFaces-1; ++f)
      s << cell.faces[f] << ", ";
    s << cell.faces[cell.numFaces-1] << ")";
    return s;
  }

  //! Array of face indices.
  int* faces;

  //! Number of faces.
  int numFaces;

  //! Centroid
  Point<D> centroid;

  //! Volume
  Real volume;

  private:

  template <int DD>
  friend void tessellate(const std::vector<Point<DD> >& points,
                         const std::vector<Point<DD> >& ghostPoints,
                         Mesh<DD>& mesh,
                         MeshDiff& meshDiff);

  // This object may not be copied.
  Cell(const Cell&);
  Cell& operator=(const Cell&);

};

} // end namespace

#endif
