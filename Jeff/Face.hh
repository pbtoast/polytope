#ifndef FACE_HH
#define FACE_HH

#include "MeshDiff.hh"
#include <set>

namespace Charybdis
{

//! \class Face
//! This class defines a Voronoi face.
template <int D>
class Face
{
  public:

  //! Default constructor - creates an unassociated Voronoi cell.
  Face(): cell1(-1), cell2(-1), nodes(0), numNodes(0), centroid(), area(0.0) {}

  //! Destructor.
  ~Face() 
  {
    if (nodes != 0)
      delete [] nodes;
  }

  //! Assign cells to the face. The face will store the cells 
  //! in ascending index order.
  void assignCells(int cell1, int cell2)
  {
    ASSERT(cell1 >= 0);
    ASSERT(cell2 >= 0);
    ASSERT(cell1 != cell2);
    this->cell1 = std::min(cell1, cell2);
    this->cell2 = std::max(cell1, cell2);
  }

  //! Assign nodes to the face. Note that data is 
  //! copied from \a nodes.
  void assignNodes(const int* nodes, int numNodes)
  {
    ASSERT(nodes != 0);
    ASSERT(numNodes > 0);
    if (this->nodes != 0)
      delete [] this->nodes;
    this->nodes = new int[numNodes];
    this->numNodes = numNodes;
    std::copy(nodes, nodes + numNodes, this->nodes);
  }

  //! Assign nodes to the face. Note that data is 
  //! copied from \a nodes.
  void assignNodes(const std::vector<int>& nodes)
  {
    assignNodes(&nodes[0], nodes.size());
  }

  //! Assign nodes to the face. Note that data is 
  //! copied from \a nodes.
  void assignNodes(const std::set<int>& nodes)
  {
    ASSERT(!nodes.empty());
    if (this->nodes != 0)
      delete [] this->nodes;
    this->nodes = new int[numNodes];
    this->numNodes = nodes.size();
    std::copy(nodes.begin(), nodes.end(), this->nodes);
  }

  //! Index of first attached cell.
  int cell1;

  //! Index of second attached cell.
  int cell2;

  //! Array of nodes attached to the faces.
  int* nodes;

  //! Number of nodes attached to the face.
  int numNodes;

  //! Face centroid.
  Point<D> centroid;

  //! Area.
  Real area;

  //! Given the index of one of the cells attached to this face, returns the 
  //! other one.
  int otherCell(int cellIndex) const 
  {
    ASSERT((cellIndex == cell1) or (cellIndex == cell2));
    if (cellIndex == cell1)
      return cell2;
    else return cell1;
  }

  friend std::ostream& operator<<(std::ostream& s, const Face& face)
  {
    s << "Face(nodes = ";
    for (int n = 0; n < face.numNodes-1; ++n)
      s << face.nodes[n] << ", ";
    s << face.nodes[face.numNodes-1] << ")";
    return s;
  }

  private:

  template <int DD>
  friend void tessellate(const std::vector<Point<DD> >& points,
                         const std::vector<Point<DD> >& ghostPoints,
                         Mesh<DD>& mesh,
                         MeshDiff& meshDiff);

  // This object may not be copied.
  Face(const Face&);
  Face& operator=(const Face&);

};

} // end namespace

#endif
