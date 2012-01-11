#include "tessellate.hh"
#include "voro_2d/voro++_2d.hh"
#include <cmath>
#include <float.h>
#include <vector>
#include <set>
#include <map>

namespace Charybdis
{

using namespace voro;
using std::min;
using std::max;
using std::vector;
using std::set;
using std::pair;
using std::swap;
using std::find;
using std::binary_search;

//-------------------------------------------------------------------
// 2D specialization
template <>
void 
tessellate(const std::vector<Point<2> >& points,
           const std::vector<Point<2> >& ghostPoints,
           Mesh<2>& mesh,
           MeshDiff& diff)
{
  diff.clear();

  // Determine a bounding box for the points.
  Real xmin = FLT_MAX, xmax = -FLT_MAX,
       ymin = FLT_MAX, ymax = -FLT_MAX;
  for (int i = 0; i < points.size(); ++i)
  {
    const Point<2>& p = points[i];
    xmin = min(xmin, p[0]), xmax = max(xmax, p[0]);
    ymin = min(ymin, p[1]), ymax = max(ymax, p[1]);
  }
  for (int i = 0; i < ghostPoints.size(); ++i)
  {
    const Point<2>& p = ghostPoints[i];
    xmin = min(xmin, p[0]), xmax = max(xmax, p[0]);
    ymin = min(ymin, p[1]), ymax = max(ymax, p[1]);
  }
  // Add a little buffering.
  {
    Real Lx = xmax - xmin,
         Ly = ymax - ymin;
    xmin -= 0.05*Lx; xmax += 0.05*Lx;
    ymin -= 0.05*Ly; ymax += 0.05*Ly;
  }

  // Make a list of the neighbors in the old mesh.
  int oldNumCells = mesh.numCells() + mesh.numGhostCells();
  vector<vector<int> > oldNeighbors(oldNumCells);
  for (int i = 0; i < oldNumCells; ++i)
  {
    // Make a list of the neighboring cells of the old cell. They should be sorted 
    // in ascending order of neighboring cell indices.
    Cell<2>* cell = mesh.m_cells[i];
    for (int f = 0; f < cell->numFaces; ++f)
      oldNeighbors[i].push_back(mesh.m_faces[cell->faces[f]]->otherCell(i));
  }

  // Make sure the new Mesh has the right number of cells.
  int newNumCells = points.size() + ghostPoints.size();
  if (oldNumCells < newNumCells)
  {
    diff.addCells(oldNumCells, newNumCells - oldNumCells);
    mesh.m_cells.resize(newNumCells);
    for (int i = oldNumCells; i < newNumCells; ++i)
      mesh.m_cells[i] = new Cell<2>();
    if (mesh.numGhostCells() < ghostPoints.size())
    {
      diff.deleteGhostCells(ghostPoints.size() - mesh.numGhostCells());
      mesh.m_numGhostCells = ghostPoints.size();
    }
    else if (mesh.numGhostCells() > ghostPoints.size())
    {
      diff.deleteGhostCells(mesh.numGhostCells() - ghostPoints.size());
      mesh.m_numGhostCells = ghostPoints.size();
    }

    // Old neighbor entries have empty sets for new cells.
    oldNeighbors.resize(newNumCells);
  }
  else if (oldNumCells > newNumCells)
  {
    diff.deleteCells(newNumCells, oldNumCells - newNumCells);
    for (int i = newNumCells; i < (oldNumCells - newNumCells); ++i)
      delete mesh.m_cells[i];
    mesh.m_cells.resize(newNumCells);
    if (mesh.numGhostCells() < ghostPoints.size())
    {
      diff.deleteGhostCells(ghostPoints.size() - mesh.numGhostCells());
      mesh.m_numGhostCells = ghostPoints.size();
    }
    else if (mesh.numGhostCells() > ghostPoints.size())
    {
      diff.deleteGhostCells(mesh.numGhostCells() - ghostPoints.size());
      mesh.m_numGhostCells = ghostPoints.size();
    }
  }
  else 
  {
    if (mesh.numGhostCells() < ghostPoints.size())
    {
      diff.deleteGhostCells(ghostPoints.size() - mesh.numGhostCells());
      mesh.m_numGhostCells = ghostPoints.size();
    }
    else if (mesh.numGhostCells() > ghostPoints.size())
    {
      diff.deleteGhostCells(mesh.numGhostCells() - ghostPoints.size());
      mesh.m_numGhostCells = ghostPoints.size();
    }
  }

  // Nodes of the mesh are completely re-indexed.
  mesh.m_nodes.clear(); 
  mesh.m_numGhostNodes = 0;

  // Now construct the new tessellation using Voro++.
  vector<vector<int> > newNeighbors(newNumCells); // A table mapping cells to their neighbors.
  vector<vector<int> > cellNodes(newNumCells); // Cell -> node mapping.
  {
    // Number of cells in the grid used by Voro for neighbor search.
    int nx = 10, ny = 10;

    // Number of points per grid cell.
    int ptsPerCell = 16;

    // Initialize the container.
    container_2d con(xmin,xmax,ymin,ymax,nx,ny,false,false,ptsPerCell);

    // Add the points to the container.
    for (int i = 0; i < points.size(); ++i)
      con.put(i, points[i][0], points[i][1]);
    for (int i = 0; i < ghostPoints.size(); ++i)
      con.put(points.size() + i, ghostPoints[i][0], ghostPoints[i][1]);
    ASSERT(con.total_particles() == points.size() + ghostPoints.size());

    voronoicell_neighbor_2d cell; // Use cells with neighbor tracking.
    c_loop_all_2d loop(con); // Loop over all cells.
    if (loop.start())
    {
      do 
      {
        // Compute the cell.
        if (con.compute_cell(cell, loop))
        {
          int i = loop.pid(); // Fetch the cell index.

          // Fetch the cell center coordinates.
          // These must be added to the node coordinates.
          double *pp = con.p[loop.ij] + con.ps*loop.q; 
          Real xc = pp[0], yc = pp[1];

          // Compute the centroid and volume of the cell.
          cell.centroid(mesh.m_cells[i]->centroid[0], mesh.m_cells[i]->centroid[1]);
          mesh.m_cells[i]->centroid[0] += xc; 
          mesh.m_cells[i]->centroid[1] += yc;
          mesh.m_cells[i]->volume = cell.area();

          vector<int>& neighbors = newNeighbors[i];

          // Record the neighbors of the cell and sort them.
          cell.neighbors(neighbors);
          sort(neighbors.begin(), neighbors.end());

          // FIXME: Voro seems to be sticking negative numbers in here
          // FIXME: in certain cases. we need to get rid of these.
          int k = 0;
          while (neighbors[k] < 0) 
            ++k;
          if (k > 0)
            neighbors.erase(neighbors.begin(), neighbors.begin() + k);

          // Record the node positions of the cell.
          for (int n = 0; n < cell.p; ++n)
          {
            // Compute the node position.
            Point<2> node(xc + 0.5*cell.pts[2*n], yc + 0.5*cell.pts[2*n+1]);

            // Check to see if any of the neighboring cells have a node in
            // the same position as this one.
            bool nodeAlreadyExists = false;
            for (int nn = 0; nn < neighbors.size(); ++nn)
            {
              int neighbor = neighbors[nn];
              for (int p = 0; p < cellNodes[neighbor].size(); ++p)
              {
                int otherNode = cellNodes[neighbor][p];
                if (distance(mesh.m_nodes[otherNode], node) < 1e-14)
                {
                  // This is the same node.
                  nodeAlreadyExists = true;
                  cellNodes[i].push_back(otherNode);
                }
              }
            }
            if (!nodeAlreadyExists)
            {
              // New node.
              cellNodes[i].push_back(mesh.m_nodes.size());
              mesh.m_nodes.push_back(node);
            }
          }
        }
      }
      while (loop.inc());
    }

    // Find nodes that are attached only to ghost cells. These are the 
    // ghost nodes, and they must be moved to the back of the nodes array.
    set<int> ghostNodes, interiorNodes;
    for (int gcell = mesh.numCells(); gcell < mesh.m_cells.size(); ++gcell)
    {
      for (int n1 = 0; n1 < cellNodes[gcell].size(); ++n1)
      {
        int node1 = cellNodes[gcell][n1];

        // Have we already determined that this is a ghost node or 
        // an interior node?
        if (ghostNodes.find(node1) != ghostNodes.end()) continue;
        if (interiorNodes.find(node1) != interiorNodes.end()) continue;

        // Go over the neighboring cells and see whether this node 
        // touches an interior cell. If it doesn't, it's a ghost node.
        for (int nc = 0; nc < newNeighbors[gcell].size(); ++nc)
        {
          int icell = newNeighbors[gcell][nc];
          if (icell >= mesh.numCells()) continue; // Consider interior cells only!
          for (int n2 = 0; n2 < cellNodes[icell].size(); ++n2)
          {
            if (node1 == cellNodes[icell][n2])
            {
              // The node touches an interior cell. It's not a ghost node.
              interiorNodes.insert(node1);
              break;
            }
          }
        }

        // If it's not marked as an interior node at this point, it's 
        // a ghost node.
        if (interiorNodes.find(node1) == interiorNodes.end())
          ghostNodes.insert(node1);
      }
    }
    mesh.m_numGhostNodes = ghostNodes.size();

    // Now that we have the ghost nodes, reorder the nodes in the mesh.
    map<int, int> swaps;
    for (int n = mesh.numNodes(); n < mesh.numNodes() + mesh.m_numGhostNodes; ++n)
    {
      // Is this a ghost node? If not, we must swap it with the minimum ghost node.
      if (!ghostNodes.empty() and (ghostNodes.find(n) == ghostNodes.end()))
      {
        int minGhostNode = *ghostNodes.begin();
        swaps[n] = minGhostNode;
        swap(mesh.m_nodes[n], mesh.m_nodes[minGhostNode]);
        ghostNodes.erase(ghostNodes.begin());
      }
    }

    // Reorder the nodes in the cellNodes array.
    for (int c = 0; c < cellNodes.size(); ++c)
    {
      for (int n = 0; n < cellNodes[c].size(); ++n)
      {
        map<int, int>::const_iterator iter = swaps.find(cellNodes[c][n]);
        if (iter != swaps.end())
          cellNodes[c][n] = iter->second;
      }
    }
  }

  // Now loop over the cells of the mesh and edit them as necessary, 
  // creating faces and such.
  set<pair<int, int> > facesCreated; // Memo for created faces.
  vector<vector<int> > cellFaces(newNumCells); // Faces attached to cells.
  set<int> ghostFaces;
  for (int i = 0; i < newNumCells; ++i)
  {
    Cell<2>& cell = *mesh.m_cells[i];
    if (cell.numFaces > 0)
    {
      cellFaces[i].insert(cellFaces[i].begin(), 
                          cell.faces, cell.faces + cell.numFaces);
    }

    int n = 0;
    for (int newn = 0; newn < newNeighbors[i].size(); ++newn)
    {
//      if (newn >= oldNeighbors[i].size()) break;

      // Index of the cell opposite this face.
      int newj = newNeighbors[i][newn];

      // If the faces of the cell i have remained the same, 
      // should be identical. Otherwise, we will insert/delete neighbors from 
      // the list of faces as necessary.
      if ((newn >= oldNeighbors[i].size()) or (newj < oldNeighbors[i][n])) // Perhaps a new face must be inserted?
      {
        // Has this face already been created?
        pair<int, int> cellPair(min(i, newj), max(i, newj));
        if (facesCreated.find(cellPair) != facesCreated.end())
        {
          ASSERT(!cellFaces[i].empty());
          continue;
        }

        // We only insert a face for (i, j) for which there exists an entry (j, i).
        if (!binary_search(newNeighbors[newj].begin(), newNeighbors[newj].end(), i))
          continue;

        int j = newNeighbors[i][newn]; // Other cell index.
        int faceIndex = mesh.m_faces.size();

        // Add a new face and attach it to its cells.
        diff.addFace(mesh.m_faces.size());
        Face<2>* newFace = new Face<2>();
        newFace->assignCells(i, newNeighbors[i][newn]);
        mesh.m_faces.push_back(newFace);

        // Is this face a ghost face? It is if both its cells are ghost cells.
        if ((cellPair.first >= mesh.numCells()) and (cellPair.second >= mesh.numCells()))
          ghostFaces.insert(faceIndex);

        // Associate the face with each of its cells, inserting it in
        // the proper sorted order.
        int index1 = std::distance(cellFaces[i].begin(), 
                                   lower_bound(cellFaces[i].begin(),
                                               cellFaces[i].end(),
                                               newj));
        cellFaces[i].insert(cellFaces[i].begin() + index1, faceIndex);
        int index2 = std::distance(cellFaces[j].begin(), 
                                   lower_bound(cellFaces[j].begin(),
                                               cellFaces[j].end(),
                                               i));
        cellFaces[j].insert(cellFaces[j].begin() + index2, faceIndex);
        facesCreated.insert(cellPair);

        // Now figure out which nodes are needed for this face. These are the nodes
        // attached to both cells.
        vector<int> faceNodes;
        for (int n1 = 0; n1 < cellNodes[i].size(); ++n1)
        {
          int node1 = cellNodes[i][n1];
          for (int n2 = 0; n2 < cellNodes[newj].size(); ++n2)
          {
            int node2 = cellNodes[newj][n2];
            if (node1 == node2) // Same node!
            {
              if (faceNodes.empty() or (faceNodes[0] != node1))
              {
                faceNodes.push_back(node1);
                if (faceNodes.size() == 2) break; // in 2D, no face has more than 2 nodes.
              }
            }
          }
          if (faceNodes.size() == 2) break;
        }
        newFace->assignNodes(faceNodes);
        ASSERT(newFace->numNodes == 2);
        ASSERT(newFace->nodes[0] != newFace->nodes[1]);
      }
      else
      {
        ASSERT(!oldNeighbors[i].empty());
        int j = oldNeighbors[i][n];
        if (j < newj) // A previous face no longer exists.
        {
          int faceIndex = cellFaces[i][j];

          // Swap this face with the last one.
          diff.swapFaces(cellFaces[i][j], mesh.m_faces.size()-1);
          swap(mesh.m_faces[cellFaces[i][j]], mesh.m_faces.back());

          // Dissociate the face from each of its cells.
          swap(cellFaces[i][j], cellFaces[i].back());
          cellFaces[i].pop_back();
          swap(*find(cellFaces[j].begin(), cellFaces[j].end(), faceIndex),
              cellFaces[j].back());
          cellFaces[j].pop_back();

          // Delete the last face.
          diff.deleteFace(mesh.m_faces.size());
          delete mesh.m_faces.back();
          mesh.m_faces.pop_back();
        }
        else // This face existed before and still exists.
        {
          // We still have to figure out which nodes this face has.
          Face<2>* face = mesh.m_faces[cellFaces[i][j]];
          vector<int> faceNodes;
          for (int n1 = 0; n1 < cellNodes[i].size(); ++n1)
          {
            int node1 = cellNodes[i][n1];
            for (int n2 = 0; n2 < cellNodes[newj].size(); ++n2)
            {
              int node2 = cellNodes[newj][n2];
              if (node1 == node2) // Same node!
              {
                if (faceNodes.empty() or (faceNodes[0] != node1))
                {
                  faceNodes.push_back(node1);
                  if (faceNodes.size() == 2) break; // in 2D, no face has more than 2 nodes.
                }
              }
            }
            if (faceNodes.size() == 2) break;
          }
          face->assignNodes(faceNodes);

          // Next neighbor of cell i.
          ++n;
        }
      }
    }

    // We should have at least 1 face.
    ASSERT(cellFaces[i].size() >= 1);

    // At this point, the number of faces should not 
    // exceed the number of neighbors.
    ASSERT(cellFaces[i].size() <= newNeighbors[i].size());
  }

  // Reorder the faces in the mesh and in the cellFaces array before we add 
  // them to the cells.
  if (!ghostFaces.empty())
  {
    map<int, int> swaps;
    for (int f = mesh.numFaces(); f < mesh.numFaces() + mesh.m_numGhostFaces; ++f)
    {
      // Is this a ghost face? If not, we must swap it with the minimum ghost face.
      if (!ghostFaces.empty() and (ghostFaces.find(f) == ghostFaces.end()))
      {
        int minGhostFace = *ghostFaces.begin();
        swaps[f] = minGhostFace;
        swap(mesh.m_faces[f], mesh.m_faces[minGhostFace]);
        ghostFaces.erase(ghostFaces.begin());
      }
    }

    // Reorder the faces in the cellFaces array.
    for (int c = 0; c < cellFaces.size(); ++c)
    {
      for (int f = 0; f < cellFaces[c].size(); ++f)
      {
        map<int, int>::const_iterator iter = swaps.find(cellFaces[c][f]);
        if (iter != swaps.end())
          cellFaces[c][f] = iter->second;
      }
    }
  }

  // Finally, add the faces to the cells.
  for (int c = 0; c < mesh.numCells() + mesh.m_numGhostCells; ++c)
    mesh.m_cells[c]->assignFaces(cellFaces[c]);
}
//-------------------------------------------------------------------

} // end namespace

#ifdef BUILD_TESTS
// Test executable.
#include "testMacros.hh"
#include "repartition.hh"
#include "writeVTKPlot.hh"
#include "writeSiloPlot.hh"
using namespace Charybdis;
using std::vector;
using std::map;

//-------------------------------------------------------------------
void 
testLattice2D(int numProcs)
{
  // Create a uniform lattice of points.
  vector<Point<2> > points;

  int Nx = 10, Ny = 10;
  Real L = 1.0, dx = L/Nx, dy = L/Ny;
  for (int i = 0; i < Nx; ++i)
  {
    for (int j = 0; j < Ny; ++j)
    {
      Point<2> p((i+0.5)*dx, (j+0.5)*dy);
      points.push_back(p);
    }
  }

  // Load-balance the points.
  repartition(MPI_COMM_WORLD, points);

  // Generate the tessellation.
  Mesh<2> mesh;
  MeshDiff diff;
  tessellate(points, mesh, diff);

  // Firstly, no node should appear in any face more than once.
  for (int i = 0; i < mesh.numCells(); ++i)
  {
    const Cell<2>* cell = mesh.cell(i);
    for (int j = 0; j < mesh.cell(i)->numFaces; ++j)
    {
      const Face<2>* face = mesh.face(cell->faces[j]);
      for (int k = 0; k < face->numNodes; ++k)
      {
        TEST_ASSERT(count(face->nodes, 
                          face->nodes + face->numNodes,
                          face->nodes[0]) == 1);
      }
    }
  }

  // Make sure none of the cells has more than 4 faces 
  // and that each face has 2 nodes. Also: count the 
  // number of cells with 2, 3, and 4 faces.
  TEST_ASSERT(mesh.numCells() == Nx*Ny);
  int twoFaced = 0, threeFaced = 0, fourFaced = 0;
  for (int i = 0; i < mesh.numCells(); ++i)
  {
    const Cell<2>* cell = mesh.cell(i);
    TEST_ASSERT(cell->numFaces <= 4);
    for (int j = 0; j < mesh.cell(i)->numFaces; ++j)
    {
      const Face<2>* face = mesh.face(cell->faces[j]);
      TEST_ASSERT(face->numNodes == 2);
    }
    if (cell->numFaces == 2)
      ++twoFaced;
    else if (cell->numFaces == 3)
      ++threeFaced;
    else
      ++fourFaced;
  }

  // There should be exactly 4 cells with 2 faces.
  TEST_ASSERT(twoFaced == 4);

  // There should be exactly 16 cells with 3 faces.
  TEST_ASSERT(threeFaced == 32);

  // There should be exactly 80 cells with 4 faces.
  TEST_ASSERT(fourFaced == 64);

  map<string, Real*> nothing;
  map<string, Vector<2>*> vnothing;
  writeVTKPlot(mesh, nothing, vnothing, "lattice-2d.vtk");
#ifdef HAVE_SILO
  writeSiloPlot(mesh, nothing, vnothing, "lattice-2d.silo");
#endif
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
int 
main(int argc, char** argv)
{
  TEST_LIBRARIES(geometry, utils, voro_2d);

  MPI_Init(&argc, &argv);

  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  ASSERT((numProcs == 1) or ((numProcs % 4) == 0));

  testLattice2D(numProcs);

  MPI_Finalize();
  return 0;
}
//-------------------------------------------------------------------
#endif


