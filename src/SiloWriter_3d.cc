#include "SiloWriter.hh"
#include <fstream>
#include <set>
#include "polygons.hh"

#ifdef HAVE_SILO
#include "silo.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "pmpio.h"
#endif

#endif

namespace Charybdis
{

#ifdef HAVE_SILO

using namespace std;

#ifdef HAVE_MPI
namespace 
{

//-------------------------------------------------------------------
PMPIO_createFile(const char* filename,
                 const char* dirname,
                 void* userData)
{
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, prefix.c_str(), 
                          driver);
  DBMkDir(file, dirname);
  DBSetDir(file, dirname);
  return (void*)file;
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void* 
PMPIO_openFile(const char* filename, 
               const char* dirname,
               PMPIO_iomode_t iomode, 
               void* userData)
{
  int driver = DB_HDF5;
  DBfile* file;
  if (iomode == PMPIO_WRITE)
  { 
    file = DBCreate(filename, 0, DB_LOCAL, prefix.c_str(), 
                    driver);
    DBMkDir(file, dirname);
    DBSetDir(file, dirname);
  }
  else
  {
    file = DBOpen(filename, driver, DB_READ);
    DBSetDir(file, dirname);
  }
  return (void*)file;
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void*
PMPIO_closeFile(void* file,
                void* userData)
{
  DBClose((DBfile*)file);
}
//-------------------------------------------------------------------

}
#endif

//-------------------------------------------------------------------
template <typename Real>
void 
SiloWriter<3, Real>::
write(const Mesh<3>& mesh, 
      const map<string, Real*>& scalarFields,
      const map<string, Vector<3>*>& vectorFields,
      const string& filePrefix,
      int cycle,
      Real time,
      MPI_Comm comm,
      int numFiles,
      int mpiTag)
{
  int nproc;
  MPI_Comm_size(comm, &nproc);
  ASSERT(numFiles <= nproc);

  // Strip .silo off of the prefix if it's there.
  string prefix = filePrefix;
  size_t index = prefix.find(".silo");
  if (index >= 0)
    prefix.erase(index);

  char filename[1024];
  if (cycle >= 0)
    snprintf(filename, 1024, "%s-%d.silo", prefix.c_str(), cycle);
  else
    snprintf(filename, 1024, "%s.silo", prefix.c_str());

  // Open a file in Silo/HDF5 format for writing.
#ifdef HAVE_MPI
  PMPIO_baton_t* baton = PMPIO_Init(numFiles, PMPIO_WRITE, comm, mpiTag, 
                                    &PMPIO_createFile, 
                                    &PMPIO_openFile, 
                                    &PMPIO_closeFile,
                                    0);
  DBfile* file = (DBfile*)PMPIO_WaitForBaton(baton, filename, "/");
#else
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, prefix.c_str(), 
                          driver);
  DBSetDir(file, "/");
#endif

  // Add cycle/time metadata if needed.
  DBoptlist* optlist = DBMakeOptlist(10);
  double dtime = static_cast<double>(time);
  if (cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
  if (dtime != -FLT_MAX)
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

  // This is optional for now, but we'll give it anyway.
  char *coordnames[3];
  coordnames[0] = (char*)"xcoords";
  coordnames[1] = (char*)"ycoords";
  coordnames[2] = (char*)"zcoords";

  // Node coordinates.
  int N = mesh.numNodes();
  vector<double> x(N), y(N), z(N);
  for (int i = 0; i < N; ++i)
  {
    const Point<3>& node = mesh.node(i);
    x[i] = node[0];
    y[i] = node[1];
    z[i] = node[2];
  }
  double* coords[3];
  coords[0] = &(x[0]);
  coords[1] = &(y[0]);
  coords[2] = &(z[0]);

  // Figure out face-node connectivity. We do this by computing centers
  // for all the cells and then using them to define face normals, 
  // from which node orderings can be determining using a convex hull 
  // determination algorithm (gift wrapping).
  vector<Point<3> > cellCenters(mesh.numCells());
  for (int c = 0; c < mesh.numCells(); ++c)
  {
    const Cell<3>& cell = *mesh.cell(c);
    int numNodes = 0;
    for (int f = 0; f < cell.numFaces; ++f)
    {
      const Face<3>& face = *mesh.face(cell.faces[f]);
      for (int n = 0; n < face.numNodes; ++n)
      {
        const Point<3>& node = mesh.node(face.nodes[n]);
        cellCenters[c][0] += node[0];
        cellCenters[c][1] += node[1];
        cellCenters[c][2] += node[2];
        ++numNodes;
      }
    }
    cellCenters[c][0] /= numNodes;
    cellCenters[c][1] /= numNodes;
    cellCenters[c][2] /= numNodes;
  }
  vector<int> faceNodeCounts(mesh.numFaces()), 
              allFaceNodes;
  for (int f = 0; f < mesh.numFaces(); ++f)
  {
    const Face<3>& face = *mesh.face(f);

    // Compute the normal vector for the face, pointing outward from 
    // its first cell.
    ASSERT(face.numNodes >= 3);
    Point<3> faceCenter;
    for (int n = 0; n < face.numNodes; ++n)
    {
      const Point<3>& node = mesh.node(face.nodes[n]);
      faceCenter[0] += node[0];
      faceCenter[1] += node[1];
      faceCenter[2] += node[2];
    }
    faceCenter[0] /= face.numNodes;
    faceCenter[1] /= face.numNodes;
    faceCenter[2] /= face.numNodes;
    Vector<3> v1(faceCenter, mesh.node(face.nodes[0]));
    Vector<3> v2, normal;
    for (int n = 1; n < face.numNodes; ++n)
    {
      v2 = Vector<3>(faceCenter, mesh.node(face.nodes[n]));
      normal = cross(v1, v2);
      if (normal.mag() > 1e-14) break;
    }
    normal /= normal.mag();
    const Point<3>& xcell = cellCenters[face.cell1];
    Vector<3> v3(xcell, faceCenter);
    if (normal.dot(v3) < 0.0)
      normal *= -1.0;

    // Now project the coordinates of the face's nodes to the plane
    // with the given normal and centered about the face center.
    vector<Point<2> > points(face.numNodes);
    Vector<3> e1 = v1.unit(), // Basis vectors in the plane.
              e2 = cross(normal, e1);
    for (int p = 0; p < points.size(); ++p)
    {
      // Compute the perpendicular component of the point
      // with location v:
      // (vPerp = v - (n o v)n and so forth).
      Vector<3> v(xcell, mesh.node(face.nodes[p]));
      Vector<3> vPerp = v - normal.dot(v)*normal;

      // Project it to the plane.
      points[p][0] = vPerp.dot(e1);
      points[p][1] = vPerp.dot(e2);
    }

    // Find the node order by traversing the convex hull of 
    // the points within the plane, appending them to allFaceNodes.
    vector<int> indices;
    polygons::traverseConvexHull(points, indices);
    faceNodeCounts[f] = indices.size();
    for (int n = 0; n < indices.size(); ++n)
      allFaceNodes.push_back(face.nodes[indices[n]]);
  }

  // Figure out cell-face connectivity.
  vector<int> cellFaceCounts(mesh.numCells()), 
              allCellFaces;
  for (int i = 0; i < mesh.numCells(); ++i)
  {
    const Cell<3>& cell = *mesh.cell(i);
    cellFaceCounts[i] = cell.numFaces;
    allCellFaces.insert(allCellFaces.end(), 
                        cell.faces, cell.faces + cell.numFaces);
  }

  // The polyhedral zone list is referred to in the options list.
  DBAddOption(optlist, DBOPT_PHZONELIST, (char*)"mesh_zonelist");

  // Write out the 3D polyhedral mesh.
  DBPutUcdmesh(file, (char*)"mesh", 3, coordnames, coords,
               mesh.numNodes(), mesh.numCells(), 0, 0,
               DB_DOUBLE, optlist); 

  // Write the connectivity information.
  DBPutPHZonelist(file, (char*)"mesh_zonelist", 
                  faceNodeCounts.size(), &faceNodeCounts[0], 
                  allFaceNodes.size(), &allFaceNodes[0], 0, 
                  cellFaceCounts.size(), &cellFaceCounts[0],
                  allCellFaces.size(), &allCellFaces[0], 
                  0, 0, mesh.numCells()-1, optlist);

  // Write out the cell-centered mesh data.

  // Scalar fields.
  for (map<string, Real*>::const_iterator iter = scalarFields.begin();
       iter != scalarFields.end(); ++iter)
  {
    DBPutUcdvar1(file, (char*)iter->first.c_str(), (char*)"mesh",
                 (void*)iter->second, mesh.numCells(), 0, 0,
                 DB_DOUBLE, DB_ZONECENT, optlist);
  }

  // Vector fields.
  {
    vector<Real> xdata(mesh.numCells()), ydata(mesh.numCells()),
                 zdata(mesh.numCells());
    char* compNames[3];
    compNames[0] = new char[1024];
    compNames[1] = new char[1024];
    compNames[2] = new char[1024];
    for (map<string, Vector<3>*>::const_iterator iter = vectorFields.begin();
         iter != vectorFields.end(); ++iter)
    {
      snprintf(compNames[0], 1024, "%s_x" , iter->first.c_str());
      snprintf(compNames[1], 1024, "%s_y" , iter->first.c_str());
      snprintf(compNames[2], 1024, "%s_z" , iter->first.c_str());
      Vector<3>* data = iter->second;
      for (int i = 0; i < mesh.numCells(); ++i)
      {
        xdata[i] = data[i][0];
        ydata[i] = data[i][1];
        zdata[i] = data[i][2];
      }
      void* vardata[3];
      vardata[0] = (void*)&xdata[0];
      vardata[1] = (void*)&ydata[0];
      vardata[2] = (void*)&zdata[0];
      DBPutUcdvar(file, (char*)iter->first.c_str(), (char*)"mesh",
                  3, compNames, vardata, mesh.numCells(), 0, 0,
                  DB_DOUBLE, DB_ZONECENT, optlist);
    }
    delete [] compNames[0];
    delete [] compNames[1];
    delete [] compNames[2];
  }

  // Clean up.
  DBFreeOptlist(optlist);

#ifdef HAVE_MPI
  PMPIO_HandOffBaton(baton, (void*)file);
  PMPIO_finish(baton);

  // Write the multi-block objects to the file.
  // FIXME
#else
  // Write the file.
  DBClose(file);
#endif
}
//-------------------------------------------------------------------

} // end namespace

