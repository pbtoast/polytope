#include "SiloWriter.hh"
#include "polygons.hh"
#include <fstream>
#include <set>

#ifdef HAVE_SILO
#include "silo.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "pmpio.h"
#endif

#endif

namespace polytope
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
SiloWriter<2, Real>::
write(const Mesh<2>& mesh, 
      const map<string, Real*>& scalarFields,
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
  char *coordnames[2];
  coordnames[0] = (char*)"xcoords";
  coordnames[1] = (char*)"ycoords";

  // Node coordinates.
  int N = mesh.numNodes();
  vector<double> x(N), y(N);
  for (int i = 0; i < N; ++i)
  {
    const Point<2>& node = mesh.node(i);
    x[i] = node[0];
    y[i] = node[1];
  }
  double* coords[2];
  coords[0] = &(x[0]);
  coords[1] = &(y[0]);

  // All zones are polygonal.
  vector<int> shapesize(mesh.numCells(), 0),
              shapetype(mesh.numCells(), DB_ZONETYPE_POLYGON),
              shapecount(mesh.numCells(), 1),
              nodeList;
  for (int i = 0; i < mesh.numCells(); ++i)
  {
    const Cell<2>& cell = *mesh.cell(i);
    vector<int> cellNodes;

    // Gather the nodes from this cell in traversal order.
    polygons::traverseNodes(mesh, cell, cellNodes);

    // Insert the cell's node connectivity into the node list.
    nodeList.push_back(cellNodes.size());
    nodeList.insert(nodeList.end(), cellNodes.begin(), cellNodes.end());
  }

  // Write out the 2D polygonal mesh.
  DBPutUcdmesh(file, (char*)"mesh", 2, coordnames, coords,
               mesh.numNodes(), mesh.numCells(), 
               (char *)"mesh_zonelist", 0, DB_DOUBLE, optlist); 
  DBPutZonelist2(file, (char*)"mesh_zonelist", mesh.numCells(),
                 2, &nodeList[0], nodeList.size(), 0, 0, 0,
                 &shapetype[0], &shapesize[0], &shapecount[0],
                 mesh.numCells(), optlist);

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
    vector<Real> xdata(mesh.numCells()), ydata(mesh.numCells());
    char* compNames[2];
    compNames[0] = new char[1024];
    compNames[1] = new char[1024];
    for (map<string, Vector<2>*>::const_iterator iter = vectorFields.begin();
         iter != vectorFields.end(); ++iter)
    {
      snprintf(compNames[0], 1024, "%s_x" , iter->first.c_str());
      snprintf(compNames[1], 1024, "%s_y" , iter->first.c_str());
      Vector<2>* data = iter->second;
      for (int i = 0; i < mesh.numCells(); ++i)
      {
        xdata[i] = data[i][0];
        ydata[i] = data[i][1];
      }
      void* vardata[2];
      vardata[0] = (void*)&xdata[0];
      vardata[1] = (void*)&ydata[0];
      DBPutUcdvar(file, (char*)iter->first.c_str(), (char*)"mesh",
                  2, compNames, vardata, mesh.numCells(), 0, 0,
                  DB_DOUBLE, DB_ZONECENT, optlist);
    }
    delete [] compNames[0];
    delete [] compNames[1];
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

