#ifdef HAVE_SILO
#include "polytope.hh"
#include <fstream>
#include <set>
#include <cstring>
#include <sys/stat.h>
#include <dirent.h>
#include "silo.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "pmpio.h"
#endif

namespace polytope
{

using namespace std;

namespace 
{

//-------------------------------------------------------------------
// Traverse the nodes of cell i within the given tessellation in 
// order, writing their indices to nodes.  We rely here on two
// assumptions:
// 1.  cellFaces are given such that the faces are in counter-clockwise
//     order around the cell. 
// 2.  if cellFaces[j] > 0, the nodes of the face cellFaces[j] are 
//     given in counter-clockwise orientation for cell i,
//     otherwise the nodes of ~cellFaces[j] (the 1s complement) 
//     are in *clockwise* order and need to be reversed.
//-------------------------------------------------------------------
template <typename RealType>
void 
traverseNodes(const Tessellation<2, RealType>& mesh,
              int i,
              vector<int>& nodes)
{
  const vector<int>& cellFaces = mesh.cells[i];
  for (int j = 0; j != cellFaces.size(); ++j) 
  {
    int k = cellFaces[j];
    nodes.push_back(k >= 0 ? mesh.faces[ k][0] :
                             mesh.faces[~k][1]);
  }
  nodes.push_back(nodes.front());

#ifndef NDEBUG
  // Make sure we don't have any garbage in our list of nodes.
  for (int n = 0; n < nodes.size(); ++n)
  {
    ASSERT(nodes[n] >= 0);
    ASSERT(nodes[n] < mesh.nodes.size()/2);
  }
#endif
}
//-------------------------------------------------------------------

#ifdef HAVE_MPI

//-------------------------------------------------------------------
void*
PMPIO_createFile(const char* filename,
                 const char* dirname,
                 void* userData)
{
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
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
    file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
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
void
PMPIO_closeFile(void* file,
                void* userData)
{
  DBClose((DBfile*)file);
}
//-------------------------------------------------------------------

#endif
}

//-------------------------------------------------------------------
template <typename RealType>
void 
SiloWriter<2, RealType>::
write(const Tessellation<2, RealType>& mesh, 
      const map<string, RealType*>& fields,
      const string& filePrefix,
      const string& directory,
      int cycle,
      RealType time,
      MPI_Comm comm,
      int numFiles,
      int mpiTag)
{
  // Strip .silo off of the prefix if it's there.
  string prefix = filePrefix;
  int index = prefix.find(".silo");
  if (index >= 0)
    prefix.erase(index);

  // Open a file in Silo/HDF5 format for writing.
  char filename[1024];
#ifdef HAVE_MPI
  int nproc = 1, rank = 0;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  if (numFiles == -1)
    numFiles = nproc;
  ASSERT(numFiles <= nproc);

  // We put the entire data set into a directory named after the 
  // prefix, and every process gets its own subdirectory therein.

  // Create the master directory if we need to.
  string masterDirName = directory;
  if (masterDirName.empty())
  {
    char dirname[1024];
    snprintf(dirname, 1024, "%s-%d", filePrefix.c_str(), nproc);
    masterDirName = dirname;
  }
  if (rank == 0)
  {
    DIR* masterDir = opendir(directory.c_str());
    if (masterDir == 0)
      mkdir((char*)masterDirName.c_str(), S_IRWXU | S_IRWXG);
    else
      closedir(masterDir);
    MPI_Barrier(comm);
  }
  else
  {
    MPI_Barrier(comm);
  }

  // Initialize poor man's I/O and figure out group ranks.
  PMPIO_baton_t* baton = PMPIO_Init(numFiles, PMPIO_WRITE, comm, mpiTag, 
                                    &PMPIO_createFile, 
                                    &PMPIO_openFile, 
                                    &PMPIO_closeFile,
                                    0);
  int groupRank = PMPIO_GroupRank(baton, rank);
  int rankInGroup = PMPIO_RankInGroup(baton, rank);

  // Create a subdirectory for each group.
  char groupdirname[1024];
  snprintf(groupdirname, 1024, "%s/%d", masterDirName.c_str(), groupRank);
  if (rankInGroup == 0)
  {
    DIR* groupDir = opendir(groupdirname);
    if (groupDir == 0)
      mkdir((char*)groupdirname, S_IRWXU | S_IRWXG);
    else
      closedir(groupDir);
    MPI_Barrier(comm);
  }
  else
  {
    MPI_Barrier(comm);
  }

  // Determine a file name.
  if (cycle >= 0)
    snprintf(filename, 1024, "%s/%s-%d.silo", groupdirname, prefix.c_str(), cycle);
  else
    snprintf(filename, 1024, "%s/%s.silo", groupdirname, prefix.c_str());

  char dirname[1024];
  snprintf(dirname, 1024, "domain_%d", rankInGroup);
  DBfile* file = (DBfile*)PMPIO_WaitForBaton(baton, filename, dirname);
#else
  string dirname = directory;
  if (dirname.empty())
    dirname = ".";

  if (cycle >= 0)
    snprintf(filename, 1024, "%s/%s-%d.silo", dirname.c_str(), prefix.c_str(), cycle);
  else
    snprintf(filename, 1024, "%s/%s.silo", dirname.c_str(), prefix.c_str());

  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
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
  int numNodes = mesh.nodes.size() / 2;
  vector<double> x(numNodes), y(numNodes);
  for (int i = 0; i < numNodes; ++i)
  {
    x[i] = mesh.nodes[2*i];
    y[i] = mesh.nodes[2*i+1];
  }
  double* coords[2];
  coords[0] = &(x[0]);
  coords[1] = &(y[0]);

  // All zones are polygonal.
  int numCells = mesh.cells.size();
  vector<int> shapesize(numCells, 0),
              shapetype(numCells, DB_ZONETYPE_POLYGON),
              shapecount(numCells, 1),
              nodeList;
  for (int i = 0; i < numCells; ++i)
  {
    // Gather the nodes from this cell in traversal order.
    vector<int> cellNodes;
    traverseNodes(mesh, i, cellNodes);
// cout << "cell " << i << ": ";
// for (int j = 0; j < cellNodes.size(); ++j)
// cout << cellNodes[j] << " ";
// cout << endl;
    // Insert the cell's node connectivity into the node list.
    nodeList.push_back(cellNodes.size());
    nodeList.insert(nodeList.end(), cellNodes.begin(), cellNodes.end());
  }

  // Write out the 2D polygonal mesh.
  DBPutUcdmesh(file, (char*)"mesh", 2, coordnames, coords,
               numNodes, numCells, 
               (char *)"mesh_zonelist", 0, DB_DOUBLE, optlist); 
  DBPutZonelist2(file, (char*)"mesh_zonelist", numCells,
                 2, &nodeList[0], nodeList.size(), 0, 0, 0,
                 &shapetype[0], &shapesize[0], &shapecount[0],
                 numCells, optlist);

  // Write out the cell-face connectivity data.
  vector<int> conn(numCells);
  int elemlengths[3];
  char* elemnames[3];
  for (int c = 0; c < numCells; ++c)
    conn[c] = mesh.cells[c].size();
  for (int c = 0; c < numCells; ++c)
  {
    for (int f = 0; f < mesh.cells[c].size(); ++f)
      conn.push_back(mesh.cells[c][f]);
  }
  for (int f = 0; f < mesh.faceCells.size(); ++f)
  {
    conn.push_back(mesh.faceCells[f][0]);
    conn.push_back(mesh.faceCells[f][0]);
  }
  elemnames[0] = strdup("ncellfaces");
  elemlengths[0] = numCells;
  elemnames[2] = strdup("facecells");
  elemlengths[2] = conn.size() - 2*mesh.faces.size();
  elemnames[1] = strdup("cellfaces");
  elemlengths[1] = conn.size() - elemlengths[2] - elemlengths[0];
  DBPutCompoundarray(file, "conn", elemnames, elemlengths, 3, 
                     (void*)&conn[0], conn.size(), DB_INT, 0);
  free(elemnames[0]);
  free(elemnames[1]);
  free(elemnames[2]);

  // Write out convex hull data.
  vector<int> hull(1+mesh.convexHull.facets.size());
  hull[0] = mesh.convexHull.facets.size();
  for (int f = 0; f < mesh.convexHull.facets.size(); ++f)
    hull[1+f] = mesh.convexHull.facets[f].size();
  for (int f = 0; f < mesh.convexHull.facets.size(); ++f)
    for (int n = 0; n < mesh.convexHull.facets[f].size(); ++n)
      hull.push_back(mesh.convexHull.facets[f][n]);
  elemnames[0] = strdup("nfacets");
  elemlengths[0] = 1;
  elemnames[1] = strdup("nfacetnodes");
  elemlengths[1] = mesh.convexHull.facets.size();
  elemnames[2] = strdup("facetnodes");
  elemlengths[2] = hull.size() - elemlengths[0] - elemlengths[1];
  DBPutCompoundarray(file, "convexhull", elemnames, elemlengths, 3, 
                     (void*)&hull[0], hull.size(), DB_INT, 0);
  free(elemnames[0]);
  free(elemnames[1]);
  free(elemnames[2]);

  // Write out the cell-centered mesh data.

  // Scalar fields.
  for (typename map<string, RealType*>::const_iterator iter = fields.begin();
       iter != fields.end(); ++iter)
  {
    DBPutUcdvar1(file, (char*)iter->first.c_str(), (char*)"mesh",
                 (void*)iter->second, numCells, 0, 0,
                 DB_DOUBLE, DB_ZONECENT, optlist);
  }

#if 0
  // Vector fields.
  {
    vector<RealType> xdata(mesh.numCells()), ydata(mesh.numCells());
    char* compNames[2];
    compNames[0] = new char[1024];
    compNames[1] = new char[1024];
    for (typename map<string, Vector<2>*>::const_iterator iter = vectorFields.begin();
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
#endif

  // Clean up.
  DBFreeOptlist(optlist);

#ifdef HAVE_MPI
  // Write the multi-block objects to the file if needed.
  int numChunks = nproc / numFiles;
  if (rankInGroup == 0)
  {
    vector<char*> meshNames(numChunks);
    vector<int> meshTypes(numChunks, DB_UCDMESH);
    vector<vector<char*> > varNames(fields.size());
    vector<int> varTypes(numChunks, DB_UCDVAR);
    for (int i = 0; i < numChunks; ++i)
    {
      // Mesh.
      char meshName[1024];
      snprintf(meshName, 1024, "domain_%d/mesh", i);
      meshNames[i] = strdup(meshName);

      // Field data.
      int fieldIndex = 0;
      for (typename map<string, RealType*>::const_iterator iter = fields.begin();
           iter != fields.end(); ++iter, ++fieldIndex)
      {
        char varName[1024];
        snprintf(varName, 1024, "domain_%d/%s", i, iter->first.c_str());
        varNames[fieldIndex].push_back(strdup(varName));
      }
    }

    // Stick cycle and time in there if needed.
    DBoptlist* optlist = DBMakeOptlist(10);
    double dtime = static_cast<double>(time);
    if (cycle >= 0)
      DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    if (dtime != -FLT_MAX)
      DBAddOption(optlist, DBOPT_DTIME, &dtime);

    // Write the mesh and variable data.
    DBSetDir(file, "/");
    DBPutMultimesh(file, "mesh", numChunks, &meshNames[0], 
                   &meshTypes[0], optlist);
    int fieldIndex = 0;
    for (typename map<string, RealType*>::const_iterator iter = fields.begin();
         iter != fields.end(); ++iter, ++fieldIndex)
    {
      DBPutMultivar(file, iter->first.c_str(), numChunks, 
                    &varNames[fieldIndex][0], &varTypes[0], optlist);
    }

    // Clean up.
    DBFreeOptlist(optlist);
    for (int i = 0; i < numChunks; ++i)
      free(meshNames[i]);
    for (int f = 0; f < varNames.size(); ++f)
      for (int i = 0; i < numChunks; ++i)
        free(varNames[f][i]);
  }

  // Write the file.
  PMPIO_HandOffBaton(baton, (void*)file);
  PMPIO_Finish(baton);

  // Finally, write the uber-master file.
  if (rank == 0)
  {
    char masterFileName[1024];
    if (cycle >= 0)
      snprintf(masterFileName, 1024, "%s-%d/%s-%d.silo", prefix.c_str(), nproc, prefix.c_str(), cycle);
    else
      snprintf(masterFileName, 1024, "%s-%d/%s.silo", prefix.c_str(), nproc, prefix.c_str());
    int driver = DB_HDF5;
    DBfile* file = DBCreate(masterFileName, DB_CLOBBER, DB_LOCAL, "Master file", driver);

    vector<char*> meshNames(numFiles*numChunks);
    vector<int> meshTypes(numFiles*numChunks, DB_UCDMESH);
    vector<vector<char*> > varNames(fields.size());
    vector<int> varTypes(numFiles*numChunks, DB_UCDVAR);
    for (int i = 0; i < numFiles; ++i)
    {
      for (int c = 0; c < numChunks; ++c)
      {
        // Mesh.
        char meshName[1024];
        if (cycle >= 0)
          snprintf(meshName, 1024, "%d/%s-%d.silo:/domain_%d/mesh", i, prefix.c_str(), cycle, c);
        else
          snprintf(meshName, 1024, "%d/%s.silo:/domain_%d/mesh", i, prefix.c_str(), c);
        meshNames[i*numChunks+c] = strdup(meshName);

        // Field data.
        int fieldIndex = 0;
        for (typename map<string, RealType*>::const_iterator iter = fields.begin();
             iter != fields.end(); ++iter, ++fieldIndex)
        {
          char varName[1024];
          if (cycle >= 0)
            snprintf(varName, 1024, "%d/%s-%d.silo:/domain_%d/%s", i, prefix.c_str(), cycle, c, iter->first.c_str());
          else
            snprintf(varName, 1024, "%d/%s.silo:/domain_%d/%s", i, prefix.c_str(), c, iter->first.c_str());
          varNames[fieldIndex].push_back(strdup(varName));
        }
      }
    }

    DBoptlist* optlist = DBMakeOptlist(10);
    double dtime = static_cast<double>(time);
    if (cycle >= 0)
      DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    if (dtime != -FLT_MAX)
      DBAddOption(optlist, DBOPT_DTIME, &dtime);

    // Write the multimesh and variable data, and close the file.
    DBPutMultimesh(file, "mesh", numFiles*numChunks, &meshNames[0], 
                   &meshTypes[0], optlist);
    int fieldIndex = 0;
    for (typename map<string, RealType*>::const_iterator iter = fields.begin();
         iter != fields.end(); ++iter, ++fieldIndex)
    {
      DBPutMultivar(file, iter->first.c_str(), numFiles*numChunks, 
                    &(varNames[fieldIndex][0]), &varTypes[0], optlist);
    }
    DBClose(file);

    // Clean up.
    DBFreeOptlist(optlist);
    for (int i = 0; i < numFiles*numChunks; ++i)
      free(meshNames[i]);
    for (int f = 0; f < varNames.size(); ++f)
      for (int i = 0; i < numFiles*numChunks; ++i)
        free(varNames[f][i]);
  }
#else
  // Write the file.
  DBClose(file);
#endif
}
//-------------------------------------------------------------------

// Explicit instantiation.
template class SiloWriter<2, double>;

} // end namespace

#endif

