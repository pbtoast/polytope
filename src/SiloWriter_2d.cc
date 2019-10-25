#include "polytope.hh"
#ifdef HAVE_SILO
#include <fstream>
#include <set>
#include <cstring>
#include <sys/stat.h>
#include <dirent.h>
#include "silo.h"

#ifdef HAVE_MPI
// extern "C" {
#include "mpi.h"
#include "pmpio.h"
// }
#else
#define MPI_Comm int
#define MPI_COMM_WORLD 0
#endif

#include "SiloUtils.hh"

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
    POLY_ASSERT(nodes[n] >= 0);
    POLY_ASSERT(nodes[n] < mesh.nodes.size()/2);
  }
#endif
}
//-------------------------------------------------------------------
template <typename RealType>
void 
kullTraverseNodes(const Tessellation<2, RealType>& mesh,
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
}
//-------------------------------------------------------------------

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
     const map<string, RealType*>& nodeFields,
     const map<string, vector<int>*>& nodeTags,
     const map<string, RealType*>& edgeFields,
     const map<string, vector<int>*>& edgeTags,
     const map<string, RealType*>& faceFields,
     const map<string, vector<int>*>& faceTags,
     const map<string, RealType*>& cellFields,
     const map<string, vector<int>*>& cellTags,
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
  POLY_ASSERT(numFiles <= nproc);

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

  // Build the list of nodes describing the boundary faces.
  int numBoundaryFaces = mesh.boundaryFaces.size();
  vector<int> boundaryNodes(2*numBoundaryFaces);
  for (int i = 0; i < numBoundaryFaces; ++i)
  {
    boundaryNodes[2*i] = mesh.faces[mesh.boundaryFaces[i]][0];
    boundaryNodes[2*i+1] = mesh.faces[mesh.boundaryFaces[i]][1];
  }

  // Write the boundary face list.
  {
    vector<int> shapesize(size_t(1), 2), shapecnt(size_t(1), numBoundaryFaces);
    DBPutFacelist(file, (char*)"boundary_faces", numBoundaryFaces,
                  2, &boundaryNodes[0], boundaryNodes.size(), 0,
                  0, &shapesize[0], &shapecnt[0], 1, 0, 0, 0);
  }

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
  
  // DEBUG 
  std::cerr << " writing out UCDMesh in ::write " << std::endl;
  // DEBUG


  DBPutUcdmesh(file, (char*)"mesh", 2, coordnames, coords,
               numNodes, numCells, 
               (char *)"mesh_zonelist", NULL, DB_DOUBLE, optlist); 
               // (char *)"mesh_zonelist", (char*)"boundary_faces", DB_DOUBLE, optlist); 
  // DEBUG 
  std::cerr << " done writing out UCDMesh in ::write " << std::endl;
  // DEBUG

  DBPutZonelist2(file, (char*)"mesh_zonelist", numCells,
      2, &nodeList[0], nodeList.size(), 0, 0, 0,
      &shapetype[0], &shapesize[0], &shapecount[0],
      numCells, optlist);

  // DEBUG 
  std::cerr << " done writing out ZoneList2 in ::write " << std::endl;
  // DEBUG

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
  elemnames[0] = strDup("ncellfaces");
  elemlengths[0] = numCells;
  elemnames[2] = strDup("facecells");
  elemlengths[2] = conn.size() - 2*mesh.faces.size();
  elemnames[1] = strDup("cellfaces");
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
  elemnames[0] = strDup("nfacets");
  elemlengths[0] = 1;
  elemnames[1] = strDup("nfacetnodes");
  elemlengths[1] = mesh.convexHull.facets.size();
  elemnames[2] = strDup("facetnodes");
  elemlengths[2] = hull.size() - elemlengths[0] - elemlengths[1];
  DBPutCompoundarray(file, "convexhull", elemnames, elemlengths, 3, 
      (void*)&hull[0], hull.size(), DB_INT, 0);
  free(elemnames[0]);
  free(elemnames[1]);
  free(elemnames[2]);

  // Write out tag information.
  writeTagsToFile(nodeTags, file, DB_NODECENT);
  writeTagsToFile(edgeTags, file, DB_EDGECENT);
  writeTagsToFile(faceTags, file, DB_FACECENT);
  writeTagsToFile(cellTags, file, DB_ZONECENT);

  // Write out the field mesh data.
  // FIXME: We really should try to use the number of edges for edge fields.
  const int numFaces = mesh.faces.size();
  writeFieldsToFile<RealType>(nodeFields, file, numNodes, DB_NODECENT, optlist);
  writeFieldsToFile<RealType>(edgeFields, file, numFaces, DB_EDGECENT, optlist);
  writeFieldsToFile<RealType>(faceFields, file, numFaces, DB_FACECENT, optlist);
  writeFieldsToFile<RealType>(cellFields, file, numCells, DB_ZONECENT, optlist);

#if 0
  // Vector fields.
  {
    vector<RealType> xdata(mesh.numCells()), ydata(mesh.numCells());
    char* compNames[2];
    compNames[0] = new char[1024];
    compNames[1] = new char[1024];
    for (typename map<string, vector<2>*>::const_iterator iter = vectorFields.begin();
        iter != vectorFields.end(); ++iter)
    {
      snprintf(compNames[0], 1024, "%s_x" , iter->first.c_str());
      snprintf(compNames[1], 1024, "%s_y" , iter->first.c_str());
      vector<2>* data = iter->second;
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
    vector<vector<char*> > varNames(nodeFields.size() +
        edgeFields.size() +
        faceFields.size() +
        cellFields.size());
    vector<int> varTypes(numChunks, DB_UCDVAR);
    for (int i = 0; i < numChunks; ++i)
    {
      // Mesh.
      char meshName[1024];
      snprintf(meshName, 1024, "domain_%d/mesh", i);
      meshNames[i] = strDup(meshName);

      // Field data.
      int fieldIndex = 0;
      appendFieldNames<RealType>(nodeFields, fieldIndex, i, varNames);
      appendFieldNames<RealType>(edgeFields, fieldIndex, i, varNames);
      appendFieldNames<RealType>(faceFields, fieldIndex, i, varNames);
      appendFieldNames<RealType>(cellFields, fieldIndex, i, varNames);
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
    putMultivarInFile<RealType>(nodeFields, fieldIndex, varNames, varTypes, file, numChunks, optlist);
    putMultivarInFile<RealType>(edgeFields, fieldIndex, varNames, varTypes, file, numChunks, optlist);
    putMultivarInFile<RealType>(faceFields, fieldIndex, varNames, varTypes, file, numChunks, optlist);
    putMultivarInFile<RealType>(cellFields, fieldIndex, varNames, varTypes, file, numChunks, optlist);

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
      snprintf(masterFileName, 1024, "%s/%s-%d.silo", masterDirName.c_str(), prefix.c_str(), cycle);
    else
      snprintf(masterFileName, 1024, "%s/%s.silo", masterDirName.c_str(), prefix.c_str());
    int driver = DB_HDF5;
    // cerr << "Opening MASTER file " << masterFileName << endl;
    DBfile* file = DBCreate(masterFileName, DB_CLOBBER, DB_LOCAL, "Master file", driver);

    vector<char*> meshNames(numFiles*numChunks);
    vector<int> meshTypes(numFiles*numChunks, DB_UCDMESH);
    vector<vector<char*> > varNames(nodeFields.size() +
        edgeFields.size() +
        faceFields.size() +
        cellFields.size());
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
        meshNames[i*numChunks+c] = strDup(meshName);

        // Field data.
        int fieldIndex = 0;
        appendFieldNames<RealType>(nodeFields, fieldIndex, i, c, cycle, prefix, varNames);
        appendFieldNames<RealType>(edgeFields, fieldIndex, i, c, cycle, prefix, varNames);
        appendFieldNames<RealType>(faceFields, fieldIndex, i, c, cycle, prefix, varNames);
        appendFieldNames<RealType>(cellFields, fieldIndex, i, c, cycle, prefix, varNames);
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
    putMultivarInFile<RealType>(nodeFields, fieldIndex, varNames, varTypes, file, numFiles*numChunks, optlist);
    putMultivarInFile<RealType>(edgeFields, fieldIndex, varNames, varTypes, file, numFiles*numChunks, optlist);
    putMultivarInFile<RealType>(faceFields, fieldIndex, varNames, varTypes, file, numFiles*numChunks, optlist);
    putMultivarInFile<RealType>(cellFields, fieldIndex, varNames, varTypes, file, numFiles*numChunks, optlist);
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
template <typename RealType>
void 
SiloWriter<2, RealType>::writeKULLMesh(const Tessellation<2, RealType>& mesh, 
                                       const string& filePrefix,
                                       const string& directory,
                                       const map<string, RealType*>& nodeFields,
                                       const map<string, vector<int>*>& nodeTags,
                                       const map<string, RealType*>& edgeFields,
                                       const map<string, vector<int>*>& edgeTags,
                                       const map<string, RealType*>& faceFields,
                                       const map<string, vector<int>*>& faceTags,
                                       const map<string, RealType*>& cellFields,
                                       const map<string, vector<int>*>& cellTags,
                                       int cycle,
                                       RealType time,
                                       MPI_Comm comm,
                                       int numFiles,
                                       int mpiTag)
{
  // DEBUG 
  std::cerr << " SiloWriter_2d::writeKULLMesh -- top of routine  " << std::endl;
#ifdef HAVE_MPI  
  std::cerr << " SiloWriter_2d::writeKULLMesh -- HAVE_MPI  " << std::endl;
#else
  std::cerr << " SiloWriter_2d::writeKULLMesh -- we do not HAVE_MPI  " << std::endl;
#endif
  // DEBUG
  
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
  POLY_ASSERT(numFiles <= nproc);

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
    {
       mkdir((char*)masterDirName.c_str(), S_IRWXU | S_IRWXG);
    }
    else
    {
       closedir(masterDir);
    }
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
  {
     // snprintf(filename, 1024, "%s/%s-%d.silo", groupdirname, prefix.c_str(), cycle);
     snprintf(filename, 1024, "%s/%s.silo", groupdirname, prefix.c_str(), cycle);
  }
  else
  {
     snprintf(filename, 1024, "%s/%s.silo", groupdirname, prefix.c_str());
  }

  char dirname[1024];
  snprintf(dirname, 1024, "domain_%d", rankInGroup);
  DBfile* file = (DBfile*)PMPIO_WaitForBaton(baton, filename, dirname);
  
#else
  string dirname = directory;
  if (dirname.empty())
    dirname = ".";

  // create directory for SILO mesh
  DIR* masterDir = opendir(dirname.c_str());
  if (masterDir == 0)
  {
     mkdir((char*)dirname.c_str(), S_IRWXU | S_IRWXG);
  }

  if (cycle >= 0)
    snprintf(filename, 1024, "%s/%s-%d.silo", dirname.c_str(), prefix.c_str(), cycle);
  else
    snprintf(filename, 1024, "%s/%s.silo", dirname.c_str(), prefix.c_str());

  // DEBUG 
  std::cerr << " SiloWriter_2d::writeKULLMesh -- dirname = " << dirname << std::endl;
  std::cerr << " SiloWriter_2d::writeKULLMesh -- filename = " << filename << std::endl;
  std::cerr << " cycle = " << cycle << std::endl; 
  std::cerr << " prefix = " << prefix << std::endl; 
  // DEBUG
  
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  DBSetDir(file, "/");
#endif

  // DEBUG 
  std::cerr << " SiloWriter_2d::writeKULLMesh -- after #have_MPI " << std::endl;
  // DEBUG
  
  // Add cycle/time metadata if needed.
  DBoptlist* optlist = DBMakeOptlist(10);
  double dtime = static_cast<double>(time);
  if (cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
  if (dtime != -FLT_MAX)
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

  // add coordinate system, hard-wired to RZ for now ( FIXME )
  // DBAddOption(optlist, DBOPT_COORDSYS, DB_CYLINDRICAL);
  // add coordinate system, hard-wired to XY for now ( FIXME )
  int coordsystem = DB_CARTESIAN;
  DBAddOption(optlist, DBOPT_COORDSYS, (void *) &coordsystem);

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
  void const * const coords[2] = {&(x[0]), &(y[0])};
  // double* coords[2];
  // coords[0] = &(x[0]);
  // coords[1] = &(y[0]);

  // Build the list of nodes describing the boundary faces.
  int numBoundaryFaces = mesh.boundaryFaces.size();
  vector<int> boundaryNodes(2*numBoundaryFaces);
  for (int i = 0; i < numBoundaryFaces; ++i)
  {
    boundaryNodes[2*i] = mesh.faces[mesh.boundaryFaces[i]][0];
    boundaryNodes[2*i+1] = mesh.faces[mesh.boundaryFaces[i]][1];
  }

  // Write the boundary face list.
  {
    vector<int> shapesize(size_t(1), 2), shapecnt(size_t(1), numBoundaryFaces);
    DBPutFacelist(file, (char*)"boundary_faces", numBoundaryFaces,
                  2, &boundaryNodes[0], boundaryNodes.size(), 0,
                  0, &shapesize[0], &shapecnt[0], 1, 0, 0, 0);
  }

  // DEBUG 
  std::cerr << " SiloWriter_2d::writeKULLMesh -- before writing zones " << std::endl;
  // DEBUG

  // All zones are polygonal.
  size_t numShapes = 1; // only one shape --- polygonal
  int numCells = mesh.cells.size();
  vector<int> shapetype(numShapes, DB_ZONETYPE_POLYGON);
  vector<int> shapecount(numShapes, numCells);
  vector<int> nodeList;
  for (int i = 0; i < numCells; ++i)
  {
    // Gather the nodes from this cell in traversal order.
    vector<int> cellNodes;
    kullTraverseNodes(mesh, i, cellNodes);
    
    // cout << "cell " << i << ": " << cellNodes.size() << " - ";
    // for (int j = 0; j < cellNodes.size(); ++j) cout << cellNodes[j] << " ";
    // cout << endl;
    
    // Insert the cell's node connectivity into the node list.
    nodeList.push_back(cellNodes.size());
    nodeList.insert(nodeList.end(), cellNodes.begin(), cellNodes.end());
  }

  // From the SILO Manual, v4.10, July 2014, p2-103, "for a sequence of consecutive zones of
  // type DB_ZONETYPE_POLYHEDRON in a zonelist, the shapesize entry is taken to be the sum
  // of all the associated positions occupied in the nodelist data."
  std::cerr << " nodeList.size() = " << nodeList.size() << std::endl;

  vector<int> shapesize(numShapes, nodeList.size());
  
  std::cerr << " shapesize.size() = " << shapesize.size() << std::endl; 
  std::cerr << " shapesize[0] = " << shapesize[0] << std::endl; 
  std::cerr << " shapecount.size() = " << shapecount.size() << std::endl; 
  std::cerr << " shapecount[0] = " << shapecount[0] << std::endl; 
  std::cerr << " shapetype.size() = " << shapetype.size() << std::endl; 
  std::cerr << " shapetype[0] = " << shapetype[0] << std::endl; 


  std::cerr << " x.size() = " << x.size() << std::endl; 
  std::cerr << " y.size() = " << y.size() << std::endl; 
  std::cerr << " numNodes = " << numNodes << std::endl; 
  std::cerr << " numCells = " << numCells << std::endl; 

  // Write out the 2D polygonal mesh.
  DBPutUcdmesh(file,             // database file pointer   
               (char*)"MESH",    // name of the mesh in the SILO file
               2,                // number of spatial dimensions
               coordnames,       // parameter is ignored, can be set as NULL
               coords,           // Array of length ndims containing pointers to the coordinate arrays. 
               numNodes,         // Number of nodes in this UCD mesh.
               numCells,         // Number of zones in this UCD mesh.
               (char *)"mesh_zonelist",   // name of the zonelist structure (written with DBPutZonelist2). 
               // (char *)"boundary_faces",  // name of the facelist structure.  If none, set to NULL
               NULL,  // name of the facelist structure.  If none, set to NULL
               DB_DOUBLE,            // datatype of the coordinate arrays
               optlist);             // pointer to structure containing additional options

  std::cerr << " wrote out ucdMesh " << std::endl;

  // write out the zone-node-list
  DBPutZonelist2(file,             // database file pointer
                 (char*)"mesh_zonelist",  // name of the zonelist structure in SILO file
                 numCells,         // number of zones
                 2,                // number of spatial dimensions
                 &nodeList[0],     // Array containing node indices describing mesh zones
                 nodeList.size(),  // length of nodeList array
                 0,                // origin for indices in the nodeList array (zero or one)
                 0,                // low offset: number of ghost zones at beginning of nodeList
                 0,                // high offset: number of ghost zones at end of nodeList
                 &shapetype[0],    // array of length nshapes containing type of each zone
                 &shapesize[0],    // array of length nshapes containing number of nodes used by each shape
                 &shapecount[0],   // array of length nshapes containing number of zones having each shape
                 numShapes,        // number of shapes
                 optlist);

  std::cerr << " wrote out zonelist " << std::endl;

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
  elemnames[0] = strDup("ncellfaces");
  elemlengths[0] = numCells;
  elemnames[2] = strDup("facecells");
  elemlengths[2] = conn.size() - 2*mesh.faces.size();
  elemnames[1] = strDup("cellfaces");
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
  elemnames[0] = strDup("nfacets");
  elemlengths[0] = 1;
  elemnames[1] = strDup("nfacetnodes");
  elemlengths[1] = mesh.convexHull.facets.size();
  elemnames[2] = strDup("facetnodes");
  elemlengths[2] = hull.size() - elemlengths[0] - elemlengths[1];
  DBPutCompoundarray(file, "convexhull", elemnames, elemlengths, 3, 
      (void*)&hull[0], hull.size(), DB_INT, 0);
  free(elemnames[0]);
  free(elemnames[1]);
  free(elemnames[2]);

  // Write out tag information.
  writeTagsToFile(nodeTags, file, DB_NODECENT);
  writeTagsToFile(edgeTags, file, DB_EDGECENT);
  writeTagsToFile(faceTags, file, DB_FACECENT);
  writeTagsToFile(cellTags, file, DB_ZONECENT);

  // Write the material information
  // (This stupidly hard-wired while we are in a development phase.  --DSM 3/28/2019)
  size_t nMats = 1;
  int matnos[1] = {1};
  int ndims = 1;
  int dims[1];
  dims[0] = numCells;
  vector<int> matlist(numCells, 1);
  // If no mixed zones, these should be NULL, won't be used.
  int* mix_next = NULL;
  int* mix_mat = NULL;
  int* mix_zone = NULL;
  float* mix_vf = NULL;
  int mixlen = 0; // length of mixed zones vectors
  
  DBPutMaterial(file, "MATERIAL", "MESH",
                nMats,         // num materials
                matnos,        // array of material numbers
                &matlist[0],   // array with dimensions from 'dims' of materials in each zone
                dims, ndims, mix_next, mix_mat, mix_zone,
                mix_vf, mixlen, DB_DOUBLE, optlist);

  // Write out the field mesh data.
  // FIXME: We really should try to use the number of edges for edge fields.
  const int numFaces = mesh.faces.size();
  writeFieldsToFile<RealType>(nodeFields, file, numNodes, DB_NODECENT, optlist);
  writeFieldsToFile<RealType>(edgeFields, file, numFaces, DB_EDGECENT, optlist);
  writeFieldsToFile<RealType>(faceFields, file, numFaces, DB_FACECENT, optlist);
  writeFieldsToFile<RealType>(cellFields, file, numCells, DB_ZONECENT, optlist);

#if 0
  // Vector fields.
  {
    vector<RealType> xdata(mesh.numCells()), ydata(mesh.numCells());
    char* compNames[2];
    compNames[0] = new char[1024];
    compNames[1] = new char[1024];
    for (typename map<string, vector<2>*>::const_iterator iter = vectorFields.begin();
        iter != vectorFields.end(); ++iter)
    {
      snprintf(compNames[0], 1024, "%s_x" , iter->first.c_str());
      snprintf(compNames[1], 1024, "%s_y" , iter->first.c_str());
      vector<2>* data = iter->second;
      for (int i = 0; i < mesh.numCells(); ++i)
      {
        xdata[i] = data[i][0];
        ydata[i] = data[i][1];
      }
      void* vardata[2];
      vardata[0] = (void*)&xdata[0];
      vardata[1] = (void*)&ydata[0];
      DBPutUcdvar(file, (char*)iter->first.c_str(), (char*)"MESH",
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
    vector<vector<char*> > varNames(nodeFields.size() +
        edgeFields.size() +
        faceFields.size() +
        cellFields.size());
    vector<int> varTypes(numChunks, DB_UCDVAR);
    for (int i = 0; i < numChunks; ++i)
    {
      // Mesh.
      char meshName[1024];
      snprintf(meshName, 1024, "domain_%d/MESH", i);
      meshNames[i] = strDup(meshName);

      // Field data.
      int fieldIndex = 0;
      appendFieldNames<RealType>(nodeFields, fieldIndex, i, varNames);
      appendFieldNames<RealType>(edgeFields, fieldIndex, i, varNames);
      appendFieldNames<RealType>(faceFields, fieldIndex, i, varNames);
      appendFieldNames<RealType>(cellFields, fieldIndex, i, varNames);
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
    DBPutMultimesh(file, "MESH", numChunks, &meshNames[0], 
        &meshTypes[0], optlist);
    int fieldIndex = 0;
    putMultivarInFile<RealType>(nodeFields, fieldIndex, varNames, varTypes, file, numChunks, optlist);
    putMultivarInFile<RealType>(edgeFields, fieldIndex, varNames, varTypes, file, numChunks, optlist);
    putMultivarInFile<RealType>(faceFields, fieldIndex, varNames, varTypes, file, numChunks, optlist);
    putMultivarInFile<RealType>(cellFields, fieldIndex, varNames, varTypes, file, numChunks, optlist);

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
      snprintf(masterFileName, 1024, "%s/%s-%d.silo", masterDirName.c_str(), prefix.c_str(), cycle);
    else
      snprintf(masterFileName, 1024, "%s/%s.silo", masterDirName.c_str(), prefix.c_str());
    int driver = DB_HDF5;
    // cerr << "Opening MASTER file " << masterFileName << endl;
    DBfile* file = DBCreate(masterFileName, DB_CLOBBER, DB_LOCAL, "Master file", driver);

    vector<char*> meshNames(numFiles*numChunks);
    vector<int> meshTypes(numFiles*numChunks, DB_UCDMESH);
    vector<vector<char*> > varNames(nodeFields.size() +
        edgeFields.size() +
        faceFields.size() +
        cellFields.size());
    vector<int> varTypes(numFiles*numChunks, DB_UCDVAR);
    for (int i = 0; i < numFiles; ++i)
    {
      for (int c = 0; c < numChunks; ++c)
      {
        // Mesh.
        char meshName[1024];
        if (cycle >= 0)
          snprintf(meshName, 1024, "%d/%s-%d.silo:/domain_%d/MESH", i, prefix.c_str(), cycle, c);
        else
          snprintf(meshName, 1024, "%d/%s.silo:/domain_%d/MESH", i, prefix.c_str(), c);
        meshNames[i*numChunks+c] = strDup(meshName);

        // Field data.
        int fieldIndex = 0;
        appendFieldNames<RealType>(nodeFields, fieldIndex, i, c, cycle, prefix, varNames);
        appendFieldNames<RealType>(edgeFields, fieldIndex, i, c, cycle, prefix, varNames);
        appendFieldNames<RealType>(faceFields, fieldIndex, i, c, cycle, prefix, varNames);
        appendFieldNames<RealType>(cellFields, fieldIndex, i, c, cycle, prefix, varNames);
      }
    }

    DBoptlist* optlist = DBMakeOptlist(10);
    double dtime = static_cast<double>(time);
    if (cycle >= 0)
      DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    if (dtime != -FLT_MAX)
      DBAddOption(optlist, DBOPT_DTIME, &dtime);

    // Write the multimesh and variable data, and close the file.
    DBPutMultimesh(file, "MESH", numFiles*numChunks, &meshNames[0], 
        &meshTypes[0], optlist);
    int fieldIndex = 0;
    putMultivarInFile<RealType>(nodeFields, fieldIndex, varNames, varTypes, file, numFiles*numChunks, optlist);
    putMultivarInFile<RealType>(edgeFields, fieldIndex, varNames, varTypes, file, numFiles*numChunks, optlist);
    putMultivarInFile<RealType>(faceFields, fieldIndex, varNames, varTypes, file, numFiles*numChunks, optlist);
    putMultivarInFile<RealType>(cellFields, fieldIndex, varNames, varTypes, file, numFiles*numChunks, optlist);
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
  char masterFileName[1024];
  char topName[11];
  snprintf(topName, 11, "OvlTop");
  
  std::cerr << "Opening MASTER file " << directory.c_str() << "/" << topName << ".silo" << std::endl;
  snprintf(masterFileName, 1024, "%s/%s.silo", directory.c_str(), topName);
  DBfile* masterfile = DBCreate(masterFileName, DB_CLOBBER, DB_LOCAL, "Master file", driver);
  std::cerr << "MASTER file " << directory.c_str() << "/" << topName << ".silo is opened" << std::endl;

  // allocate memory
  numFiles = 1;
  vector<char*> meshNames(numFiles);
  vector<char*> matNames(numFiles);
  vector<int> meshTypes(numFiles, DB_UCDMESH);
  DBoptlist* masterOptlist = DBMakeOptlist(10);

  std::cerr << "allocated memory for MASTER file " << directory.c_str() << "/" << topName << ".silo" << std::endl;
  
  for (int i = 0; i < numFiles; ++i)
  {
     // Mesh.
     char meshName[1024];
     snprintf(meshName, 1024, "%s/%s.silo:MESH", directory.c_str(), prefix.c_str());
     meshNames[i] = strDup(meshName);
  }
  std::cerr << "wrote names for MASTER file " << directory.c_str() << "/" << topName << ".silo" << std::endl;
  
  DBPutMultimesh(masterfile, "MMESH", numFiles, &meshNames[0], &meshTypes[0], masterOptlist);
  std::cerr << "putMultiMesh names into MASTER file " << directory.c_str() << "/" << topName << ".silo" << std::endl;

  // hard-coded for development, SORRY, FIXME, should not be hard-wiring materials in 
  char matName[10];
  snprintf(matName, 10, "%s", "air");
  matNames[0] = strDup(matName);
  size_t numMats = 1;
  DBPutMultimat(masterfile, "MMATERIAL", numMats, &matNames[0], masterOptlist);

  // free memory
  DBFreeOptlist(masterOptlist);
  for (int i = 0; i < numFiles; ++i)
     free(meshNames[i]);
  std::cerr << "freed memory" << std::endl;
  
  DBClose(masterfile);
  DBClose(file);
#endif
  // DEBUG 
  std::cerr << " SiloWriter_2d::writeKULLMesh -- end of routine ################## " << std::endl << std::endl;
  // DEBUG
  
}
//-------------------------------------------------------------------

// Explicit instantiation.
template class SiloWriter<2, double>;

} // end namespace

#endif

