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
SiloReader<3, RealType>::
read(Tessellation<3, RealType>& mesh, 
     map<string, RealType*>& fields,
     const string& filePrefix,
     const string& directory,
     int cycle,
     RealType& time,
     MPI_Comm comm,
     int numFiles,
     int mpiTag)
{
  // Strip .silo off of the prefix if it's there.
  string prefix = filePrefix;
  int index = prefix.find(".silo");
  if (index >= 0)
    prefix.erase(index);

  // Open a file in Silo/HDF5 format for reading.
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

  // Check for the existence of the master directory.
  string masterDirName = directory;
  if (masterDirName.empty())
  {
    char dirname[1024];
    snprintf(dirname, 1024, "%s-%d", filePrefix.c_str(), nproc);
    masterDirName = dirname;
  }
  DIR* masterDir = opendir(directory.c_str());
  if (masterDir == 0)
  {
    char err[1024];
    snprintf(err, 1024, "Could not find the directory %s", masterDirName.c_str());
    error(err);
  }

  // Initialize poor man's I/O and figure out group ranks.
  PMPIO_baton_t* baton = PMPIO_Init(numFiles, PMPIO_READ, comm, mpiTag, 
                                    &PMPIO_createFile, 
                                    &PMPIO_openFile, 
                                    &PMPIO_closeFile,
                                    0);
  int groupRank = PMPIO_GroupRank(baton, rank);
  int rankInGroup = PMPIO_RankInGroup(baton, rank);

  // Figure out the subdirectory for this group.
  char groupdirname[1024];
  snprintf(groupdirname, 1024, "%s/%d", masterDirName.c_str(), groupRank);
  DIR* groupDir = opendir(groupdirname);
  if (groupDir == 0)
  {
    char err[1024];
    snprintf(err, 1024, "Could not find the directory %s", groupdirname);
    error(err);
  }

  // Determine the file name.
  if (cycle >= 0)
    snprintf(filename, 1024, "%s/%s-%d.silo", groupdirname, prefix.c_str(), cycle);
  else
    snprintf(filename, 1024, "%s/%s.silo", groupdirname, prefix.c_str());

  char dirname[1024];
  snprintf(dirname, 1024, "domain_%d", rankInGroup);
  DBfile* file = (DBfile*)PMPIO_WaitForBaton(baton, filename, dirname);
  DBSetDir(file, dirname);
#else
  string dirname = directory;
  if (dirname.empty())
    dirname = ".";

  if (cycle >= 0)
    snprintf(filename, 1024, "%s/%s-%d.silo", dirname.c_str(), prefix.c_str(), cycle);
  else
    snprintf(filename, 1024, "%s/%s.silo", dirname.c_str(), prefix.c_str());

  int driver = DB_HDF5;
  DBfile* file = DBOpen(filename, driver, DB_READ);
  DBSetDir(file, "/");
#endif

  // Retrieve the mesh. Note that we must deallocate the storage 
  // for this object after we're through!
  DBucdmesh* dbmesh = DBGetUcdmesh(file, "mesh");

  // Extract time.
  time = dbmesh->dtime;

  // Node coordinates.
  mesh.nodes.resize(3*dbmesh->nnodes);
  for (int i = 0; i < dbmesh->nnodes; ++i)
  {
    mesh.nodes[3*i]  = ((RealType*)(dbmesh->coords[i]))[0];
    mesh.nodes[3*i+1] = ((RealType*)(dbmesh->coords[i]))[1];
    mesh.nodes[3*i+2] = ((RealType*)(dbmesh->coords[i]))[2];
  }

  // Reconstruct the faces.
  mesh.faces.resize(dbmesh->faces->nfaces);
  int noffset = 0;
  for (int f = 0; f < mesh.faces.size(); ++f)
  {
    mesh.faces[f].resize(dbmesh->faces->shapesize[f]);
    for (int n = 0; n < mesh.faces[f].size(); ++n, ++noffset)
      mesh.faces[f][n] = dbmesh->faces->nodelist[noffset];
  }

  // Reconstruct the "zones" (cells) in terms of the faces.
  mesh.cells.resize(dbmesh->zones->nzones);
  // FIXME
  DBFreeUcdmesh(dbmesh);
  
  // Retrieve the fields.
  DBtoc* contents = DBGetToc(file);
  for (int f = 0; f < contents->nucdvar; ++f)
  {
    string varname = contents->ucdvar_names[f];
    DBucdvar* dbvar = DBGetUcdvar(file, varname.c_str());
    fields[varname] = new RealType[dbvar->nels];
    copy((RealType*)(dbvar->vals[0]), (RealType*)(dbvar->vals[0]) + dbvar->nels,
         fields[varname]);

    // Clean up.
    DBFreeUcdvar(dbvar);
  }

  // Clean up.
  DBClose(file);
}
//-------------------------------------------------------------------

// Explicit instantiation.
template class SiloReader<3, double>;

} // end namespace

#endif

