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
// Traverse the given points of a polygonal facet along their convex
// hull, writing their indices to indices.
//-------------------------------------------------------------------
template <typename RealType>
void 
traverseConvexHull(const vector<RealType>& points,
                   vector<int>& indices)
{
  // Find the "lowest" point in the set.
  RealType ymin = FLT_MAX;
  int index0 = -1;
  int numPoints = points.size() / 2;
  for (int p = 0; p < numPoints; ++p)
  {
    if (ymin > points[2*p+1])
    {
      ymin = points[2*p+1];
      index0 = p;
    }
  }

  // We start with this point and a horizontal angle.
  RealType thetaPrev = 0.0;
  indices.push_back(index0);

  // Now start gift wrapping.
  int i = index0;
  do 
  {
    RealType dthetaMin = 2.0*M_PI;
    int jMin = -1;
    for (int j = 0; j < numPoints; ++j)
    {
      if (j != i)
      {
        RealType dx = points[2*j] - points[2*i],
             dy = points[2*j+1] - points[2*i+1];
        RealType theta = atan2(dy, dx);
        RealType dtheta = theta - thetaPrev;
        if (dtheta < 0.0)
          dtheta += 2.0*M_PI;
        if (dthetaMin > dtheta)
        {
          dthetaMin = dtheta;
          jMin = j;
        }
      }
    }
    if (jMin != index0)
      indices.push_back(jMin);
    thetaPrev += dthetaMin;
    i = jMin;
  }
  while (i != index0);

  // The convex hull should be a polygon unless the input points 
  // don't form a polygon.
  ASSERT((numPoints <= 2) or 
         ((numPoints > 2) and (indices.size() > 2)));
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
SiloWriter<3, RealType>::
write(const Tessellation<3, RealType>& mesh, 
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
    DIR* masterDir = opendir(masterDirName.c_str());
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
    snprintf(filename, 1024, "%s/%s-d.silo", groupdirname, prefix.c_str(), cycle);
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
  char *coordnames[3];
  coordnames[0] = (char*)"xcoords";
  coordnames[1] = (char*)"ycoords";
  coordnames[2] = (char*)"zcoords";

  // Node coordinates.
  int numNodes = mesh.nodes.size() / 3;
  vector<double> x(numNodes), y(numNodes), z(numNodes);
  for (int i = 0; i < numNodes; ++i)
  {
    x[i] = mesh.nodes[3*i];
    y[i] = mesh.nodes[3*i+1];
    z[i] = mesh.nodes[3*i+2];
  }
  double* coords[3];
  coords[0] = &(x[0]);
  coords[1] = &(y[0]);
  coords[2] = &(z[0]);

  // Figure out face-node connectivity. We do this by computing centers
  // for all the cells and then using them to define face normals, 
  // from which node orderings can be determining using a convex hull 
  // determination algorithm (gift wrapping).
  int numCells = mesh.cells.size();
  int numFaces = mesh.faces.size();
  vector<RealType> cellCenters(numCells);
  for (int c = 0; c < numCells; ++c)
  {
    const vector<int>& cellFaces = mesh.cells[c];
    int numNodes = 0;
    for (int f = 0; f < cellFaces.size(); ++f)
    {
      const vector<unsigned>& faceNodes = mesh.faces[cellFaces[f]];
      for (int n = 0; n < faceNodes.size(); ++n)
      {
        cellCenters[3*c]   += mesh.nodes[3*faceNodes[n]];
        cellCenters[3*c+1] += mesh.nodes[3*faceNodes[n]+1];
        cellCenters[3*c+2] += mesh.nodes[3*faceNodes[n]+2];
        ++numNodes;
      }
    }
    cellCenters[3*c+0] /= numNodes;
    cellCenters[3*c+1] /= numNodes;
    cellCenters[3*c+2] /= numNodes;
  }
  vector<int> faceNodeCounts(numFaces), 
              allFaceNodes;
  for (int f = 0; f < mesh.faces.size(); ++f)
  {
    const vector<unsigned>& faceNodes = mesh.faces[f];

    // Compute the normal vector for the face, pointing outward from 
    // its first cell.
    ASSERT(faceNodes.size() >= 3);
    double faceCenter[3];
    for (int n = 0; n < faceNodes.size(); ++n)
    {
      faceCenter[0] += mesh.nodes[3*faceNodes[n]];
      faceCenter[1] += mesh.nodes[3*faceNodes[n]+1];
      faceCenter[2] += mesh.nodes[3*faceNodes[n]+2];
    }
    faceCenter[0] /= faceNodes.size();
    faceCenter[1] /= faceNodes.size();
    faceCenter[2] /= faceNodes.size();

    // Construct vectors v1, v2, and v3, where v1 is the vector pointing from the 
    // face center to the first face node, v2 is a vector pointing from the face 
    // center to any other face, node, and v3 is their cross product.
    double v1[3];
    for (int d = 0; d < 3; ++d)
      v1[d] = mesh.nodes[3*faceNodes[0]+d] - faceCenter[d];

    double v2[3], normal[3], normalMag;
    for (int n = 1; n < faceNodes.size(); ++n)
    {
      for (int d = 0; d < 3; ++d)
        v2[d] = mesh.nodes[3*faceNodes[n]+d] - faceCenter[d];

      // normal = v1 x v2.
      normal[0] = v1[1]*v2[2] - v1[2]*v2[1];
      normal[1] = v1[2]*v2[0] - v1[0]*v2[2];
      normal[1] = v1[0]*v2[1] - v1[1]*v2[0];
      normalMag = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
      if (normalMag > 1e-14) break;
    }
    normal[0] /= normalMag; normal[1] /= normalMag; normal[2] /= normalMag;

    double v3[3], cellCenter[3];
    for (int d = 0; d < 3; ++d)
    {
      cellCenter[d] = cellCenters[3*mesh.faceCells[f][0]+d];
      v3[d] = faceCenter[d] - cellCenter[d];
    }
    if ((normal[0]*v3[0] + normal[1]*v3[1] + normal[2]*v3[2]) < 0.0)
    {
      normal[0] *= -1.0; normal[1] *= -1.0; normal[2] *= -1.0;
    }

    // Now project the coordinates of the face's nodes to the plane
    // with the given normal and centered about the face center.
    vector<RealType> points(2*faceNodes.size()); // NOTE: planar coordinates (2D)
    double e1[3], e2[3]; // Basis vectors in the plane.
    double v1Mag = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    for (int d = 0; d < 3; ++d)
      e1[d] = v1[d] / v1Mag;

    // e2 = normal x e1.
    e2[0] = normal[1]*e1[2] - normal[2]*e1[1];
    e2[1] = normal[2]*e1[0] - normal[0]*e1[2];
    e2[1] = normal[0]*e1[1] - normal[1]*e1[0];
    for (int p = 0; p < points.size(); ++p)
    {
      // v = node center - cell center.
      double v[3];
      for (int d = 0; d < 3; ++d)
        v[d] = mesh.nodes[3*faceNodes[p]+d] - cellCenter[d];

      // Compute the perpendicular component of the point
      // with location v:
      // vPerp = v - (n o v)n.
      double vPerp[3];
      for (int d = 0; d < 3; ++d)
        vPerp[d] = v[d] - (normal[0]*v[0] + normal[1]*v[1] + normal[2]*v[2]) * normal[d];

      // Project it to the plane.
      points[2*p]   = vPerp[0]*e1[0] + vPerp[1]*e1[1] + vPerp[2]*e1[2];
      points[2*p+1] = vPerp[0]*e2[0] + vPerp[1]*e2[1] + vPerp[2]*e2[2];
    }

    // Find the node order by traversing the convex hull of 
    // the points within the plane, appending them to allFaceNodes.
    vector<int> indices;
    traverseConvexHull(points, indices);
    faceNodeCounts[f] = indices.size();
    for (int n = 0; n < indices.size(); ++n)
      allFaceNodes.push_back(faceNodes[indices[n]]);
  }

  // Figure out cell-face connectivity.
  vector<int> cellFaceCounts(numCells), 
              allCellFaces;
  for (int c = 0; c < numCells; ++c)
  {
    const vector<int>& cellFaces = mesh.cells[c];
    cellFaceCounts[c] = cellFaces.size();
    allCellFaces.insert(allCellFaces.end(), 
                        cellFaces.begin(), cellFaces.end());
  }

  // The polyhedral zone list is referred to in the options list.
  DBAddOption(optlist, DBOPT_PHZONELIST, (char*)"mesh_zonelist");

  // Write out the 3D polyhedral mesh.
  DBPutUcdmesh(file, (char*)"mesh", 3, coordnames, coords,
               numNodes, numCells, 0, 0,
               DB_DOUBLE, optlist); 

  // Write the connectivity information.
  DBPutPHZonelist(file, (char*)"mesh_zonelist", 
                  faceNodeCounts.size(), &faceNodeCounts[0], 
                  allFaceNodes.size(), &allFaceNodes[0], 0, 
                  cellFaceCounts.size(), &cellFaceCounts[0],
                  allCellFaces.size(), &allCellFaces[0], 
                  0, 0, numCells-1, optlist);

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
    vector<RealType> xdata(mesh.numCells()), ydata(mesh.numCells()),
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
template class SiloWriter<3, double>;

} // end namespace

#endif
