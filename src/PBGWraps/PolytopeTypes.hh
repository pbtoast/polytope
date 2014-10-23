#ifndef __PBGWRAPS_POLYTOPETYPES__
#define __PBGWRAPS_POLYTOPETYPES__

#include "Python.h"
#include "polytope.hh"

#include "PLC.hh"
#include "ReducedPLC.hh"
#include "Tessellation.hh"
#include "Tessellator.hh"
#include "TriangleTessellator.hh"
#include "BoostTessellator.hh"
#include "TetgenTessellator.hh"
#include "MeshEditor.hh"
#include "SiloWriter.hh"

#if HAVE_MPI
#include "DistributedTessellator.hh"
#include "SerialDistributedTessellator.hh"
#endif

//#include "VoroPP_2d.hh"
//#include "VoroPP_3d.hh"

namespace polytope {

//------------------------------------------------------------------------------
// PLC names
//------------------------------------------------------------------------------
typedef PLC<2, double> PLC2d;
typedef PLC<3, double> PLC3d;

typedef ReducedPLC<2, double> ReducedPLC2d;
typedef ReducedPLC<3, double> ReducedPLC3d;

//------------------------------------------------------------------------------
// Tessellation names
//------------------------------------------------------------------------------
typedef Tessellation<2, double> Tessellation2d;
typedef Tessellation<3, double> Tessellation3d;

//------------------------------------------------------------------------------
// Tessellator names
//------------------------------------------------------------------------------
typedef Tessellator<2, double> Tessellator2d;
typedef Tessellator<3, double> Tessellator3d;

#if HAVE_TRIANGLE
typedef TriangleTessellator<double> TriangleTessellator2d;
#endif

#if HAVE_TETGEN
typedef TetgenTessellator TetgenTessellator3d;
#endif

#if HAVE_BOOST_VORONOI
typedef BoostTessellator<double> BoostTessellator2d;
#endif

#if HAVE_MPI
typedef DistributedTessellator<2, double> DistributedTessellator2d;
typedef DistributedTessellator<3, double> DistributedTessellator3d;
typedef SerialDistributedTessellator<2, double> SerialDistributedTessellator2d;
typedef SerialDistributedTessellator<3, double> SerialDistributedTessellator3d;
#endif

// typedef VoroPP_2d<double> VoroTessellator2d;
// typedef VoroPP_3d<double> VoroTessellator3d;

//------------------------------------------------------------------------------
// MeshEditor names
//------------------------------------------------------------------------------
typedef MeshEditor<2, double> MeshEditor2d;
typedef MeshEditor<3, double> MeshEditor3d;


//------------------------------------------------------------------------------
// PLC Getters and Setters
//------------------------------------------------------------------------------

// --- Facets
template<int Dimension, typename RealType>
std::vector<std::vector<int> >* getfacets( PLC<Dimension, RealType>& self )
{ return &self.facets; }
template<int Dimension, typename RealType>
void setfacets( PLC<Dimension, RealType>& self, 
                std::vector<std::vector<int> >* facetsIn )
{ self.facets = *facetsIn; }

// --- Holes
template<int Dimension, typename RealType>
std::vector<std::vector<std::vector<int> > >* getholes( PLC<Dimension, RealType>& self )
{ return &self.holes; }
template<int Dimension, typename RealType>
void setholes( PLC<Dimension, RealType>& self, 
               std::vector<std::vector<std::vector<int> > >* holesIn )
{ self.holes = *holesIn; }


//------------------------------------------------------------------------------
// ReducedPLC Getters and Setters
//------------------------------------------------------------------------------

// --- points
template<int Dimension, typename RealType>
std::vector<RealType>* getpoints( ReducedPLC<Dimension, RealType>& self )
{ return &self.points; }
template<int Dimension, typename RealType>
void setpoints( ReducedPLC<Dimension, RealType>& self, 
                std::vector<RealType>* pointsIn )
{ self.points = *pointsIn; }


//------------------------------------------------------------------------------
// Tessellation Getters and Setters
//------------------------------------------------------------------------------

// --- Nodes
template<int Dimension, typename RealType>
std::vector<RealType>* getnodes( Tessellation<Dimension, RealType>& self )
{ return &self.nodes; }
template<int Dimension, typename RealType>
void setnodes( Tessellation<Dimension, RealType>& self, 
               std::vector<RealType>* nodesIn )
{ self.nodes = *nodesIn; }


// --- Cells
template<int Dimension, typename RealType>
std::vector<std::vector<int> >* getcells( Tessellation<Dimension, RealType>& self )
{ return &self.cells; }
template<int Dimension, typename RealType>
void setcells( Tessellation<Dimension, RealType>& self, 
               std::vector<std::vector<int> >* cellsIn )
{ self.cells = *cellsIn; }


// --- Faces
template<int Dimension, typename RealType>
std::vector<std::vector<unsigned> >* getfaces( Tessellation<Dimension, RealType>& self )
{ return &self.faces; }
template<int Dimension, typename RealType>
void setfaces( Tessellation<Dimension, RealType>& self, 
               std::vector<std::vector<unsigned> >* facesIn )
{ self.faces = *facesIn; }


// --- InfNodes
template<int Dimension, typename RealType>
std::vector<unsigned>* getinfNodes( Tessellation<Dimension, RealType>& self )
{ return &self.infNodes; }
template<int Dimension, typename RealType>
void setinfNodes( Tessellation<Dimension, RealType>& self, 
                  std::vector<unsigned>* infNodesIn )
{ self.infNodes = *infNodesIn; }


// --- faceCells
template<int Dimension, typename RealType>
std::vector<std::vector<int> >* getfaceCells( Tessellation<Dimension, RealType>& self )
{ return &self.faceCells; }
template<int Dimension, typename RealType>
void setfaceCells( Tessellation<Dimension, RealType>& self, 
                   std::vector<std::vector<int> >* faceCellsIn )
{ self.faceCells = *faceCellsIn; }


// --- convexHull
template<int Dimension, typename RealType>
PLC<Dimension, RealType>* getconvexHull( Tessellation<Dimension, RealType>& self )
{ return &self.convexHull; }
template<int Dimension, typename RealType>
void setconvexHull( Tessellation<Dimension, RealType>& self, 
                    PLC<Dimension,RealType>* convexHullIn )
{ self.convexHull = *convexHullIn; }


// --- neighborDomains
template<int Dimension, typename RealType>
std::vector<unsigned>* getneighborDomains( Tessellation<Dimension, RealType>& self )
{ return &self.neighborDomains; }
template<int Dimension, typename RealType>
void setneighborDomains( Tessellation<Dimension, RealType>& self, 
                         std::vector<unsigned>* neighborDomainsIn )
{ self.neighborDomains = *neighborDomainsIn; }


// --- SharedNodes
template<int Dimension, typename RealType>
std::vector<std::vector<unsigned> >* getsharedNodes( Tessellation<Dimension, RealType>& self )
{ return &self.sharedNodes; }
template<int Dimension, typename RealType>
void setsharedNodes( Tessellation<Dimension, RealType>& self, 
                     std::vector<std::vector<unsigned> >* sharedNodesIn )
{ self.sharedNodes = *sharedNodesIn; }


// --- SharedFaces
template<int Dimension, typename RealType>
std::vector<std::vector<unsigned> >* getsharedFaces( Tessellation<Dimension, RealType>& self )
{ return &self.sharedFaces; }
template<int Dimension, typename RealType>
void setsharedFaces( Tessellation<Dimension, RealType>& self, 
                     std::vector<std::vector<unsigned> >* sharedFacesIn )
{ self.sharedFaces = *sharedFacesIn; }

//------------------------------------------------------------------------------
// Wrap writing silo files.
//------------------------------------------------------------------------------
// Start out with a few helper methods in the unnamed namespace.
namespace {
  template<typename RealType>
  std::map<std::string, RealType*>
  buildMapFromPyObject(PyObject* obj, 
                       const unsigned expectedSize) {
    std::map<std::string, RealType*> result;
    if (obj == Py_None) return result;
    if (not PyDict_Check(obj)) {
      std::cerr << "ERROR : need to pass a dictionary of lists with string keys." << std::endl;
      return result;
    }
    char* label;
    PyObject *key, *item, *val;
    Py_ssize_t i = 0, j;
    while (PyDict_Next(obj, &i, &key, &item)) {
      if (not PyString_Check(key)) {
        std::cerr << "ERROR : keys need to be strings." << std::endl;
        return result;
      }
      label = PyString_AsString(key);
      result[std::string(label)] = new RealType[expectedSize];
      if (not PyList_Check(item)) {
        std::cerr << "ERROR : values need to be lists of floats: " << std::string(label) << std::endl;
        return result;
      }
      if (PyList_Size(item) != expectedSize) {
        std::cerr << "ERROR : list of values incorrect size : " << PyList_Size(item) << " != " << expectedSize << std::endl;
        return result;
      }
      for (j = 0; j != expectedSize; ++j) {
        val = PyList_GetItem(item, j);
        if (not PyFloat_Check(val)) {
          std::cerr << "ERROR : list of values should be convertable to floats: " << std::string(label) << std::endl;
          return result;
        }
        result[std::string(label)][j] = PyFloat_AsDouble(val);
      }
    }
    return result;
  }

  // Clean up memory.
  template<typename RealType>
  void
  deleteMemory(std::map<std::string, RealType*>& stuff) {
    for (typename std::map<std::string, RealType*>::iterator itr = stuff.begin();
         itr != stuff.end();
         ++itr) {
      delete [] (itr->second);
    }
  }
}

// Write with just a file prefix.
template<int Dimension, typename RealType>
inline
void
writeTessellation(const Tessellation<Dimension, RealType>& mesh,
                  std::string filePrefix,
                  PyObject* nodeFieldsDict,
                  PyObject* edgeFieldsDict,
                  PyObject* faceFieldsDict,
                  PyObject* cellFieldsDict,
                  int cycle,
                  RealType time) {

#if HAVE_SILO
  // Extract the various optional fields.
  std::map<std::string, RealType*> nodeFields = buildMapFromPyObject<RealType>(nodeFieldsDict, mesh.nodes.size()/Dimension);
  std::map<std::string, RealType*> faceFields = buildMapFromPyObject<RealType>(faceFieldsDict, mesh.faces.size());
  std::map<std::string, RealType*> cellFields = buildMapFromPyObject<RealType>(cellFieldsDict, mesh.cells.size());

  // The edge fields don't currently work in 3D, so screen that case out.
  std::map<std::string, RealType*> edgeFields;
  if (Dimension == 3) {
    if (edgeFieldsDict != Py_None) {
      std::cerr << "WARNING : writing edge fields currently not supported in 3D.  Ignoring passed edgeFields set." << std::endl;
    }
  } else {
    edgeFields = buildMapFromPyObject<RealType>(edgeFieldsDict, mesh.faces.size());
  }

  // Do the deed.
  SiloWriter<Dimension, RealType>::write(mesh,
                                         nodeFields,
                                         edgeFields,
                                         faceFields,
                                         cellFields,
                                         filePrefix,
                                         cycle,
                                         time);

  // Clean up our memory.
  deleteMemory<RealType>(nodeFields);
  deleteMemory<RealType>(edgeFields);
  deleteMemory<RealType>(faceFields);
  deleteMemory<RealType>(cellFields);

#else
  std::cerr << "WARNING : apparently polytope was built without silo support, so cannot write silo files." << std::endl;
#endif
}
  
// Write with a file prefix and directory.
template<int Dimension, typename RealType>
inline
void
writeTessellation(const Tessellation<Dimension, RealType>& mesh,
                  std::string filePrefix,
                  std::string directory,
                  PyObject* nodeFieldsDict,
                  PyObject* edgeFieldsDict,
                  PyObject* faceFieldsDict,
                  PyObject* cellFieldsDict,
                  int cycle,
                  RealType time) {

#if HAVE_SILO
  // Extract the various optional fields.
  std::map<std::string, RealType*> nodeFields = buildMapFromPyObject<RealType>(nodeFieldsDict, mesh.nodes.size()/Dimension);
  std::map<std::string, RealType*> faceFields = buildMapFromPyObject<RealType>(faceFieldsDict, mesh.faces.size());
  std::map<std::string, RealType*> cellFields = buildMapFromPyObject<RealType>(cellFieldsDict, mesh.cells.size());

  // The edge fields don't currently work in 3D, so screen that case out.
  std::map<std::string, RealType*> edgeFields;
  if (Dimension == 3) {
    if (edgeFieldsDict != Py_None) {
      std::cerr << "WARNING : writing edge fields currently not supported in 3D.  Ignoring passed edgeFields set." << std::endl;
    }
  } else {
    edgeFields = buildMapFromPyObject<RealType>(edgeFieldsDict, mesh.faces.size());
  }

  // Do the deed.
  SiloWriter<Dimension, RealType>::write(mesh,
                                         nodeFields,
                                         edgeFields,
                                         faceFields,
                                         cellFields,
                                         filePrefix,
                                         directory,
                                         cycle,
                                         time);

  // Clean up our memory.
  deleteMemory<RealType>(nodeFields);
  deleteMemory<RealType>(edgeFields);
  deleteMemory<RealType>(faceFields);
  deleteMemory<RealType>(cellFields);

#else
  std::cerr << "WARNING : apparently polytope was built without silo support, so cannot write silo files." << std::endl;
#endif
}
  
}

#endif
