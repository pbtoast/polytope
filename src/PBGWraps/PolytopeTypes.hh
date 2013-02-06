#ifndef __PBGWRAPS_POLYTOPETYPES__
#define __PBGWRAPS_POLYTOPETYPES__

#include "PLC.hh"
#include "Tessellation.hh"
#include "Tessellator.hh"

// TODO: Get rid of this line once I'm building with cmake
#define HAVE_TRIANGLE 1
#define HAVE_TETGEN 1

#if HAVE_TRIANGLE
#include "TriangleTessellator.hh"
#endif

#if HAVE_TETGEN
#include "TetgenTessellator.hh"
#endif

//#include "VoroPP_2d.hh"
//#include "VoroPP_3d.hh"

namespace polytope {

//------------------------------------------------------------------------------
// PLC names
//------------------------------------------------------------------------------
typedef PLC<2, double> PLC2d;
typedef PLC<3, double> PLC3d;

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

// typedef VoroPP_2d<double> VoroTessellator2d;
// typedef VoroPP_3d<double> VoroTessellator3d;

//------------------------------------------------------------------------------
// Tessellation Getters and Setters
//------------------------------------------------------------------------------

// --- Nodes
template<int Dimension, typename RealType>
std::vector<RealType>* getNodes( Tessellation<Dimension, RealType>& self )
{ return &self.nodes; }
template<int Dimension, typename RealType>
void setNodes( Tessellation<Dimension, RealType>& self, 
               std::vector<RealType>* nodesIn )
{ self.nodes = *nodesIn; }


// --- Cells
template<int Dimension, typename RealType>
std::vector<std::vector<int> >* getCells( Tessellation<Dimension, RealType>& self )
{ return &self.cells; }
template<int Dimension, typename RealType>
void setCells( Tessellation<Dimension, RealType>& self, 
               std::vector<std::vector<int> >* cellsIn )
{ self.cells = *cellsIn; }


// --- Faces
template<int Dimension, typename RealType>
std::vector<std::vector<unsigned> >* getFaces( Tessellation<Dimension, RealType>& self )
{ return &self.faces; }
template<int Dimension, typename RealType>
void setFaces( Tessellation<Dimension, RealType>& self, 
               std::vector<std::vector<unsigned> >* facesIn )
{ self.faces = *facesIn; }


// --- InfNodes
template<int Dimension, typename RealType>
std::vector<unsigned>* getInfNodes( Tessellation<Dimension, RealType>& self )
{ return &self.infNodes; }
template<int Dimension, typename RealType>
void setInfNodes( Tessellation<Dimension, RealType>& self, 
                  std::vector<unsigned>* infNodesIn )
{ self.infNodes = *infNodesIn; }


// --- faceCells
template<int Dimension, typename RealType>
std::vector<std::vector<int> >* getFaceCells( Tessellation<Dimension, RealType>& self )
{ return &self.faceCells; }
template<int Dimension, typename RealType>
void setFaceCells( Tessellation<Dimension, RealType>& self, 
                   std::vector<std::vector<int> >* faceCellsIn )
{ self.faceCells = *faceCellsIn; }


// --- convexHull
template<int Dimension, typename RealType>
PLC<Dimension, RealType>* getConvexHull( Tessellation<Dimension, RealType>& self )
{ return &self.convexHull; }
template<int Dimension, typename RealType>
void setConvexHull( Tessellation<Dimension, RealType>& self, 
                    PLC<Dimension,RealType>* convexHullIn )
{ self.convexHull = *convexHullIn; }


// --- neighborDomains
template<int Dimension, typename RealType>
std::vector<unsigned>* getNeighborDomains( Tessellation<Dimension, RealType>& self )
{ return &self.neighborDomains; }
template<int Dimension, typename RealType>
void setNeighborDomains( Tessellation<Dimension, RealType>& self, 
                         std::vector<unsigned>* neighborDomainsIn )
{ self.neighborDomains = *neighborDomainsIn; }


// --- SharedNodes
template<int Dimension, typename RealType>
std::vector<std::vector<unsigned> >* getSharedNodes( Tessellation<Dimension, RealType>& self )
{ return &self.sharedNodes; }
template<int Dimension, typename RealType>
void setSharedNodes( Tessellation<Dimension, RealType>& self, 
                     std::vector<std::vector<unsigned> >* sharedNodesIn )
{ self.sharedNodes = *sharedNodesIn; }


// --- SharedFaces
template<int Dimension, typename RealType>
std::vector<std::vector<unsigned> >* getSharedFaces( Tessellation<Dimension, RealType>& self )
{ return &self.sharedFaces; }
template<int Dimension, typename RealType>
void setSharedFaces( Tessellation<Dimension, RealType>& self, 
                     std::vector<std::vector<unsigned> >* sharedFacesIn )
{ self.sharedFaces = *sharedFacesIn; }


}

#endif
