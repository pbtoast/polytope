#ifndef __PBGWRAPS_POLYTOPETYPES__
#define __PBGWRAPS_POLYTOPETYPES__

#include "polytope.hh"

#include "PLC.hh"
#include "Tessellation.hh"
#include "Tessellator.hh"
#include "TriangleTessellator.hh"
#include "TetgenTessellator.hh"

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

}

#endif
