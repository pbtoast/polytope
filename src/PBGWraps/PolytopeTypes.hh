#ifndef __PBGWRAPS_POLYTOPETYPES__
#define __PBGWRAPS_POLYTOPETYPES__

#include "PLC.hh"
#include "Tessellation.hh"
#include "Tessellator.hh"

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
// A New Section
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
std::vector<RealType>*
getNodes( Tessellation<Dimension, RealType>& self ){
   return &self.nodes;
}

template<int Dimension, typename RealType>
void
setNodes( Tessellation<Dimension, RealType>& self, std::vector<RealType>& nodesIn ){
   self.nodes = nodesIn;
}


}

#endif
