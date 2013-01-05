#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "Boundary2D.hh"
#include "Generators.hh"

using namespace std;
using namespace polytope;

// -------------------------------------------------------------------- //
template <class T>
inline std::string to_string(const T& t){
   std::stringstream ss;
   ss << t;
   return ss.str();
}
// -------------------------------------------------------------------- //



// -----------------------------------------------------------------------
// outputResult
// -----------------------------------------------------------------------
void outputResult(Tessellator<2,double>& tessellator,
		  int bType,
		  unsigned nPoints)
{
   Boundary2D<double> boundary;

   boundary.computeDefaultBoundary(bType);
   Generators<2,double> generators( boundary );

   generators.randomPoints( nPoints );

   Tessellation<2,double> mesh;
   tessellator.tessellate(generators.mGenerators,
			  boundary.mGens,
			  boundary.mPLC,
			  mesh);
   POLY_ASSERT( mesh.cells.size() == nPoints );
   
#if HAVE_SILO
   std::string name = "test_RandomPoints_boundary" + 
      to_string(bType) + "_" + to_string(nPoints) + "points";
   vector<double> index( mesh.cells.size());
   for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
   map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
   cellFields["cell_index"] = &index[0];
   polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                          faceFields, cellFields, name);
#endif
}


// -----------------------------------------------------------------------
// tessellate
// -----------------------------------------------------------------------
void tessellate(Boundary2D<double>& boundary,
                Generators<2,double>& generators,
                Tessellator<2,double>& tessellator,
                Tessellation<2,double>& mesh)
{
   if( tessellator.handlesPLCs() ){
      tessellator.tessellate(generators.mGenerators,
                             boundary.mGens, boundary.mPLC, mesh);
   }else{
      tessellator.tessellate(generators.mGenerators,
                             boundary.mLow, boundary.mHigh, mesh);
   }
}


// -----------------------------------------------------------------------
// testBoundary
// -----------------------------------------------------------------------
void testBoundary(Boundary2D<double>& boundary,
                  Tessellator<2,double>& tessellator)
{
   Generators<2,double> generators( boundary );
   unsigned nPoints = 1;
   Tessellation<2,double> mesh;
   for( unsigned n = 0; n < 3; ++n ){
      POLY_ASSERT( mesh.empty() );
      nPoints = nPoints * 10;
      cout << nPoints << " points...";

      generators.randomPoints( nPoints );      
      tessellate(boundary,generators,tessellator,mesh);

      POLY_ASSERT( mesh.cells.size() == nPoints );
      cout << "PASS" << endl;
// #if HAVE_SILO
//       std::string name = "test_RandomPoints_boundary" 
//          + to_string(boundary.mType) + "_" + to_string(nPoints) + "points";
//       vector<double> index( mesh.cells.size());
//       for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
//       map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
//       cellFields["cell_index"] = &index[0];
//       polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
//                                              faceFields, cellFields, name);
// #endif
      mesh.clear();   
   }
}




// -----------------------------------------------------------------------
// testAllBoundaries
// -----------------------------------------------------------------------
void testAllBoundaries(Tessellator<2,double>& tessellator)
{
   for (int bid = 0; bid < 7; ++bid){
      cout << "Testing boundary type " << bid << endl;
      Boundary2D<double> boundary;
      boundary.computeDefaultBoundary(bid);
      testBoundary( boundary, tessellator );
   }
}



// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
   MPI_Init(&argc, &argv);
#endif

   cout << "\nTriangle Tessellator:\n" << endl;
   TriangleTessellator<double> triangle;
   testAllBoundaries(triangle);
   
   cout << "\nVoro 2D Tessellator:\n" << endl;
   VoroPP_2d<double> voro;
   testAllBoundaries(voro);
   
   // int nPoints = 1;
   // for (int i=0; i<4; ++i){
   //    nPoints = nPoints*10;
   //    outputResult(triangle,3,nPoints);
   // }
   

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
