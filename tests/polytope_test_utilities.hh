//------------------------------------------------------------------------------
// A collection of random stuff useful for testing in polytope.
//------------------------------------------------------------------------------
#ifndef __polytope_test_utilities__
#define __polytope_test_utilities__

#include "polytope.hh"
#include "Boundary2D.hh"

//------------------------------------------------------------------------------
// A macro for checking true/false test conditions.
//------------------------------------------------------------------------------
#define POLY_CHECK(x) if (!(x)) { cout << "FAIL: " << #x << endl; exit(-1); }
#define POLY_CHECK2(x, msg) if (!(x)) { cout << "FAIL: " << #x << endl << msg << endl; exit(-1); }

//------------------------------------------------------------------------------
// Define our own local random number generator wrapping the standard srand &
// rand methods.
//------------------------------------------------------------------------------
double random01() {
  return double(rand())/RAND_MAX;
}

//------------------------------------------------------------------------------
// Wrapper to tessellate a 2D boundary for both VoroPP_2D and Triangle
//------------------------------------------------------------------------------
template<typename RealType>
void tessellate2D(std::vector<RealType>& points,
                  Boundary2D<RealType>& boundary,
                  Tessellator<2,RealType>& tessellator,
                  Tessellation<2,RealType>& mesh){
   if( tessellator.handlesPLCs() ){
      tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
   }else{
      tessellator.tessellate(points, boundary.mLow, boundary.mHigh, mesh);
   }
}


#endif
