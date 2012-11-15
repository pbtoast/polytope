//------------------------------------------------------------------------------
// A collection of random stuff useful for testing in polytope.
//------------------------------------------------------------------------------
#ifndef __polytope_test_utilities__
#define __polytope_test_utilities__

//------------------------------------------------------------------------------
// A macro for checking true/false test conditions.
//------------------------------------------------------------------------------
#define CHECK(x) if (!(x)) { cout << "FAIL: " << #x << endl; exit(-1); }
#define CHECK2(x, msg) if (!(x)) { cout << "FAIL: " << #x << endl << msg << endl; exit(-1); }

//------------------------------------------------------------------------------
// Define our own local random number generator wrapping the standard srand &
// rand methods.
//------------------------------------------------------------------------------
double random01() {
  return double(rand())/RAND_MAX;
}

#endif
