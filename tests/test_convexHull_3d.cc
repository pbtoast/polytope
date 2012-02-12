// Try tessellating a simple lattice of generators in a box in parallel.
// We use randomly chosen seed locations to divide up the generators
// between processors.

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <limits>
#include <sstream>

#include "polytope.hh"
#include "convexHull_3d.hh"

#define CHECK(x) if (!(x)) { cout << "FAIL: " << #x << endl; exit(-1); }

using namespace std;

//------------------------------------------------------------------------------
// Define our own local random number generator wrapping the standard srand &
// rand methods.
//------------------------------------------------------------------------------
double random01() {
  return double(rand())/RAND_MAX;
}

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  int n, i;  
  n = 100;

  polytope::convexHull_helpers::Point<double> *P = new polytope::convexHull_helpers::Point<double>[n];  // input
  for (i = 0; i < n; i++) { 
    P[i].x = double(rand())/RAND_MAX;
    P[i].y = double(rand())/RAND_MAX;
    P[i].z = double(rand())/RAND_MAX;
  }

  polytope::convexHull_helpers::Point<double> *list = sort(P, n);
  polytope::convexHull_helpers::Point<double> **A = new polytope::convexHull_helpers::Point<double> *[2*n], **B = new polytope::convexHull_helpers::Point<double> *[2*n];
  polytope::convexHull_helpers::lowerHull(list, n, A, B);

  for (i = 0; A[i] != polytope::convexHull_helpers::Point<double>::NIL; A[i++]->act())  // output 
    std::cout << A[i]->prev-P << " " << A[i]-P << " " << A[i]->next-P << "\n";
  delete A;  delete B;  delete P;

  cout << "PASS" << endl;
  return 0;
}
