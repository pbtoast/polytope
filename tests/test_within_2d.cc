// 2D within on PLC unit test.

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <limits>
#include <sstream>
#include <ctime>

#include "polytope.hh"
#include "within.hh"
#include "polytope_test_utilities.hh"

using namespace std;

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  // Create a PLC with one hole:
  //  ________    ________
  //  |       \  /       |
  //  |        \/        |
  //  |      __  __      |
  //  |      | \/ |      |
  //  |      |____|      |
  //  |__________________|
  //
  const unsigned numVertices = 10;
  double vertices[20] = {-5.0, -5.0,
                          5.0, -5.0,
                          5.0,  5.0,
                          0.0,  0.0,
                         -5.0,  5.0,
                         -1.0, -4.0,
                          1.0, -4.0,
                          1.0, -1.0,
                          0.0, -2.5,
                         -1.0, -1.0};
  polytope::PLC<2, double> plc;
  plc.facets.resize(5);
  plc.holes.resize(1);
  plc.holes[0].resize(5);
  for (unsigned i = 0; i != 5; ++i) {
    plc.facets[i].resize(2);
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i + 1) % 5;
    plc.holes[0][i].resize(2);
    plc.holes[0][i][0] = 5 + i;
    plc.holes[0][i][1] = 5 + (i + 1) % 5;
  }

  // Test some points.  

  { // Inside: fully
    double p[2] = {-4.0, 1.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK( inside == answer );
  }

  { // Outside: fully
    double p[2] = {-6.0, -3.0};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK( inside == answer );
  }

  { // Outside: within hole
    double p[2] = {0.0, -3.5};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK( inside == answer );
  }

  { // Inside: tangent to vertex
    double p[2] = {3.0, 0.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK( inside == answer );
  }
  
  { // Outside: tangent to vertex
    double p[2] = {0.0, 5.0};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK( inside == answer );
  }

  { // Inside: tangent to hole vertex
    double p[2] = {2.0, -2.5};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK( inside == answer );
  }

  { // Outside: within hole, tangent to hole vertex
    double p[2] = {-0.5, -2.5};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK( inside == answer );
  }

  { // Inside: tangent to face of hole
    double p[2] = {-1.0, 0.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK( inside == answer );
  }

  cout << "PASS" << endl;
  return 0;
}
