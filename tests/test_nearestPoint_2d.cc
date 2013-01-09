// 2D nearestPoint on PLC unit test.

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <limits>
#include <sstream>
#include <ctime>

#include "polytope.hh"
#include "nearestPoint.hh"

#define POLY_CHECK(x) if (!(x)) { cout << "FAIL: " << #x << endl; exit(-1); }

using namespace std;

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  // Create a square PLC with one hole.
  const unsigned numVertices = 8;
  double vertices[16] = {-5.0, -5.0,
                          5.0, -5.0,
                          5.0,  5.0,
                         -5.0,  5.0,
                         -1.0, -1.0,
                          1.0, -1.0,
                          1.0,  1.0,
                         -1.0,  1.0};
  PLC<2, double> plc;
  plc.facets.resize(4);
  plc.holes.resize(1);
  for (unsigned i = 0; i != 4; ++i) {
    plc.facets[i].resize(2);
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i + 1) % 4;
    plc.holes[0][i].resize(2);
    plc.holes[0][i][0] = 4 + i;
    plc.holes[0][i][1] = 4 + (i + 1) % 4;
  }

  // Test some points.
  const double tol = 1.0e-10;
  double[2] p = {-10.0, -10.0}, result[2], answer[2];
  nearestPoint(p, numVertices, vertices, plc, result);
  answer = {-5.0, -5.0};
  POLY_CHECK(distance<2, double>(result, answer) < tol);

  cout << "PASS" << endl;
  return 0;
}
