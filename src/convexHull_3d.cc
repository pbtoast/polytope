//------------------------------------------------------------------------------
// Some explicit instantiations for convexHull_3d.
//------------------------------------------------------------------------------
#include <stdint.h>
#include <stdlib.h>
#include <limits>
#include <iostream>

#include "convexHull_3d.hh"

using namespace std;

namespace polytope {
namespace convexHull_helpers {

// double
template<> double        Point<double>::INF = numeric_limits<double>::max();
template<> Point<double> Point<double>::nil = Point<double>(Point<double>::INF,
                                                            Point<double>::INF,
                                                            Point<double>::INF,
                                                            0,
                                                            0);
template<> Point<double>* Point<double>::NIL = &Point<double>::nil;

}
}
