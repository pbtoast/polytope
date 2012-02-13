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

template<> const std::vector<ConvexHull3d<double>::Point>* ConvexHull3d<double>::Face::points = 0;
template<> const std::vector<ConvexHull3d<int64_t>::Point>* ConvexHull3d<int64_t>::Face::points = 0;

}
}
