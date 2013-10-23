#ifndef POLYTOPE_C_TESSELLATOR_H
#define POLYTOPE_C_TESSELLATOR_H

#include "polytope_c.h"

#ifdef __cplusplus
extern "C"
{
#endif

// Forward declaration.
//struct polytope_tessellation_t;

// This struct represents a tessellator object that corresponds to the 
// C++ polytope::Tessellator class.
typedef struct polytope_tessellator_t polytope_tessellator_t;

// This function should be called to destroy a tessellator.
void polytope_tessellator_free(polytope_tessellator_t* tessellator);

// Creates an unbounded tessellation of the given points. The coordinates of 
// these points are stored in point-major order and the 0th component of the 
// ith point appears in points[Dimension*i].
void polytope_tessellator_tessellate_unbounded(polytope_tessellator_t* tessellator,
                                               polytope_real_t* points, int num_points,
                                               polytope_tessellation_t* mesh);

// Creates a tessellation bounded by the given bounding box with corners
// (low[0], ..., low[N]) and (high[0], ..., high[N]) for the given set of 
// generator points. The coordinates of these points are stored in point-major 
// order and the 0th component of the ith point appears in points[Dimension*i].
void polytope_tessellator_tessellate_in_box(polytope_tessellator_t* tessellator,
                                            polytope_real_t* points, int num_points,
                                            polytope_real_t* low,
                                            polytope_real_t* high,
                                            polytope_tessellation_t* mesh);

// Generate a Voronoi-like tessellation for the given set of generator 
// points and a description of the geometry (a piecewise linear complex, or 
// PLC) in which they exist. The coordinates of these points are stored in 
// point-major order and the 0th component of the ith point appears in 
// points[Dimension*i].
void polytope_tessellator_tessellate_in_plc(polytope_tessellator_t* tessellator,
                                            polytope_real_t* points, int num_points,
                                            polytope_real_t* plc_points, int num_plc_points,
                                            polytope_plc_t* piecewise_linear_complex,
                                            polytope_tessellation_t* mesh);

// Returns true if the tessellator can tessellate points within a piecewise 
// linear complex (PLC), false if not.
bool polytope_tessellator_handles_plcs(polytope_tessellator_t* tessellator);

// Returns an internal string containing the name of the tessellator.
const char* polytope_tessellator_name(polytope_tessellator_t* tessellator);

// Returns the accuracy to which the tessellator can distinguish coordinates.
polytope_real_t polytope_tessellator_degeneracy(polytope_tessellator_t* tessellator);

// Returns the dimension of the tessellator.
int polytope_tessellator_dimension(polytope_tessellator_t* tessellator);

// Tessellator factories.

#ifdef HAVE_BOOST
// Creates a 2D tessellator using the Boost.Polyhedron implementation.
polytope_tessellator_t* boost_tessellator_new();
#endif

#ifdef HAVE_TRIANGLE
// Creates a 2D tessellator using Triangle.
polytope_tessellator_t* triangle_tessellator_new();
#endif

#ifdef HAVE_TETGEN
// Creates a 3D tessellator using TetGen.
polytope_tessellator_t* tetgen_tessellator_new();
#endif

// Creates a tessellator that uses Voro++ in the given dimension.
polytope_tessellator_t* voroplusplus_tessellator_new(int dimension);

#ifdef __cplusplus
}
#endif

#endif
