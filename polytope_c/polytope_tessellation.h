#ifndef POLYTOPE_C_TESSELLATION_H
#define POLYTOPE_C_TESSELLATION_H

#include "polytope_c.h"
#include "polytope_plc.h"

#ifdef __cplusplus
extern "C"
{
#endif

// This struct represents a tessellation corresponding to the C++ 
// polytope::Tessellation class.
typedef struct polytope_tessellation_t polytope_tessellation_t;

// This function creates an empty tessellation with the given dimension.
polytope_tessellation_t* polytope_tessellation_new(int dimension);

// This function should be called to destroy a tessellation.
void polytope_tessellation_free(polytope_tessellation_t* tessellation);

// Returns the dimension of the tessellation.
int polytope_tessellation_dimension(polytope_tessellation_t* tessellation);

// Clear the tessellation, emptying it of all data.
void polytope_tessellation_clear(polytope_tessellation_t* tessellation);

// Returns true if the tessellation is empty (not defined), false otherwise.
bool polytope_tessellation_empty(polytope_tessellation_t* tessellation);

// Returns the number of nodes in the tessellation.
int polytope_tessellation_num_nodes(polytope_tessellation_t* tessellation);

// Retrieves the array of (Dimension*numNodes) node position coordinates for 
// the tessellation. The components are stored in node-major order and 
// the 0th component of the ith node appears in node_coords[Dimension*i].
void polytope_tessellation_get_nodes(polytope_tessellation_t* tessellation,
                                     polytope_real_t* node_coords);

// Returns the number of cells in the tessellation.
int polytope_tessellation_num_cells(polytope_tessellation_t* tessellation);

// Retrieves the cell-face connectivity of the tessellation in a compressed-
// row-storage format.
void polytope_get_num_cell_faces(polytope_tessellation_t* tessellation,
                                 int* num_cell_faces);

void polytope_get_cell_faces(polytope_tessellation_t* tessellation,
                             int* cell_face_indices);

// Writes a human-readable representation of the tessellation to the given file.
void polytope_tessellation_fprintf(polytope_tessellation_t* tessellation, FILE* stream);

int polytope_tessellation_num_inf_nodes(polytope_tessellation_t* tessellation);

// Retrieves the array of node indices: 0 for interior nodes and 1 for nodes at
// "infinity" if this is an unbounded tessellation. The infinite node
// is the termination point on a spherical surface of a ray going
// out to infinity.
void polytope_tessellation_get_inf_nodes(polytope_tessellation_t* tessellation, unsigned* inf_nodes);

int polytope_tessellation_num_inf_faces(polytope_tessellation_t* tessellation);

// Returns the array of face indices: 0 for interior faces and 1 for faces at
// "infinity" if this is an unbounded tessellation. The infinite face connects
// the collection of infinite nodes for a given unbounded cell.
void polytope_tessellation_get_inf_faces(polytope_tessellation_t* tessellation, unsigned* inf_faces);

// Retrieves the array of cell indices for each face, i.e., the cells that share
// the face. There are two entries per face: the first and second cells attached 
// to the face. An index of -1 as the second cell index indicates that the 
// face lies on a boundary of the tessellation.
void polytope_tessellation_get_face_cells(polytope_tessellation_t* tessellation, int* face_cells);

// Returns a newly-created piecewise linear complex containing the 
// convex hull of the point distribution. Not all Tessellators hand back the convex 
// hull, so this may be empty, in which case you must compute the convex 
// hull yourself.
polytope_plc_t* polytope_tessellation_convex_hull(polytope_tessellation_t* tessellation);

// In the case of a parallel calculation, this function returns the number of 
// neighbor domains this portion of the tessellation is in contact with. In a 
// serial calculation, this returns 0.
int polytope_tessellation_num_neighbor_domains(polytope_tessellation_t* tessellation);

// In the case of a parallel calculation, this function retrieves an 
// array containing the set of neighbor domains this portion of the 
// tessellation is in contact with. In a serial calculation, this 
// returns NULL.
void polytope_tessellation_get_neighbor_domains(polytope_tessellation_t* tessellation, unsigned* domains);

// In the case of a parallel calculation, this function returns the number of 
// nodes this tessellation interacts with on the given neighboring domain.
int polytope_tessellation_num_neighbor_nodes(polytope_tessellation_t* tessellation, unsigned domain);

// In the case of a parallel calculation, this function returns an internal
// array containing the indices of nodes that are shared with the given 
// neighbor domain.
void polytope_tessellation_get_neighbor_nodes(polytope_tessellation_t* tessellation, unsigned domain, unsigned* nodes);

// In the case of a parallel calculation, this function returns the number of 
// faces this tessellation interacts with on the given neighboring domain.
int polytope_tessellation_num_neighbor_faces(polytope_tessellation_t* tessellation, unsigned domain);

// In the case of a parallel calculation, this function returns an internal
// array containing the indices of faces that are shared with the given 
// neighbor domain.
void polytope_tessellation_get_neighbor_faces(polytope_tessellation_t* tessellation, unsigned domain, unsigned* faces);

// Retrieves arrays containing the indices of the cells 
// that touch the given mesh node. (CRS)
void polytope_get_num_node_cells(polytope_tessellation_t* tessellation,
                                 int* num_node_cells);
void polytope_get_node_cells(polytope_tessellation_t* tessellation,
                             int* node_cells);

#ifdef __cplusplus
}
#endif

#endif
