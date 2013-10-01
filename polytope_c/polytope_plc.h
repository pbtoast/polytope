#ifndef POLYTOPE_C_PLC_H
#define POLYTOPE_C_PLC_H

#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif

// This struct represents a piecewise linear complex corresponding to the 
// C++ class polytope::PLC.
typedef struct polytope_plc_t polytope_plc_t;

// Creates a new empty PLC with the given dimension.
polytope_plc_t* polytope_plc_new(int dimension);

// This function should be called to destroy a PLC.
void polytope_plc_free(polytope_plc_t* plc);

// Adds a facet to the PLC, returning its index.
int polytope_plc_add_facet(polytope_plc_t* plc);

// Adds a node index to the facet to the PLC.
void polytope_plc_add_facet_node(polytope_plc_t* plc, int facet, int node);

// Adds a hole index to the PLC, returning its index.
int polytope_plc_add_hole(polytope_plc_t* plc);

// Adds a facet to the given hole within the PLC, returning its index.
int polytope_plc_add_hole_facet(polytope_plc_t* plc, int hole);

// Adds a node index to the given facet in the given hole within the PLC.
void polytope_plc_add_hole_facet_node(polytope_plc_t* plc, int hole, int facet, int node);

// Returns the number of nodes in the given facet within the PLC.
int polytope_plc_num_facet_nodes(polytope_plc_t* plc, int facet);

// Retrieves the indices of the nodes attached to the given facet.
void polytope_plc_get_facet_nodes(polytope_plc_t* plc, int facet, int* facet_nodes);

// Returns the number of facets in the given hole within the PLC.
int polytope_plc_num_hole_facets(polytope_plc_t* plc, int hole);

// Returns the number of nodes in the given facet of the given hole within the PLC.
int polytope_plc_num_hole_facet_nodes(polytope_plc_t* plc, int hole, int hole_facet);

// Retrieves the indices of nodes attached to the given facet in the given hole within the PLC.
void polytope_plc_get_hole_facet_nodes(polytope_plc_t* plc, int hole, int hole_facet, int* hole_facet_nodes);

// Clears facets and holes to empty the PLC
void polytope_plc_clear(polytope_plc_t* plc);

// Returns true if this PLC is empty, false otherwise.
bool polytope_plc_empty(polytope_plc_t* plc);

// Returns true if this PLC is valid (at first glance), false if it 
// is obviously invalid. This is not a rigorous check!
bool polytope_plc_valid(polytope_plc_t* plc);

// Writes a human-readable representation of the PLC to the given string.
void polytope_plc_snprintf(polytope_plc_t* plc, char* str, int n);

// Writes a human-readable representation of the PLC to the given file.
void polytope_plc_fprintf(polytope_plc_t* plc, FILE* stream);

#ifdef __cplusplus
}
#endif

#endif
