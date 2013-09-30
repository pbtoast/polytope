#ifndef POLYTOPE_C_TESSELLATION_H
#define POLYTOPE_C_TESSELLATION_H

#include "polytope_c.h"

#ifdef __cplusplus
extern "C"
{
#endif

// This struct represents a proxy for the C++ polytope::Tessellation class
// that can be queried accordingly.
typedef struct polytope_tessellation_t polytope_tessellation_t;

// This struct holds the data from a tessellation created in polytope, and 
// can be emitted by a polytope_tessellation object. All of its data is 
// copied, so it must be destroyed accordingly.
typedef struct 
{
  // The dimension of the tessellation.
  int dimension;

  // An array of (Dimension*numNodes) values containing components of 
  // node positions. The components are stored in node-major order and 
  // the 0th component of the ith node appears in nodes[Dimension*i].
  int num_nodes; 
  polytope_real_t* nodes;

  // These arrays define the cell-face topology of the mesh in Compressed 
  // Row Storage format. That is, cell_offsets[i] gives the offset within 
  // the cell_faces array that denotes the first face of the ith cell.
  int num_cells;
  int* cell_offsets;
  int* cell_faces;

  // These arrays defines the topology of the faces of the mesh in Compressed
  // Row Storage format. A face has an arbitrary number of nodes in 3D and 2 
  // nodes in 2D. face_offsets[i] gives the offset within the face_nodes 
  // array that denotes the first node of the ith face. Nodes for a given face 
  // are arranged counterclockwise around the face viewed from the 
  // "positive" (outside) direction. 
  int num_faces;
  int* face_offsets;
  unsigned* face_nodes;

  // Array of node indices: 0 for interior nodes and 1 for nodes at
  // "infinity" if this is an unbounded tessellation. The infinite node
  // is the termination point on a spherical surface of a ray going
  // out to infinity.
  int num_inf_nodes;
  unsigned* inf_nodes;

  // Array of face indices: 0 for interior faces and 1 for faces at
  // "infinity" for unbounded tessellations. The infinite face connects
  // the collection of infinite nodes for a given unbounded cell.
  int num_inf_faces;
  unsigned* inf_faces;

  // An array of cell indices for each face, i.e., the cells that share
  // the face. face_cells[2*i] is the index of the first cell for the ith 
  // face, and face_cells[2*i+1] is the index of the second cell, or -1 if 
  // that face lies on the boundary of the tessellation.
  int* face_cells;

  // The convex hull, if it was computed by the tessellator.
  polytope_plc_t* convex_hull;

  // In parallel calculations, the set of neighbor domains this portion of
  // the tessellation is in contact with.
  int num_neighbor_domains;
  unsigned* neighbor_domains;

  // In parallel calculations, the nodes and faces this domain shares with
  // each neighbor domain. These arrays are stored in compressed row storage
  // format like those above.
  // NOTE: we implicitly assume that any domains of rank less than ours we
  //       are receiving from, while any domains of greater rank we send
  //       to.
  int* shared_node_domain_offsets;
  unsigned* shared_nodes;
  int* shared_face_domain_offsets;
  unsigned* shared_faces;

} polytope_tessellation_data_t;

// This function should be called to destroy a tessellation_data object.
void polytope_tessellation_data_free(polytope_tessellation_data_t* tessellation_data);

// Writes a human-readable representation of the tessellation to the given file.
void polytope_tessellation_fprintf(polytope_tessellation_t* tessellation, FILE* stream);

// This function should be called to destroy a tessellation.
void polytope_tessellation_free(polytope_tessellation_t* tessellation);

// Emits a newly-allocated struct containing the data within the tessellation. 
// This struct must be freed with tessellation_data_free().
polytope_tessellation_data_t* polytope_tessellation_data(polytope_tessellation_t* tessellation);

// Retrieves arrays containing the indices of the cells 
// that touch the given mesh node. (CRS format)
void polytope_get_num_node_cells(polytope_tessellation_t* tessellation,
                                 int* num_node_cells);
void polytope_get_node_cells(polytope_tessellation_t* tessellation,
                             int* node_cells);

#ifdef __cplusplus
}
#endif

#endif
