#ifndef POLYTOPE_C_TESSELLATION_H
#define POLYTOPE_C_TESSELLATION_H

#ifdef __cplusplus
extern "C"
{
#endif

// This struct holds the data from a tessellation created in polytope.
// It owns all of its data, so it must be destroyed with polytope_tessellation_free.
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

  // Array of node indices: 0 for interior nodes and 1 for nodes on the boundary.
  int num_boundary_nodes;
  unsigned* boundary_nodes;

  // Array of face indices: 0 for interior faces and 1 for faces on the boundadry.
  int num_boundary_faces;
  unsigned* boundary_faces;

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

  // These arrays store the cells connected to nodes in compressed row storage
  // format IF this data was requested from the tessellator--otherwise they 
  // are NULL.
  int* node_cell_offsets;
  int* node_cells;

} polytope_tessellation_t;

// Creates an empty container that will store a tessellation with the given 
// dimension.
polytope_tessellation_t* polytope_tessellation_new(int dimension);

// Writes a human-readable representation of the tessellation to the given string, truncating after n characters.
void polytope_tessellation_snprintf(polytope_tessellation_t* tessellation, char* str, int n);

// Writes a human-readable representation of the tessellation to the given file.
void polytope_tessellation_fprintf(polytope_tessellation_t* tessellation, FILE* stream);

// Clears a tessellation, emptying its contents.
void polytope_tessellation_clear(polytope_tessellation_t* tessellation);

// This function should be called to destroy a tessellation.
void polytope_tessellation_free(polytope_tessellation_t* tessellation);

#ifdef __cplusplus
}
#endif

#endif
