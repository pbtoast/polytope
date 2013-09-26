#ifndef POLYTOPE_C_SILO_WRITER_H
#define POLYTOPE_C_SILO_WRITER_H

#if HAVE_SILO

#ifdef __cplusplus
extern "C"
{
#endif

// Writes an arbitrary polyhedral mesh and an associated set of 
// (node, edge, face, cell)-centered fields to a SILO file in the given directory.
static void polytope_write_silo(const Tessellation<3, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    RealType time,
                    MPI_Comm comm = MPI_COMM_WORLD,
                    int numFiles = -1,
                    int mpiTag = 0);

  //! Write an arbitrary polyhedral mesh and an associated set of 
  //! (node, edge, face, cell)-centered fields to a SILO file. This version generates a 
  //! directory name automatically. For parallel runs, the directory 
  //! name is filePrefix-nproc. For serial runs, the directory is 
  //! the current working directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void write(const Tessellation<3, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::string& filePrefix,
                    int cycle,
                    RealType time,
                    MPI_Comm comm = MPI_COMM_WORLD,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, "", cycle, time, comm, numFiles, mpiTag);
  }

  //! This version of write omits the cycle and time arguments.
  static void write(const Tessellation<3, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    MPI_Comm comm = MPI_COMM_WORLD,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, directory, -1, -FLT_MAX,
          comm, numFiles, mpiTag);
  }

  //! This version of write omits the cycle and time arguments and 
  //! automatically generates the directory name from the file prefix.
  static void write(const Tessellation<3, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::string& filePrefix,
                    MPI_Comm comm = MPI_COMM_WORLD,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, "", -1, -FLT_MAX,
          comm, numFiles, mpiTag);
  }

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
void polytope_plc_add_hole_facet_node(polytope_plc_t* plc, int facet, int node);

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

// Writes a human-readable representation of the PLC to the given file.
void polytope_plc_fprintf(polytope_plc_t* plc, FILE* stream);

#ifdef __cplusplus
}
#endif

#endif

#endif
