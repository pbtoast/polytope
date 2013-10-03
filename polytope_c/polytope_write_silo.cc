#include "polytope_write_silo.h"
#include "polytope_c.h"
#include "polytope.hh"

#if HAVE_SILO

namespace polytope
{

// Creates a C++ tessellation from a C one. This template is defined in 
// polytope_tessellator.cc.
template <int Dimension>
void fill_tessellation(polytope_tessellation_t* tess, Tessellation<Dimension, polytope_real_t>& t);

}

using namespace std;
using namespace polytope;

extern "C"
{

void polytope_write_silo(polytope_tessellation_t* mesh, 
                         int num_node_fields,
                         char** node_field_names,
                         polytope_real_t** node_fields,
                         int num_edge_fields,
                         char** edge_field_names,
                         polytope_real_t** edge_fields,
                         int num_face_fields,
                         char** face_field_names,
                         polytope_real_t** face_fields,
                         int num_cell_fields,
                         char** cell_field_names,
                         polytope_real_t** cell_fields,
                         const char* file_prefix,
                         const char* directory,
                         int cycle,
                         polytope_real_t time,
                         MPI_Comm comm,
                         int num_files,
                         int mpi_tag)
{
  POLY_ASSERT(mesh != NULL);
  map<string, polytope_real_t*> cxxNodeFields, cxxEdgeFields, cxxFaceFields, cxxCellFields;
  string cxxPrefix = file_prefix;
  string cxxDir = directory;

  // Fill fields.
  for (int i = 0; i < num_node_fields; ++i)
    cxxNodeFields[node_field_names[i]] = node_fields[i];
  for (int i = 0; i < num_edge_fields; ++i)
    cxxEdgeFields[edge_field_names[i]] = edge_fields[i];
  for (int i = 0; i < num_face_fields; ++i)
    cxxFaceFields[face_field_names[i]] = face_fields[i];
  for (int i = 0; i < num_cell_fields; ++i)
    cxxCellFields[cell_field_names[i]] = cell_fields[i];

  if (mesh->dimension == 2)
  {
    Tessellation<2, polytope_real_t> cxxMesh;
    fill_tessellation(mesh, cxxMesh);
    SiloWriter<2, polytope_real_t>::write(cxxMesh, cxxNodeFields, cxxEdgeFields, cxxFaceFields, cxxCellFields,
                                          cxxPrefix, cxxDir, cycle, time, comm, num_files, mpi_tag);
  }
  else
  {
    POLY_ASSERT(mesh->dimension == 3);
    Tessellation<3, polytope_real_t> cxxMesh;
    fill_tessellation(mesh, cxxMesh);
    SiloWriter<3, polytope_real_t>::write(cxxMesh, cxxNodeFields, cxxEdgeFields, cxxFaceFields, cxxCellFields,
                                          cxxPrefix, cxxDir, cycle, time, comm, num_files, mpi_tag);
  }
}

}

#endif

