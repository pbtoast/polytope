#include "polytope_write_silo.h"
#include "polytope_c.h"
#include "polytope.hh"

#ifdef HAVE_SILO

#ifdef HAVE_MPI
#include <mpi.h>
#else
#ifndef MPI_Comm
#define MPI_Comm int
#define MPI_COMM_WORLD 0
#endif
#endif 

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

void polytope_write_silo_with_tags(polytope_tessellation_t* mesh, 
                                   int num_node_fields,
                                   char** node_field_names,
                                   polytope_real_t** node_fields,
                                   int num_node_tags,
                                   char** node_tag_names,
                                   int* node_tag_sizes,
                                   int** node_tags,
                                   int num_edge_fields,
                                   char** edge_field_names,
                                   polytope_real_t** edge_fields,
                                   int num_edge_tags,
                                   char** edge_tag_names,
                                   int* edge_tag_sizes,
                                   int** edge_tags,
                                   int num_face_fields,
                                   char** face_field_names,
                                   polytope_real_t** face_fields,
                                   int num_face_tags,
                                   char** face_tag_names,
                                   int* face_tag_sizes,
                                   int** face_tags,
                                   int num_cell_fields,
                                   char** cell_field_names,
                                   polytope_real_t** cell_fields,
                                   int num_cell_tags,
                                   char** cell_tag_names,
                                   int* cell_tag_sizes,
                                   int** cell_tags,
                                   const char* file_prefix,
                                   const char* directory,
                                   int cycle,
                                   polytope_real_t time,
                                   MPI_Comm comm,
                                   int num_files,
                                   int mpi_tag)
{
  POLY_ASSERT(mesh != NULL);
  string cxxPrefix = file_prefix;
  string cxxDir = directory;

  // Fill fields.
  map<string, polytope_real_t*> cxxNodeFields, cxxEdgeFields, cxxFaceFields, cxxCellFields;
  for (int i = 0; i < num_node_fields; ++i)
    cxxNodeFields[node_field_names[i]] = node_fields[i];
  for (int i = 0; i < num_edge_fields; ++i)
    cxxEdgeFields[edge_field_names[i]] = edge_fields[i];
  for (int i = 0; i < num_face_fields; ++i)
    cxxFaceFields[face_field_names[i]] = face_fields[i];
  for (int i = 0; i < num_cell_fields; ++i)
    cxxCellFields[cell_field_names[i]] = cell_fields[i];

  // Fill tags.
  map<string, vector<int>*> cxxNodeTags, cxxEdgeTags, cxxFaceTags, cxxCellTags;
  for (int i = 0; i < num_node_tags; ++i)
  {
    vector<int>* tag = new vector<int>(node_tag_sizes[i]);
    copy(node_tags[i], node_tags[i] + node_tag_sizes[i], tag->begin());
    cxxNodeTags[node_tag_names[i]] = tag;
  }
  for (int i = 0; i < num_edge_tags; ++i)
  {
    vector<int>* tag = new vector<int>(edge_tag_sizes[i]);
    copy(edge_tags[i], edge_tags[i] + edge_tag_sizes[i], tag->begin());
    cxxEdgeTags[edge_tag_names[i]] = tag;
  }
  for (int i = 0; i < num_face_tags; ++i)
  {
    vector<int>* tag = new vector<int>(face_tag_sizes[i]);
    copy(face_tags[i], face_tags[i] + face_tag_sizes[i], tag->begin());
    cxxFaceTags[face_tag_names[i]] = tag;
  }
  for (int i = 0; i < num_cell_tags; ++i)
  {
    vector<int>* tag = new vector<int>(cell_tag_sizes[i]);
    copy(cell_tags[i], cell_tags[i] + cell_tag_sizes[i], tag->begin());
    cxxCellTags[cell_tag_names[i]] = tag;
  }

  if (mesh->dimension == 2)
  {
    Tessellation<2, polytope_real_t> cxxMesh;
    fill_tessellation(mesh, cxxMesh);
    SiloWriter<2, polytope_real_t>::write(cxxMesh, 
                                          cxxNodeFields, cxxNodeTags, 
                                          cxxEdgeFields, cxxEdgeTags,
                                          cxxFaceFields, cxxFaceTags,
                                          cxxCellFields, cxxCellTags,
                                          cxxPrefix, cxxDir, cycle, 
                                          time, comm, num_files, mpi_tag);
  }
  else
  {
    POLY_ASSERT(mesh->dimension == 3);
    Tessellation<3, polytope_real_t> cxxMesh;
    fill_tessellation(mesh, cxxMesh);
    SiloWriter<3, polytope_real_t>::write(cxxMesh, 
                                          cxxNodeFields, cxxNodeTags, 
                                          cxxEdgeFields, cxxEdgeTags,
                                          cxxFaceFields, cxxFaceTags,
                                          cxxCellFields, cxxCellTags,
                                          cxxPrefix, cxxDir, cycle, 
                                          time, comm, num_files, mpi_tag);
  }

  // Clean up.
  for (map<string, vector<int>*>::iterator
       iter = cxxNodeTags.begin(); iter != cxxNodeTags.end(); ++iter)
  {
    delete iter->second;
  }
  for (map<string, vector<int>*>::iterator
       iter = cxxEdgeTags.begin(); iter != cxxEdgeTags.end(); ++iter)
  {
    delete iter->second;
  }
  for (map<string, vector<int>*>::iterator
       iter = cxxFaceTags.begin(); iter != cxxFaceTags.end(); ++iter)
  {
    delete iter->second;
  }
  for (map<string, vector<int>*>::iterator
       iter = cxxCellTags.begin(); iter != cxxCellTags.end(); ++iter)
  {
    delete iter->second;
  }
}

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
  polytope_write_silo_with_tags(mesh, 
                                num_node_fields, node_field_names, node_fields,
                                0, NULL, NULL, NULL,
                                num_edge_fields, edge_field_names, edge_fields,
                                0, NULL, NULL, NULL,
                                num_face_fields, face_field_names, face_fields,
                                0, NULL, NULL, NULL,
                                num_cell_fields, cell_field_names, cell_fields,
                                0, NULL, NULL, NULL,
                                file_prefix, directory, cycle, time, comm,
                                num_files, mpi_tag);
}


}

#endif

