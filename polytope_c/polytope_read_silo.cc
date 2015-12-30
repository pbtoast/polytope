#include <cstring>
#include "polytope_read_silo.h"
#include "polytope_c.h"
#include "polytope.hh"

#ifdef HAVE_SILO

namespace polytope
{

// Creates a C tessellation from a C++ one. This template is defined in 
// polytope_tessellator.cc.
template <int Dimension>
void fill_tessellation(const Tessellation<Dimension, polytope_real_t>& t, polytope_tessellation_t* tess);

}

using namespace std;
using namespace polytope;

extern "C"
{

void polytope_read_silo_with_tags(polytope_tessellation_t* mesh, 
                                  int* num_fields,
                                  char*** field_names,
                                  polytope_real_t*** fields,
                                  int* num_node_tags,
                                  char*** node_tag_names,
                                  int** node_tag_sizes,
                                  int*** node_tags,
                                  int* num_edge_tags,
                                  char*** edge_tag_names,
                                  int** edge_tag_sizes,
                                  int*** edge_tags,
                                  int* num_face_tags,
                                  char*** face_tag_names,
                                  int** face_tag_sizes,
                                  int*** face_tags,
                                  int* num_cell_tags,
                                  char*** cell_tag_names,
                                  int** cell_tag_sizes,
                                  int*** cell_tags,
                                  const char* file_prefix,
                                  const char* directory,
                                  int cycle,
                                  polytope_real_t* time,
                                  MPI_Comm comm,
                                  int num_files,
                                  int mpi_tag)
{
  POLY_ASSERT(mesh != NULL);
  map<string, vector<polytope_real_t> > cxxFields;
  map<string, vector<int> > cxxNodeTags, cxxEdgeTags, cxxFaceTags, cxxCellTags;
  string cxxPrefix = file_prefix;
  string cxxDir = directory;

  if (mesh->dimension == 2)
  {
    Tessellation<2, polytope_real_t> cxxMesh;
    SiloReader<2, polytope_real_t>::read(cxxMesh, cxxFields, 
                                         cxxNodeTags, cxxEdgeTags, cxxFaceTags, cxxCellTags,
                                         cxxPrefix, cxxDir,
                                         cycle, *time, comm, num_files, mpi_tag);
    fill_tessellation(cxxMesh, mesh);                                           
  }
  else
  {
    POLY_ASSERT(mesh->dimension == 3);
    Tessellation<3, polytope_real_t> cxxMesh;
    SiloReader<3, polytope_real_t>::read(cxxMesh, cxxFields, 
                                         cxxNodeTags, cxxEdgeTags, cxxFaceTags, cxxCellTags,
                                         cxxPrefix, cxxDir,
                                         cycle, *time, comm, num_files, mpi_tag);
    fill_tessellation(cxxMesh, mesh);                                           
  }

  // Parse the field information.
  int N = (int)cxxFields.size();
  *num_fields = N;
  *field_names = (char**)malloc(sizeof(char*) * N);
  *fields = (polytope_real_t**)malloc(sizeof(polytope_real_t*) * N);
  int index = 0;
  for (map<string, vector<polytope_real_t> >::const_iterator 
       iter = cxxFields.begin(); iter != cxxFields.end(); ++iter)
  {
    *field_names[index] = (char*)malloc(sizeof(char) * (iter->first.length() + 1));
    strcpy(*field_names[index], iter->first.c_str());
    *fields[index] = (polytope_real_t*)malloc(sizeof(polytope_real_t) * iter->second.size());
    copy(iter->second.begin(), iter->second.end(), *fields[index]);
    ++index;
  }

  // Parse the tag information.
  if (node_tags != NULL)
  {
    int N = cxxNodeTags.size();
    *num_node_tags = N;
    *node_tag_names = (char**)malloc(sizeof(char*) * N);
    *node_tag_sizes = (int*)malloc(sizeof(int) * N);
    *node_tags = (int**)malloc(sizeof(int*) * N);
    int index = 0;
    for (map<string, vector<polytope_real_t> >::const_iterator 
         iter = cxxFields.begin(); iter != cxxFields.end(); ++iter)
    {
      *node_tag_names[index] = (char*)malloc(sizeof(char) * (iter->first.length() + 1));
      strcpy(*node_tag_names[index], iter->first.c_str());
      *node_tag_sizes[index] = (int)iter->second.size();
      *node_tags[index] = (int*)malloc(sizeof(int) * iter->second.size());
      copy(iter->second.begin(), iter->second.end(), *node_tags[index]);
      ++index;
    }
  }

  if (edge_tags != NULL)
  {
    int N = cxxEdgeTags.size();
    *num_edge_tags = N;
    *edge_tag_names = (char**)malloc(sizeof(char*) * N);
    *edge_tag_sizes = (int*)malloc(sizeof(int) * N);
    *edge_tags = (int**)malloc(sizeof(int*) * N);
    int index = 0;
    for (map<string, vector<polytope_real_t> >::const_iterator 
         iter = cxxFields.begin(); iter != cxxFields.end(); ++iter)
    {
      *edge_tag_names[index] = (char*)malloc(sizeof(char) * (iter->first.length() + 1));
      strcpy(*edge_tag_names[index], iter->first.c_str());
      *edge_tag_sizes[index] = (int)iter->second.size();
      *edge_tags[index] = (int*)malloc(sizeof(int) * iter->second.size());
      copy(iter->second.begin(), iter->second.end(), *edge_tags[index]);
      ++index;
    }
  }

  if (face_tags != NULL)
  {
    int N = cxxFaceTags.size();
    *num_face_tags = N;
    *face_tag_names = (char**)malloc(sizeof(char*) * N);
    *face_tag_sizes = (int*)malloc(sizeof(int) * N);
    *face_tags = (int**)malloc(sizeof(int*) * N);
    int index = 0;
    for (map<string, vector<polytope_real_t> >::const_iterator 
         iter = cxxFields.begin(); iter != cxxFields.end(); ++iter)
    {
      *face_tag_names[index] = (char*)malloc(sizeof(char) * (iter->first.length() + 1));
      strcpy(*face_tag_names[index], iter->first.c_str());
      *face_tag_sizes[index] = (int)iter->second.size();
      *face_tags[index] = (int*)malloc(sizeof(int) * iter->second.size());
      copy(iter->second.begin(), iter->second.end(), *face_tags[index]);
      ++index;
    }
  }

  if (cell_tags != NULL)
  {
    int N = cxxCellTags.size();
    *num_cell_tags = N;
    *cell_tag_names = (char**)malloc(sizeof(char*) * N);
    *cell_tag_sizes = (int*)malloc(sizeof(int) * N);
    *cell_tags = (int**)malloc(sizeof(int*) * N);
    int index = 0;
    for (map<string, vector<polytope_real_t> >::const_iterator 
         iter = cxxFields.begin(); iter != cxxFields.end(); ++iter)
    {
      *cell_tag_names[index] = (char*)malloc(sizeof(char) * (iter->first.length() + 1));
      strcpy(*cell_tag_names[index], iter->first.c_str());
      *cell_tag_sizes[index] = (int)iter->second.size();
      *cell_tags[index] = (int*)malloc(sizeof(int) * iter->second.size());
      copy(iter->second.begin(), iter->second.end(), *cell_tags[index]);
      ++index;
    }
  }
}

void polytope_read_silo(polytope_tessellation_t* mesh, 
                        int* num_fields,
                        char*** field_names,
                        polytope_real_t*** fields,
                        const char* file_prefix,
                        const char* directory,
                        int cycle,
                        polytope_real_t* time,
                        MPI_Comm comm,
                        int num_files,
                        int mpi_tag)
{
  polytope_read_silo_with_tags(mesh, num_fields, field_names, fields,
                               NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL,
                               file_prefix, directory, cycle, time,
                               comm, num_files, mpi_tag);
}

}

#endif

