#include <cstring>
#include "polytope_read_silo.h"
#include "polytope_c.h"
#include "polytope.hh"

#if HAVE_SILO

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

void polytope_read_silo(polytope_tessellation_t* mesh, 
                        int* num_fields,
                        char*** field_names,
                        polytope_real_t*** fields,
                        const char* file_prefix,
                        const char* directory,
                        int cycle,
                        polytope_real_t time,
                        MPI_Comm comm,
                        int num_files,
                        int mpi_tag)
{
  POLY_ASSERT(mesh != NULL);
  map<string, vector<polytope_real_t> > cxxFields;
  string cxxPrefix = file_prefix;
  string cxxDir = directory;

  if (mesh->dimension == 2)
  {
    Tessellation<2, polytope_real_t> cxxMesh;
    SiloReader<2, polytope_real_t>::read(cxxMesh, cxxFields, cxxPrefix, cxxDir,
                                         cycle, time, comm, num_files, mpi_tag);
    fill_tessellation(cxxMesh, mesh);                                           
  }
  else
  {
    POLY_ASSERT(mesh->dimension == 3);
    Tessellation<3, polytope_real_t> cxxMesh;
    SiloReader<3, polytope_real_t>::read(cxxMesh, cxxFields, cxxPrefix, cxxDir,
                                         cycle, time, comm, num_files, mpi_tag);
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
}

}

#endif

