#ifndef POLYTOPE_SILO_WRITER_HH
#define POLYTOPE_SILO_WRITER_HH

#ifdef HAVE_SILO

#ifdef HAVE_MPI
#include <mpi.h>
#else
#define MPI_Comm int
#define MPI_COMM_WORLD 0
#endif 

#include "Tessellation.hh"
#include <string>
#include <float.h>
#include <map>

namespace polytope
{

//! \class SiloWriter
//! This class provides a static interface for writing Silo files 
//! containing tessellations made by polytope.
template <int Dimension, typename RealType>
class SiloWriter
{
  // No general recipe
};

//! Partial specialization for 2D tessellations.
template <typename RealType>
class SiloWriter<2, RealType>
{
  public:

  //! Write an arbitrary polyhedral mesh and an associated set of 
  //! cell-centered fields to a SILO file. 
  static void write(const Tessellation<2, RealType>& mesh, 
                    const std::map<std::string, RealType*>& fields,
                    const std::string& filePrefix,
                    int cycle,
                    RealType time,
                    MPI_Comm comm = MPI_COMM_WORLD,
                    int numFiles = 1,
                    int mpiTag = 0);

  //! This version of write omits the cycle and time arguments.
  static void write(const Tessellation<2, RealType>& mesh, 
                    const std::map<std::string, RealType*>& fields,
                    const std::string& filePrefix,
                    MPI_Comm comm = MPI_COMM_WORLD,
                    int numFiles = 1,
                    int mpiTag = 0)
  {
    write(mesh, fields, filePrefix, -1, -FLT_MAX,
          comm, numFiles, mpiTag);
  }

};

//! Partial specialization for 3D tessellations.
template <typename RealType>
class SiloWriter<3, RealType>
{
  public:

  //! Write an arbitrary polyhedral mesh and an associated set of 
  //! cell-centered fields to a SILO file. 
  static void write(const Tessellation<3, RealType>& mesh, 
                    const std::map<std::string, RealType*>& fields,
                    const std::string& filePrefix,
                    int cycle,
                    RealType time,
                    MPI_Comm comm = MPI_COMM_WORLD,
                    int numFiles = 1,
                    int mpiTag = 0);

  //! This version of write omits the cycle and time arguments.
  static void write(const Tessellation<3, RealType>& mesh, 
                    const std::map<std::string, RealType*>& fields,
                    const std::string& filePrefix,
                    MPI_Comm comm = MPI_COMM_WORLD,
                    int numFiles = 1,
                    int mpiTag = 0)
  {
    write(mesh, fields, filePrefix, -1, -FLT_MAX,
          comm, numFiles, mpiTag);
  }

};

} // end namespace

#endif

#endif
