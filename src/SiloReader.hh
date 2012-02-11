#ifndef POLYTOPE_SILO_READER_HH
#define POLYTOPE_SILO_READER_HH

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

//! \class SiloReader
//! This class provides a static interface for reading Silo files 
//! containing tessellations made by polytope.
template <int Dimension, typename RealType>
class SiloReader
{
  // No general recipe
};

//! Partial specialization for 2D tessellations.
template <typename RealType>
class SiloReader<2, RealType>
{
  public:

  //! Read an arbitrary polygonal mesh and an associated set of 
  //! cell-centered fields to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void read(Tessellation<2, RealType>& mesh, 
                   std::map<std::string, RealType*>& fields,
                   const std::string& filePrefix,
                   const std::string& directory,
                   int cycle,
                   RealType& time,
                   MPI_Comm comm = MPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0);

  //! Read an arbitrary polygonal mesh and an associated set of 
  //! cell-centered fields to a SILO file. This version generates a 
  //! directory name automatically. For parallel runs, the directory 
  //! name is filePrefix-nproc. For serial runs, the directory is 
  //! the current working directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void read(Tessellation<2, RealType>& mesh, 
                   std::map<std::string, RealType*>& fields,
                   const std::string& filePrefix,
                   int cycle,
                   RealType& time,
                   MPI_Comm comm = MPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    read(mesh, fields, filePrefix, "", cycle, time, comm, numFiles, mpiTag);
  }

  //! This version of read omits the cycle and time arguments.
  static void read(Tessellation<2, RealType>& mesh, 
                   std::map<std::string, RealType*>& fields,
                   const std::string& filePrefix,
                   const std::string& directory,
                   MPI_Comm comm = MPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    RealType time;
    read(mesh, fields, filePrefix, directory, -1, time,
         comm, numFiles, mpiTag);
  }

  //! This version of read omits the cycle and time arguments and 
  //! automatically generates the directory name from the file prefix.
  static void read(Tessellation<2, RealType>& mesh, 
                   std::map<std::string, RealType*>& fields,
                   const std::string& filePrefix,
                   MPI_Comm comm = MPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    RealType time;
    read(mesh, fields, filePrefix, "", -1, time,
         comm, numFiles, mpiTag);
  }

};

//! Partial specialization for 3D tessellations.
template <typename RealType>
class SiloReader<3, RealType>
{
  public:

  //! Read an arbitrary polyhedral mesh and an associated set of 
  //! cell-centered fields to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void read(Tessellation<3, RealType>& mesh, 
                   std::map<std::string, RealType*>& fields,
                   const std::string& filePrefix,
                   const std::string& directory,
                   int cycle,
                   RealType& time,
                   MPI_Comm comm = MPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0);

  //! Read an arbitrary polyhedral mesh and an associated set of 
  //! cell-centered fields to a SILO file. This version generates a 
  //! directory name automatically. For parallel runs, the directory 
  //! name is filePrefix-nproc. For serial runs, the directory is 
  //! the current working directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void read(Tessellation<3, RealType>& mesh, 
                   std::map<std::string, RealType*>& fields,
                   const std::string& filePrefix,
                   int cycle,
                   RealType& time,
                   MPI_Comm comm = MPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    read(mesh, fields, filePrefix, "", cycle, time, comm, numFiles, mpiTag);
  }

  //! This version of read omits the cycle and time arguments.
  static void read(Tessellation<3, RealType>& mesh, 
                   std::map<std::string, RealType*>& fields,
                   const std::string& filePrefix,
                   const std::string& directory,
                   MPI_Comm comm = MPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    RealType time;
    read(mesh, fields, filePrefix, directory, -1, time,
         comm, numFiles, mpiTag);
  }

  //! This version of read omits the cycle and time arguments and 
  //! automatically generates the directory name from the file prefix.
  static void read(Tessellation<3, RealType>& mesh, 
                   std::map<std::string, RealType*>& fields,
                   const std::string& filePrefix,
                   MPI_Comm comm = MPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    RealType time;
    read(mesh, fields, filePrefix, "", -1, time,
         comm, numFiles, mpiTag);
  }

};

} // end namespace

#endif

#endif
