#ifndef POLYTOPE_SILO_WRITER_HH
#define POLYTOPE_SILO_WRITER_HH

#ifdef HAVE_SILO

#include <mpi.h>
#include "Tessellation.hh"
#include <string>
#include <float.h>
#include <map>

namespace polytope
{

//! \class SiloWriter
//! This class provides a static interface for writing Silo files 
//! containing tessellations made by polytope.
template <int Dimension, typename Real>
class SiloWriter
{
  // No general recipe
};

//! Partial specialization for 2D tessellations.
template <typename Real>
class SiloWriter<2, Real>
{
  //! Write an arbitrary polyhedral mesh and an associated set of 
  //! cell-centered fields to a SILO file. 
  static void write(const Tessellation<Dimension>& mesh, 
                    const std::map<std::string, Real*>& fields,
                    const std::string& filePrefix,
                    int cycle,
                    Real time,
                    MPI_Comm comm = MPI_COMM_WORLD,
                    int numFiles = 1,
                    int mpiTag = 0);

  //! This version of write omits the cycle and time arguments.
  static void writeSiloPlot(const Mesh<Dimension>& mesh, 
                            const std::map<std::string, Real*>& scalarFields,
                            const std::map<std::string, Vector<Dimension>*>& vectorFields,
                            const std::string& filePrefix,
                            MPI_Comm comm = MPI_COMM_WORLD,
                            int numFiles = 1,
                            int mpiTag = 0)
  {
    writeSiloPlot(mesh, scalarFields, vectorFields, filePrefix, -1, -FLT_MAX,
                  comm, numFiles, mpiTag);
  }

};

} // end namespace

#endif
#endif
