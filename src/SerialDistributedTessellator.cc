//----------------------------------------------------------------------------//
// SerialDistributedTessellator
//----------------------------------------------------------------------------//
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>
#include "mpi.h"

#include "polytope.hh"
#include "deleteCells.hh"
#include "polytope_serialize.hh"
#include "polytope_parallel_utilities.hh"
#include "checkDistributedTessellation.hh"

using namespace std;
using std::min;
using std::max;
using std::abs;

namespace {  // We hide internal functions in an anonymous namespace

//------------------------------------------------------------------------------
// Return the positive index.
//------------------------------------------------------------------------------
inline
int
positiveID(const int x) {
  return x >= 0 ? x : ~x;
}

//------------------------------------------------------------------------------
// Extract the key values from a std::map.
//------------------------------------------------------------------------------
struct ExtractKey {
  template<typename Pair>
  typename Pair::first_type operator()(const Pair& value) {
    return value.first;
  }
};

} // end anonymous namespace

namespace polytope {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
SerialDistributedTessellator<Dimension, RealType>::
SerialDistributedTessellator(Tessellator<Dimension, RealType>* tessellator,
                             bool assumeControl,
                             bool buildCommunicationInfo):
  DistributedTessellator<Dimension, RealType>(tessellator, assumeControl, buildCommunicationInfo) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
SerialDistributedTessellator<Dimension, RealType>::
~SerialDistributedTessellator() {
}

//------------------------------------------------------------------------------
// This method encapsulates the actual distributed tessellation algorithm.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
SerialDistributedTessellator<Dimension, RealType>::
computeDistributedTessellation(const vector<RealType>& points,
                               Tessellation<Dimension, RealType>& mesh) const {

  // Some spiffy shorthand typedefs.
  typedef typename DimensionTraits<Dimension, RealType>::ConvexHull ConvexHull;
  typedef typename DimensionTraits<Dimension, RealType>::CoordHash CoordHash;
  typedef typename DimensionTraits<Dimension, RealType>::Point Point;
  typedef KeyTraits::Key Key;
  const double degeneracy = 1.0e-12;
  
  // Parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Get the full global set of generators from everyone.
  vector<RealType> generators;
  vector<unsigned> genProcOffsets(1, 0);
  const unsigned nlocal = points.size() / Dimension;
  {
    vector<char> localBuffer;
    serialize(points, localBuffer);
    for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
      unsigned bufSize = localBuffer.size();
      MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
      vector<char> buffer = localBuffer;
      buffer.resize(bufSize);
      MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, MPI_COMM_WORLD);
      vector<char>::const_iterator bufItr = buffer.begin();
      deserialize(generators, bufItr, buffer.end());
      POLY_ASSERT(bufItr == buffer.end());
      genProcOffsets.push_back(generators.size()/Dimension);
    }
  }
  POLY_ASSERT(genProcOffsets.size() == numProcs + 1);
  POLY_ASSERT(genProcOffsets.back() == generators.size()/Dimension);

  // Do the global tessellation.
  this->tessellationWrapper(generators, mesh);

  // Find the set of node we share with other processors.
  vector<set<unsigned> > nodeCells = mesh.computeNodeCells();
  map<unsigned, vector<unsigned> > sharedNodes;
  for (unsigned icell = genProcOffsets[rank]; icell != genProcOffsets[rank+1]; ++icell) {
    for (vector<int>::const_iterator faceItr = mesh.cells[icell].begin();
         faceItr != mesh.cells[icell].end();
         ++faceItr) {
      const unsigned iface = positiveID(*faceItr);
      for (vector<unsigned>::const_iterator nodeItr = mesh.faces[iface].begin();
           nodeItr != mesh.faces[iface].end();
           ++nodeItr) {
        for (set<unsigned>::const_iterator otherCellItr = nodeCells[*nodeItr].begin();
             otherCellItr != nodeCells[*nodeItr].end();
             ++otherCellItr) {
          const vector<unsigned>::iterator offItr = lower_bound(genProcOffsets.begin(),
                                                                genProcOffsets.end(),
                                                                *otherCellItr);
          POLY_ASSERT(offItr < genProcOffsets.end());
          const unsigned otherRank = distance(genProcOffsets.begin(), offItr);
          if (otherRank != rank) sharedNodes[otherRank].push_back(*nodeItr);
        }
      }
    }
  }

  // Copy to the mesh shared nodes.
  transform(sharedNodes.begin(), sharedNodes.end(), back_inserter(mesh.neighborDomains), ExtractKey());
  sort(mesh.neighborDomains.begin(), mesh.neighborDomains.end());
  for (typename vector<unsigned>::const_iterator neighborItr = mesh.neighborDomains.begin();
       neighborItr != mesh.neighborDomains.end();
       ++neighborItr) {
    sort(sharedNodes[*neighborItr].begin(), sharedNodes[*neighborItr].end());
    vector<unsigned>::iterator uniqueItr = unique(sharedNodes[*neighborItr].begin(), sharedNodes[*neighborItr].end());
    mesh.sharedNodes.push_back(vector<unsigned>(sharedNodes[*neighborItr].begin(), uniqueItr));
  }

  // Find the set of faces we share with other processors.
  mesh.sharedFaces.resize(mesh.neighborDomains.size());
  for (unsigned icell = genProcOffsets[rank]; icell != genProcOffsets[rank+1]; ++icell) {
    for (vector<int>::const_iterator faceItr = mesh.cells[icell].begin();
         faceItr != mesh.cells[icell].end();
         ++faceItr) {
      const unsigned iface = positiveID(*faceItr);
      POLY_ASSERT(mesh.faceCells[iface].size() == 1 or
                  mesh.faceCells[iface].size() == 2);
      POLY_ASSERT(positiveID(mesh.faceCells[iface][0]) == icell or
                  (mesh.faceCells[iface].size() == 2 and positiveID(mesh.faceCells[iface][1]) == icell));
      if (mesh.faceCells[iface].size() == 2) {
        const unsigned otherCell = (positiveID(mesh.faceCells[iface][0]) == icell ?
                                    positiveID(mesh.faceCells[iface][1]) :
                                    positiveID(mesh.faceCells[iface][0]));
        const vector<unsigned>::iterator offItr = lower_bound(genProcOffsets.begin(),
                                                              genProcOffsets.end(),
                                                              otherCell);
        POLY_ASSERT(offItr < genProcOffsets.end());
        const unsigned otherRank = distance(genProcOffsets.begin(), offItr);
        if (otherRank != rank) {
          const vector<unsigned>::iterator procItr = lower_bound(mesh.neighborDomains.begin(),
                                                                 mesh.neighborDomains.end(),
                                                                 otherRank);
          POLY_ASSERT(procItr < mesh.neighborDomains.end());
          const unsigned j = distance(mesh.neighborDomains.begin(), procItr);
          mesh.sharedFaces[j].push_back(iface);
        }
      }
    }
  }

  // Sort the mesh shared faces in global index order.
  for (unsigned j = 0; j != mesh.sharedFaces.size(); ++j) {
    sort(mesh.sharedFaces[j].begin(), mesh.sharedFaces[j].end());
    POLY_ASSERT(unique(mesh.sharedFaces[j].begin(), mesh.sharedFaces[j].end()) == mesh.sharedFaces[j].end());
  }

  // Now cull the mesh down to just our domain local cells.
  vector<unsigned> mask(generators.size(), 0);
  fill(mask.begin() + genProcOffsets[rank],
       mask.begin() + genProcOffsets[rank + 1],
       1);
  deleteCells(mesh, mask);

  // Post-conditions.
#ifndef NDEBUG
  const string msg = checkDistributedTessellation(mesh);
  if (msg != "ok" and rank == 0) cerr << msg;
  POLY_ASSERT(msg == "ok");
#endif
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class SerialDistributedTessellator<2, double>;
template class SerialDistributedTessellator<3, double>;

}
