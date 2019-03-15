#ifndef POLYTOPE_CHECKDISTRIBUTEDTESSELLATION_HH
#define POLYTOPE_CHECKDISTRIBUTEDTESSELLATION_HH

// This header provides a validity checker for the parallel data stored in a 
// tessellation.  This is not necessarily cheap to call!
// Adapated from the Spheral consistency methods.
// We assume MPI is available here, so don't include this header if that's not true.
#include <vector>
#include <sstream>
#include <algorithm>

#include "polytope.hh"
#include "Point.hh"
#include "polytope_geometric_utilities.hh"
#include "polytope_parallel_utilities.hh"
#include "DimensionTraits.hh"

#include "mpi.h"

namespace { // unnamed

//------------------------------------------------------------------------------
// Look for anyone who has a non-empty string.
//------------------------------------------------------------------------------
inline
std::string
reduceToMaxString(const std::string& x,
                  const unsigned rank,
                  const unsigned numDomains) {
  unsigned badRank = polytope::allReduce((x.size() == 0 ? numDomains : rank),
                                         MPI_MIN,
                                         MPI_COMM_WORLD);
  if (badRank == numDomains) {
    return x;
  } else {
    unsigned size = x.size();
    MPI_Bcast(&size, 1, MPI_UNSIGNED, badRank, MPI_COMM_WORLD);
    POLY_ASSERT(size > 0);
    std::string result = x;
    result.resize(size, '0');
    MPI_Bcast(&result[0], size, MPI_CHAR, badRank, MPI_COMM_WORLD);
    return result;
  }
}

//------------------------------------------------------------------------------
// Check that all domains agree about the set of shared elements between them.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
std::string
checkConsistentCommInfo(const std::string& label,
                        const unsigned rank,
                        const unsigned numDomains,
                        const std::vector<typename polytope::DimensionTraits<Dimension, RealType>::IntPoint>& hashes,
                        const std::vector<unsigned>& neighborDomains,
                        const std::vector<std::vector<unsigned> >& sharedIDs) {

  typedef typename polytope::DimensionTraits<Dimension, RealType>::IntPoint Point;

  std::string result = "";

  // Go over all domains!
  for (unsigned sendProc = 0; sendProc != numDomains; ++sendProc) {

    // Tell everyone how many domains we're checking.
    unsigned numChecks = neighborDomains.size();
    MPI_Bcast(&numChecks, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);

    // For each neighbor to be checked:
    for (unsigned k = 0; k != numChecks; ++k) {

      // ID of the neighbor domain to be checked.
      unsigned recvProc;
      if (rank == sendProc) recvProc = neighborDomains[k];
      MPI_Bcast(&recvProc, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);

      // Send the hashes the send proc thinks it shares with the receive.
      unsigned bufSize;
      std::vector<char> buffer;
      if (rank == sendProc) {
        polytope::serialize(unsigned(sharedIDs[k].size()), buffer);
        for (unsigned i = 0; i != sharedIDs[k].size(); ++i) {
          POLY_ASSERT(sharedIDs[k][i] < hashes.size());
          polytope::serialize(hashes[sharedIDs[k][i]], buffer);
        }
        bufSize = buffer.size();
      }
      MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
      std::vector<Point> recvHashes;
      buffer.resize(bufSize);
      MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, MPI_COMM_WORLD);
      std::vector<char>::const_iterator itr = buffer.begin();
      unsigned nother;
      polytope::deserialize(nother, itr, buffer.end());
      recvHashes.resize(nother);
      for (unsigned i = 0; i != nother; ++i) polytope::deserialize(recvHashes[i], itr, buffer.end());
      POLY_ASSERT(itr == buffer.end());

      // If we're the receive processor, do we agree?
      if (rank == recvProc and result == "") {
        std::vector<unsigned>::const_iterator neighbItr = find(neighborDomains.begin(), neighborDomains.end(), sendProc);
        if (neighbItr == neighborDomains.end()) {
          std::stringstream os;
          os << "Domain " << sendProc << " thinks it is sending to " << rank << ", but " << rank << " disagrees.";
          result = os.str();
        } else {
          const unsigned kk = std::distance(neighborDomains.begin(), neighbItr);
          if (recvHashes.size() != sharedIDs[kk].size()) {
            std::stringstream os;
            os << "Domain " << sendProc << " -> " << rank << " disagree about the number of shared hashes:  " << nother << " != " << sharedIDs[k].size();
            result = os.str();
          } else {
            unsigned i = 0;
            while (i != nother and result == "") {
              if (recvHashes[i] != hashes[sharedIDs[kk][i]]) {
                std::stringstream os;
                os << "Domain " << sendProc << " -> " << rank << " disagree about shared hashes." << std::endl;
                for (unsigned ii = 0; ii != nother; ++ii) {
                  os << "    " << recvHashes[ii] << " " << hashes[sharedIDs[kk][ii]] << std::endl;
                }
                result = os.str();
              }
              ++i;
            }
          }
        }
      }
    }

    if (result != "") std::cout << result << std::endl;

    // Globally reduce the message so far.
    result = reduceToMaxString(result, rank, numDomains);
    if (result != "") return label + " : " + result;
  }
  
  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Check that all shared elements have been found.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
std::string
checkAllSharedElementsFound(const std::string& label,
                            const unsigned rank,
                            const unsigned numDomains,
                            const std::vector<typename polytope::DimensionTraits<Dimension, RealType>::IntPoint>& hashes,
                            const std::vector<unsigned>& neighborDomains,
                            const std::vector<std::vector<unsigned> >& sharedIDs,
                            const bool doNotAllowEmptySets) {

  typedef typename polytope::DimensionTraits<Dimension, RealType>::IntPoint Point;

  std::string result = "";

  // Create the reverse mapping of hashes to local IDs.
  std::map<Point, unsigned> hash2id;
  for (unsigned i = 0; i != hashes.size(); ++i) hash2id[hashes[i]] = i;

  // Sort and pack up our hashes.
  std::vector<Point> sortedHashes = hashes;
  std::sort(sortedHashes.begin(), sortedHashes.end());
  std::vector<char> localBuffer;
  polytope::serialize(sortedHashes, localBuffer);

  // Go over all domains!
  for (unsigned sendProc = 0; sendProc != numDomains; ++sendProc) {

    // Broadcast the send processors full set of hashes to everyone.
    unsigned bufSize = localBuffer.size();
    MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
    std::vector<char> buffer = localBuffer;
    buffer.resize(bufSize);
    MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, MPI_COMM_WORLD);
    std::vector<Point> otherHashes;
    std::vector<char>::const_iterator bufItr = buffer.begin();
    polytope::deserialize(otherHashes, bufItr, buffer.end());

    if (rank != sendProc) {

      // Find the intersection of the remote set of hashes with our own.
      std::vector<Point> commonHashes;
      std::set_intersection(sortedHashes.begin(), sortedHashes.end(),
                            otherHashes.begin(), otherHashes.end(),
                            std::back_inserter(commonHashes));

      // Do we have an entry for this neighbor?
      std::vector<unsigned>::const_iterator neighbItr = find(neighborDomains.begin(), neighborDomains.end(), sendProc);

      // Did we erroneously entirely miss this neighbor?
      if (result == "" and
          neighbItr == neighborDomains.end() and commonHashes.size() > 0) {
        std::stringstream os;
        os << "Domain " << sendProc << " should be sending to " << rank << ", but " << rank << " is missing it.";
        result = os.str();
      }

      // Do we erroneously list intersections with this neighbor?
      if (result == "" and 
          doNotAllowEmptySets and neighbItr != neighborDomains.end() and commonHashes.size() == 0) {
        std::stringstream os;
        os << "Domain " << sendProc << " is  sending to " << rank << ", but there are no common hashes.";
        result = os.str();
      }

      // The rest of our tests only apply if there are common hashes.
      if (result == "" and
          neighbItr != neighborDomains.end() and commonHashes.size() > 0) {
        const unsigned k = std::distance(neighborDomains.begin(), neighbItr);
        
        // Did we find the right number of common points?
        if (result == "" and
            sharedIDs[k].size() != commonHashes.size()) {
          std::stringstream os;
          os << "Domain " << sendProc << " shares " << commonHashes.size()
             << " positions overlaps with domain " << rank << ", but " << rank << " has found "
             << sharedIDs[k].size() << " values.  Common hashes: ";
          for (typename std::vector<Point>::const_iterator itr = commonHashes.begin();
               itr != commonHashes.end();
               ++itr) {
            os << " " << (*itr);
          }
          result = os.str();
        }

        // Are the correct set of hashes present?
        if (result == "") {
          std::vector<Point> myHashes;
          for (std::vector<unsigned>::const_iterator itr = sharedIDs[k].begin();
               itr != sharedIDs[k].end();
               ++itr) myHashes.push_back(hashes[*itr]);
          std::sort(myHashes.begin(), myHashes.end());
          if (myHashes != commonHashes) {
            std::stringstream os;
            os << "Domains " << sendProc << " <-> " << rank << " do not have the correct set of shared hashes." << std::endl
               << sendProc << " : ";
            std::copy(sortedHashes.begin(), sortedHashes.end(), std::ostream_iterator<Point>(os, " "));
            os << std::endl << rank << " : ";
            std::copy(commonHashes.begin(), commonHashes.end(), std::ostream_iterator<Point>(os, " "));
            os << std::endl;
            result = os.str();
          }
        }
      }
    }

    // Globally reduce the message so far.
    result = reduceToMaxString(result, rank, numDomains);
    if (result != "") return label + " : " + result;
  }

  // That's it.
  return result;
}

}           // unnamed

namespace polytope {

//------------------------------------------------------------------------------
//! \fn string checkDistributedTessellation(const Tessellation& mesh)
//! \brief Checks the parallel data structures in a tessellation for consistency.
//! \param mesh The tessellation we're checking.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
std::string
checkDistributedTessellation(const Tessellation<Dimension, RealType>& mesh) {

  typedef DimensionTraits<Dimension, RealType> Traits;
  typedef typename Traits::IntPoint Point;

  // Parallel configuration.
  int rank, numDomains;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numDomains);

  // Prepare our result optimistically.
  std::string result = "";

  // Compute the bounding box for normalizing our coordinates.
  RealType xmin[Dimension], xmax[Dimension];
  geometry::computeBoundingBox<Dimension, RealType>(mesh.nodes, true, xmin, xmax);
  const RealType dxhash = Traits::maxLength(xmin, xmax)/(KeyTraits::maxKey1d - KeyTraits::minKey1d);

  // First check that all processors agree about who is talking to whom.
  for (unsigned sendProc = 0; sendProc != numDomains; ++sendProc) {
    unsigned numNeighbors = mesh.neighborDomains.size();
    std::vector<unsigned> checkNeighbors = mesh.neighborDomains;
    MPI_Bcast(&numNeighbors, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
    if (numNeighbors > 0) {
      checkNeighbors.resize(numNeighbors);
      MPI_Bcast(&checkNeighbors.front(), numNeighbors, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
      if (not (std::binary_search(checkNeighbors.begin(), checkNeighbors.end(), rank) == 
               std::binary_search(mesh.neighborDomains.begin(), mesh.neighborDomains.end(), sendProc))) {
        result = "Processors don't agree about who is talking to whom!";
      }
    }
  }
  result = reduceToMaxString(result, rank, numDomains);

  // Hash the mesh node positions.
  std::vector<Point> nodeHashes;
  POLY_ASSERT(mesh.nodes.size() % Dimension == 0);
  const unsigned numNodes = mesh.nodes.size()/Dimension;
  for (unsigned i = 0; i != numNodes; ++i) nodeHashes.push_back(Traits::constructPoint(&mesh.nodes[Dimension*i], xmin, dxhash, i));

  // Check the communicated nodes.
  if (result == "") result = checkConsistentCommInfo<Dimension, RealType>("Node", rank, numDomains, nodeHashes, mesh.neighborDomains, mesh.sharedNodes);
  if (result == "") result = checkAllSharedElementsFound<Dimension, RealType>("Node", rank, numDomains, nodeHashes, mesh.neighborDomains, mesh.sharedNodes, true);

  // Something weird about hashing the faces, so I'm suspending those checks for now.

  // // Hash the mesh face positions.
  // std::vector<Point> faceHashes;
  // const unsigned numFaces = mesh.faces.size();
  // for (unsigned i = 0; i != numFaces; ++i) faceHashes.push_back(Traits::faceCentroid(mesh, i, xmin, dxhash));

  // // Check the communicated faces.
  // if (result == "ok") result = checkConsistentCommInfo<Dimension, RealType>("Face", rank, numDomains, faceHashes, mesh.neighborDomains, mesh.sharedFaces);
  // if (result == "ok") result = checkAllSharedElementsFound<Dimension, RealType>("Face", rank, numDomains, faceHashes, mesh.neighborDomains, mesh.sharedFaces, false);

  if( result == "" ) result = "ok";
  return result;
}

}

#endif
