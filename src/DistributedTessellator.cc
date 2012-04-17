//----------------------------------------------------------------------------//
// DistributedTessellator
//----------------------------------------------------------------------------//
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>
#include <limits>
#include "mpi.h"

#include "polytope.hh"
#include "Point.hh"
#include "polytope_serialize.hh"
#include "ReducedPLC.hh"
#include "convexHull_2d.hh"
#include "convexHull_3d.hh"
#include "convexIntersect.hh"
#include "deleteCells.hh"
#include "bisectSearch.hh"

using namespace std;
using std::min;
using std::max;
using std::abs;

namespace {  // We hide internal functions in an anonymous namespace

//------------------------------------------------------------------------------
// Hide dimensional specializations.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType> struct DimensionTraits {};

// 2D
template<typename RealType>
struct DimensionTraits<2, RealType> {
  typedef typename polytope::ReducedPLC<2, RealType> ConvexHull;
  typedef uint64_t CoordHash;
  typedef polytope::Point2<CoordHash> Point;

  static ConvexHull convexHull(const vector<RealType>& points, 
                               const RealType* low,
                               const RealType& dx) { 
    return ConvexHull(polytope::convexHull_2d(points, low, dx), points);
  }
  static Point constructPoint(const RealType* ri,
                              const RealType* rlow,
                              const RealType& dx,
                              const size_t i) {
    return Point(ri[0] - rlow[0], ri[1] - rlow[1], dx, i);
  }
};

// 3D
template<typename RealType>
struct DimensionTraits<3, RealType> {
  typedef typename polytope::ReducedPLC<3, RealType> ConvexHull;
  typedef uint64_t CoordHash;
  typedef polytope::Point3<CoordHash> Point;

  static ConvexHull convexHull(const vector<RealType>& points, 
                               const RealType* low,
                               const RealType& dx) { 
    return ConvexHull(polytope::convexHull_3d(points, low, dx), points);
  }
  static Point constructPoint(const RealType* ri,
                              const RealType* rlow,
                              const RealType& dx,
                              const size_t i) {
    return Point(ri[0] - rlow[0], ri[1] - rlow[1], ri[2] - rlow[2], dx, i);
  }
};

// Traits to handle mapping RealType -> MPI data type.
template<typename RealType> struct DataTypeTraits {};

template<> struct DataTypeTraits<float> { 
  static MPI_Datatype MpiDataType() { return MPI_FLOAT; }
};

template<> struct DataTypeTraits<double> { 
  static MPI_Datatype MpiDataType() { return MPI_DOUBLE; }
};

//------------------------------------------------------------------------------
// Comparator to compare std::pair's by their first or second element.
//------------------------------------------------------------------------------
template<typename T1, typename T2>
struct ComparePairByFirstElement {
  bool operator()(const std::pair<T1, T2>& lhs, const std::pair<T1, T2>& rhs) const {
    return lhs.first < rhs.first;
  }
};

template<typename T1, typename T2>
struct ComparePairBySecondElement {
  bool operator()(const std::pair<T1, T2>& lhs, const std::pair<T1, T2>& rhs) const {
    return lhs.second < rhs.second;
  }
};

} // end anonymous namespace

namespace polytope {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
DistributedTessellator<Dimension, RealType>::
DistributedTessellator(Tessellator<Dimension, RealType>* tessellator,
                       bool assumeControl,
                       bool buildCommunicationInfo):
  mSerialTessellator(tessellator),
  mAssumeControl(assumeControl),
  mBuildCommunicationInfo(buildCommunicationInfo),
  mType(unbounded),
  mLow(0),
  mHigh(0),
  mPLCptr(0) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
DistributedTessellator<Dimension, RealType>::
~DistributedTessellator() {
  if (mAssumeControl)
    delete mSerialTessellator;
}

//------------------------------------------------------------------------------
// Compute an unbounded tessellation.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
tessellate(const vector<RealType>& points,
           Tessellation<Dimension, RealType>& mesh) const {
  mType = unbounded;
  this->computeDistributedTessellation(points, mesh);
}

//------------------------------------------------------------------------------
// Compute the tessellation in the box.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
tessellate(const vector<RealType>& points,
           RealType* low,
           RealType* high,
           Tessellation<Dimension, RealType>& mesh) const {
  mType = box;
  mLow = low;
  mHigh = high;
  this->computeDistributedTessellation(points, mesh);
}

//------------------------------------------------------------------------------
// Compute the tessellation in a PLC.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
tessellate(const vector<RealType>& points,
           const PLC<Dimension, RealType>& geometry,
           Tessellation<Dimension, RealType>& mesh) const {
  mType = plc;
  mPLCptr = &geometry;
  this->computeDistributedTessellation(points, mesh);
}

//------------------------------------------------------------------------------
// This method encapsulates the actual distributed tessellation algorithm.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
computeDistributedTessellation(const vector<RealType>& points,
                               Tessellation<Dimension, RealType>& mesh) const {

  // Some spiffy shorthand typedefs.
  typedef typename DimensionTraits<Dimension, RealType>::ConvexHull ConvexHull;
  typedef typename DimensionTraits<Dimension, RealType>::CoordHash CoordHash;
  typedef typename DimensionTraits<Dimension, RealType>::Point Point;
  const double degeneracy = 1.0e-10;
  
  // Parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Compute the bounding box for normalizing our coordinates.
  RealType rlow[Dimension], rhigh[Dimension];
  this->computeBoundingBox(points, rlow, rhigh);
  RealType fscale = 0;
  for (unsigned i = 0; i != Dimension; ++i) fscale = max(fscale, rhigh[i] - rlow[i]);
  ASSERT(fscale > 0);
  fscale = 1.0/fscale;

  // Our goal is to build up the necessary set of generators to completely 
  // specify the tessellation for all the points passed on this processor.
  // Start by copying the input set in normalized coordinates.
  const unsigned nlocal = points.size() / Dimension;
  vector<RealType> generators;
  generators.reserve(points.size());
  for (unsigned i = 0; i != nlocal; ++i) {
    for (unsigned j = 0; j != Dimension; ++j) {
      generators.push_back((points[Dimension*i + j] - rlow[j])*fscale);
    }
  }
  ASSERT(generators.size() == points.size());
  ASSERT(nlocal == 0 or *min_element(generators.begin(), generators.end()) >= 0);
  ASSERT(nlocal == 0 or *max_element(generators.begin(), generators.end()) <= 1);

  // We can skip a lot of work if there's only one domain!
  if (numProcs > 1) {

    // Compute the convex hull of each domain and distribute them to all processes.
    const ConvexHull localHull = DimensionTraits<Dimension, RealType>::convexHull(generators, rlow, degeneracy);
    vector<ConvexHull> domainHulls(numProcs);
    vector<unsigned> domainCellOffset(1, 0);
    {
      vector<char> localBuffer;
      serialize(localHull, localBuffer);
      for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
        vector<char> buffer = localBuffer;
        unsigned bufSize = localBuffer.size();
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
        buffer.resize(bufSize);
        MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, MPI_COMM_WORLD);
        vector<char>::const_iterator itr = buffer.begin();
        deserialize(domainHulls[sendProc], itr, buffer.end());
        ASSERT(itr == buffer.end());
        domainCellOffset.push_back(domainCellOffset.back() + domainHulls[sendProc].points.size()/Dimension);
      }
    }
    ASSERT(domainHulls.size() == numProcs);
    ASSERT(domainCellOffset.size() == numProcs + 1);

    // Create a tessellation of the hull vertices for all domains.
    vector<RealType> hullGenerators;
    for (unsigned i = 0; i != numProcs; ++i) {
      copy(domainHulls[i].points.begin(), domainHulls[i].points.end(), back_inserter(hullGenerators));
    }
    ASSERT(hullGenerators.size()/Dimension == domainCellOffset.back());
    Tessellation<Dimension, RealType> hullMesh;
    this->tessellationWrapper(hullGenerators, hullMesh);

    // Find the set of domains we need to communicate with according to two criteria:
    // 1.  Any domain hull that intersects our own.
    // 2.  Any domain hull that has elements adjacent to one of ours in the hullMesh.
    set<unsigned> neighborSet;

    // First any hulls that intersect ours.
    for (unsigned otherProc = 0; otherProc != numProcs; ++otherProc) {
      if (otherProc != rank and
          convexIntersect(domainHulls[otherProc], domainHulls[rank])) neighborSet.insert(otherProc);
    }

    // Now any hulls that have elements adjacent to ours in the hull mesh.
    // *NOTE*
    // In the original Spheral++ version of this algorithm I used cell adjacency through nodes 
    // as the criterion here, which is more conservative than adjacency through faces.  We don't have
    // node->cell mapping in polytope, however, so for expediency I'm just implementing the face
    // adjacency check for now.  I'm not sure this won't ever miss anything though, so it bears
    // more thought/testing.
    for (unsigned icell = domainCellOffset[rank]; icell != domainCellOffset[rank + 1]; ++icell) {
      for (typename vector<int>::const_iterator faceItr = hullMesh.cells[icell].begin();
           faceItr != hullMesh.cells[icell].end(); ++faceItr) {
        const unsigned iface = (*faceItr >= 0 ? *faceItr : ~(*faceItr));
        ASSERT(iface < hullMesh.faceCells.size());
        for (vector<unsigned>::const_iterator otherCellItr = hullMesh.faceCells[iface].begin();
             otherCellItr != hullMesh.faceCells[iface].end(); ++otherCellItr) {
          const unsigned jcell = *otherCellItr;
          if (jcell != icell) {
            const unsigned otherProc = bisectSearch(domainCellOffset, jcell);
            ASSERT(jcell >= domainCellOffset[otherProc] and
                   jcell <  domainCellOffset[otherProc + 1]);
            neighborSet.insert(otherProc);
          }
        }
      }
    }

    // Copy the neighbor set information the mesh.
    mesh.neighborDomains = vector<unsigned>();
    copy(neighborSet.begin(), neighborSet.end(), back_inserter(mesh.neighborDomains));
    sort(mesh.neighborDomains.begin(), mesh.neighborDomains.end());

#ifndef NDEBUG
    // Make sure everyone is consistent about who talks to whom.
    // This is potentially an expensive check on massively parallel systems!
    {
      for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
        unsigned numOthers = mesh.neighborDomains.size();
        vector<unsigned> otherNeighbors(mesh.neighborDomains);
        MPI_Bcast(&numOthers, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
        if (numOthers > 0) {
          otherNeighbors.resize(numOthers);
          MPI_Bcast(&(otherNeighbors.front()), numOthers, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
          ASSERT(rank == sendProc or
                 count(mesh.neighborDomains.begin(), mesh.neighborDomains.end(), sendProc) == 
                 count(otherNeighbors.begin(), otherNeighbors.end(), rank));
        }
      }
    }
#endif

    // For now we're not going to clever, just send all our generators to our 
    // potential neighbor set.  Pack 'em up.
    vector<char> localBuffer;
    serialize(generators, localBuffer);
    unsigned localBufferSize = localBuffer.size();

    // Fire off our sends asynchronously.
    vector<MPI_Request> sendRequests;
    sendRequests.reserve(2*neighborSet.size());
    for (set<unsigned>::const_iterator otherItr = neighborSet.begin();
         otherItr != neighborSet.end();
         ++otherItr) {
      const unsigned otherProc = *otherItr;
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&localBufferSize, 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &sendRequests.back());
      if (localBufferSize > 0) {
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&localBuffer.front(), localBufferSize, MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &sendRequests.back());
      }
    }
    ASSERT(sendRequests.size() <= 2*neighborSet.size());

    // Get the info from each of our neighbors and append it to the result.
    for (set<unsigned>::const_iterator otherItr = neighborSet.begin();
         otherItr != neighborSet.end();
         ++otherItr) {
      const unsigned otherProc = *otherItr;
      unsigned bufSize;
      MPI_Status recvStatus1, recvStatus2;
      MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &recvStatus1);
      if (bufSize > 0) {
        vector<char> buffer(bufSize, '\0');
        MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &recvStatus2);
        vector<char>::const_iterator itr = buffer.begin();
        vector<RealType> otherGenerators;
        deserialize(otherGenerators, itr, buffer.end());
        ASSERT(itr == buffer.end());
        copy(otherGenerators.begin(), otherGenerators.end(), back_inserter(generators));
      }
    }

    // Make sure all our sends are completed.
    vector<MPI_Status> sendStatus(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &sendRequests.front(), &sendStatus.front());
  }

  // Denormalize the generator positions.
  const unsigned ntotal = generators.size() / Dimension;
  ASSERT(ntotal >= nlocal);
  for (unsigned i = 0; i != ntotal; ++i) {
    for (unsigned j = 0; j != Dimension; ++j) {
      generators[Dimension*i + j] = generators[Dimension*i + j]/fscale + rlow[j];
    }
  }

  // Construct the tessellation including the other domains' generators.
  this->tessellationWrapper(generators, mesh);

  // Now we need to remove the elements of the tessellation corresponding to
  // other domains generators, and renumber the resulting elements.
  vector<unsigned> cellMask(mesh.cells.size(), 1);
  fill(cellMask.begin() + nlocal, cellMask.end(), 0);
  cerr << "Removing " << (mesh.cells.size() - nlocal) << " cells." << endl;
  deleteCells(mesh, cellMask);

  // If requested, build the communication info for the shared nodes & faces
  // with our neighbor domains.
  if (mBuildCommunicationInfo and numProcs > 1) {

    // Hash the positions of our nodes and faces.
    const unsigned numNodes = mesh.nodes.size()/Dimension;
    const unsigned numFaces = mesh.faces.size();
    vector<Point> localNodeHashes, localFaceHashes;
    localNodeHashes.reserve(numNodes);
    localFaceHashes.reserve(numFaces);
    for (unsigned i = 0; i != numNodes; ++i) {
      localNodeHashes.push_back(DimensionTraits<Dimension, RealType>::constructPoint(&mesh.nodes[Dimension*i], rlow, degeneracy, i));
    }
    ASSERT(localNodeHashes.size() == mesh.nodes.size()/Dimension);
    for (unsigned i = 0; i != numFaces; ++i) {
      Point rf;
      rf.index = i;
      for (typename vector<unsigned>::const_iterator itr = mesh.faces[i].begin();
           itr != mesh.faces[i].end(); 
           ++itr) rf += localNodeHashes[*itr];
      rf /= mesh.faces[i].size();
      localFaceHashes.push_back(rf);
    }
    ASSERT(localFaceHashes.size() == numFaces);

    // Serialize our node and face hashes into one spicy meatball.
    vector<char> localBuffer;
    serialize(localNodeHashes, localBuffer);
    serialize(localFaceHashes, localBuffer);
    unsigned localBufferSize = localBuffer.size();

    // Send our node and face hashes to all our neighbors.
    vector<MPI_Request> sendRequests;
    sendRequests.reserve(2*mesh.neighborDomains.size());
    for (vector<unsigned>::const_iterator otherItr = mesh.neighborDomains.begin();
         otherItr != mesh.neighborDomains.end();
         ++otherItr) {
      const unsigned otherProc = *otherItr;
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&localBufferSize, 1, MPI_UNSIGNED, otherProc, 10, MPI_COMM_WORLD, &sendRequests.back());
      if (localBufferSize > 0) {
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&localBuffer.front(), localBufferSize, MPI_CHAR, otherProc, 11, MPI_COMM_WORLD, &sendRequests.back());
      }
    }
    ASSERT(sendRequests.size() <= 2*mesh.neighborDomains.size());

    // Construct sets of our nodes and faces for fast testing.
    const set<Point> myNodes(localNodeHashes.begin(), localNodeHashes.end()),
                     myFaces(localFaceHashes.begin(), localFaceHashes.end());

    // Get the node and face hashes from each of our neighbors.
    vector<pair<unsigned, unsigned> > indices;
    for (vector<unsigned>::const_iterator otherItr = mesh.neighborDomains.begin();
         otherItr != mesh.neighborDomains.end();
         ++otherItr) {
      const unsigned otherProc = *otherItr;
      unsigned bufSize;
      MPI_Status recvStatus1, recvStatus2;
      MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 10, MPI_COMM_WORLD, &recvStatus1);
      if (bufSize > 0) {
        vector<char> buffer(bufSize, '\0');
        MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 11, MPI_COMM_WORLD, &recvStatus2);
        vector<Point> otherNodeHashes, otherFaceHashes;
        vector<char>::const_iterator itr = buffer.begin();
        deserialize(otherNodeHashes, itr, buffer.end());
        deserialize(otherFaceHashes, itr, buffer.end());
        ASSERT(itr == buffer.end());
        ASSERT(otherNodeHashes.size() % Dimension == 0);
        ASSERT(otherFaceHashes.size() % Dimension == 0);

        // Look for any nodes we share with the other domain.
        const unsigned numOtherNodes = otherNodeHashes.size();
        indices = vector<pair<unsigned, unsigned> >();
        for (unsigned i = 0; i != numOtherNodes; ++i) {
          const typename set<Point>::const_iterator itr = myNodes.find(otherNodeHashes[i]);
          if (itr != myNodes.end()) indices.push_back(make_pair(otherNodeHashes[i].index, itr->index));
        }

        // Sort by the order of the minimum domain.
        if (otherProc < rank) {
          sort(indices.begin(), indices.end(), ComparePairByFirstElement<unsigned, unsigned>());
        } else {
          sort(indices.begin(), indices.end(), ComparePairBySecondElement<unsigned, unsigned>());
        }

        // Extract to the node data to the mesh.
        mesh.sharedNodes.push_back(vector<unsigned>());
        mesh.sharedNodes.reserve(indices.size());
        for (unsigned i = 0; i != indices.size(); ++i) mesh.sharedNodes.back().push_back(indices[i].second);
        ASSERT(mesh.sharedNodes.back().size() == indices.size());

        // Look for any faces we share with other domain.
        const unsigned numOtherFaces = otherFaceHashes.size();
        indices = vector<pair<unsigned, unsigned> >();
        for (unsigned i = 0; i != numOtherFaces; ++i) {
          const typename set<Point>::const_iterator itr = myFaces.find(otherFaceHashes[i]);
          if (itr != myFaces.end()) indices.push_back(make_pair(otherFaceHashes[i].index, itr->index));
        }

        // Sort by the order of the minimum domain.
        if (otherProc < rank) {
          sort(indices.begin(), indices.end(), ComparePairByFirstElement<unsigned, unsigned>());
        } else {
          sort(indices.begin(), indices.end(), ComparePairBySecondElement<unsigned, unsigned>());
        }

        // Extract to the face data to the mesh.
        mesh.sharedFaces.push_back(vector<unsigned>());
        mesh.sharedFaces.reserve(indices.size());
        for (unsigned i = 0; i != indices.size(); ++i) mesh.sharedFaces.back().push_back(indices[i].second);
        ASSERT(mesh.sharedFaces.back().size() == indices.size());
      }
    }
  }
}

//------------------------------------------------------------------------------
// Dispatch a serial tessellation call appropriately.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
tessellationWrapper(const vector<RealType>& points,
                    Tessellation<Dimension, RealType>& mesh) const {
  switch (mType) {
  case unbounded:
    mSerialTessellator->tessellate(points, mesh);
    break;

  case box:
    ASSERT(mLow != 0);
    ASSERT(mHigh != 0);
    mSerialTessellator->tessellate(points, mLow, mHigh, mesh);
    break;

  case plc:
    ASSERT(mPLCptr != 0);
    mSerialTessellator->tessellate(points, *mPLCptr, mesh);
    break;
  }
}

//------------------------------------------------------------------------------
// Compute the bouding box to encompass all the points.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
computeBoundingBox(const vector<RealType>& points,
                   RealType* rlow,
                   RealType* rhigh) const {

  const unsigned n = points.size() / Dimension;
  for (unsigned i = 0; i != Dimension; ++i) {
    rlow[i] = numeric_limits<RealType>::max();
    rhigh[i] = (numeric_limits<RealType>::is_signed ? -rlow[i] : numeric_limits<RealType>::min());
  }

  // Find the local min & max.
  switch (mType) {
  case unbounded:
    for (unsigned i = 0; i != n; ++i) {
      for (unsigned j = 0; j != Dimension; ++j) {
        rlow[j] =  min(rlow[j],  points[Dimension*i + j]);
        rhigh[j] = max(rhigh[j], points[Dimension*i + j]);
      }
    }
    break;

  case box:
    ASSERT(mLow != 0);
    ASSERT(mHigh != 0);
    for (unsigned j = 0; j != Dimension; ++j) {
      rlow[j] = mLow[j];
      rhigh[j] = mHigh[j];
    }
    return;

  case plc:
    ASSERT(mPLCptr != 0);
    for (vector<vector<int> >::const_iterator facetItr = mPLCptr->facets.begin();
         facetItr != mPLCptr->facets.end();
         ++facetItr) {
      for (vector<int>::const_iterator iItr = facetItr->begin();
           iItr != facetItr->end();
           ++iItr) {
        const unsigned i = *iItr;
        for (unsigned j = 0; j != Dimension; ++j) {
          rlow[j] =  min(rlow[j],  points[Dimension*i + j]);
          rhigh[j] = max(rhigh[j], points[Dimension*i + j]);
        }
      }
    }
    break;
  }

  // Find the global results.
  for (unsigned j = 0; j != Dimension; ++j) {
    RealType tmp = rlow[j];
    MPI_Allreduce(&tmp, &rlow[j], 1, DataTypeTraits<RealType>::MpiDataType(), MPI_MIN, MPI_COMM_WORLD);
    tmp = rhigh[j];
    MPI_Allreduce(&tmp, &rhigh[j], 1, DataTypeTraits<RealType>::MpiDataType(), MPI_MAX, MPI_COMM_WORLD);
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class DistributedTessellator<2, double>;
template class DistributedTessellator<3, double>;

}
