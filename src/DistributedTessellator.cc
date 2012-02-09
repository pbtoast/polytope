//----------------------------------------------------------------------------//
// DistributedTessellator
//----------------------------------------------------------------------------//
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>
#include <limits>
extern "C" {
#include "mpi.h"
}

#include "polytope.hh"
#include "polytope_serialize.hh"
#include "ReducedPLC.hh"
#include "convexHull_2d.hh"
#include "convexIntersect.hh"
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

  static ConvexHull convexHull(const vector<RealType>& points, 
                               const RealType* low,
                               const RealType& dx) { 
    return ConvexHull(polytope::convexHull_2d(points, low, dx), points);
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

} // end anonymous namespace

namespace polytope {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
DistributedTessellator<Dimension, RealType>::
DistributedTessellator(const Tessellator<Dimension, RealType>& tessellator):
  mSerialTessellator(tessellator),
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

  // Parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Compute the bounding box for normalizing our coordinates.
  RealType rlow[Dimension], rhigh[Dimension];
  computeBoundingBox(points, rlow, rhigh);
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
    const ConvexHull localHull = DimensionTraits<Dimension, RealType>::convexHull(generators, rlow, 1.0e-14);
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
    ASSERT(domainCellOffset.size() == numProcs);

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
          convexIntersect(domainHulls[otherProc], domainHulls[rank],
                          domainHulls[otherProc].points, domainHulls[rank].points)) neighborSet.insert(otherProc);
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
          ASSERT(jcell < numProcs);
          if (jcell != icell) {
            const unsigned otherProc = bisectSearch(domainCellOffset, jcell);
            ASSERT(jcell >= domainCellOffset[otherProc] and
                   jcell <  domainCellOffset[otherProc + 1]);
            neighborSet.insert(otherProc);
          }
        }
      }
    }

#ifndef NDEBUG
    // Make sure everyone is consistent about who talks to whom.
    // This is potentially an expensive check on massively parallel systems!
    {
      vector<unsigned> neighborDomains;
      copy(neighborSet.begin(), neighborSet.end(), back_inserter(neighborDomains));
      sort(neighborDomains.begin(), neighborDomains.end());
      for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
        unsigned numOthers = neighborDomains.size();
        vector<unsigned> otherNeighbors(neighborDomains);
        MPI_Bcast(&numOthers, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
        if (numOthers > 0) {
          otherNeighbors.resize(numOthers);
          MPI_Bcast(&(otherNeighbors.front()), numOthers, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
          ASSERT(rank == sendProc or
                 count(neighborDomains.begin(), neighborDomains.end(), sendProc) == 
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
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&localBuffer.front(), localBufferSize, MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &sendRequests.back());
    }
    ASSERT(sendRequests.size() == 2*neighborSet.size());

    // Get the info from each of our neighbors and append it to the result.
    for (set<unsigned>::const_iterator otherItr = neighborSet.begin();
         otherItr != neighborSet.end();
         ++otherItr) {
      const unsigned otherProc = *otherItr;
      unsigned bufSize;
      MPI_Status recvStatus1, recvStatus2;
      MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &recvStatus1);
      vector<char> buffer(bufSize, '\0');
      MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &recvStatus2);
      if (bufSize > 0) {
        vector<char>::const_iterator itr = buffer.begin();
        deserialize(generators, itr, buffer.end());
        ASSERT(itr == buffer.end());
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
      generators[2*i + j] = generators[2*i + j]/fscale + rlow[j];
    }
  }

  // Construct the tessellation including the other domains' generators.
  this->tessellationWrapper(generators, mesh);

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
    mSerialTessellator.tessellate(points, mesh);
    break;

  case box:
    ASSERT(mLow != 0);
    ASSERT(mHigh != 0);
    mSerialTessellator.tessellate(points, mLow, mHigh, mesh);
    break;

  case plc:
    ASSERT(mPLCptr != 0);
    mSerialTessellator.tessellate(points, *mPLCptr, mesh);
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
    rlow = mLow;
    rhigh = mHigh;
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

}
