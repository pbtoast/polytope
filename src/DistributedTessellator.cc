//----------------------------------------------------------------------------//
// DistributedTessellator
//----------------------------------------------------------------------------//
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>
#include <list>
#include <limits>
#include "mpi.h"

#include "polytope.hh"
#include "Point.hh"
#include "ReducedPLC.hh"
#include "convexHull_2d.hh"
#include "convexHull_3d.hh"
#include "convexIntersect.hh"
#include "deleteCells.hh"
#include "bisectSearch.hh"
#include "polytope_serialize.hh"
#include "mortonOrderIndices.hh"

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
  typedef polytope::KeyTraits::Key CoordHash;
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
    return Point(ri[0], ri[1], 
                 rlow[0], rlow[1], 
                 dx, i);
  }
  static Point faceCentroid(const polytope::Tessellation<2, RealType>& mesh,
                            const unsigned iface,
                            const RealType* rlow,
                            const RealType& dx) {
    ASSERT(iface < mesh.faces.size());
    ASSERT(mesh.faces[iface].size() == 2);
    RealType pface[2];
    const unsigned n1 = mesh.faces[iface][0], n2 = mesh.faces[iface][1];
    ASSERT(n1 < mesh.nodes.size()/2);
    ASSERT(n2 < mesh.nodes.size()/2);
    pface[0] = 0.5*(mesh.nodes[2*n1]     + mesh.nodes[2*n2]);
    pface[1] = 0.5*(mesh.nodes[2*n1 + 1] + mesh.nodes[2*n2 + 1]);
    return constructPoint(pface, rlow, dx, iface);
  }
  static RealType maxLength(const RealType* low, const RealType* high) {
    return std::max(high[0] - low[0], high[1] - low[1]);
  }
  static vector<RealType> extractCoords(const vector<RealType>& allCoords,
                                        const vector<unsigned>& indices) {
    vector<RealType> result;
    result.reserve(2*indices.size());
    for (vector<unsigned>::const_iterator itr = indices.begin();
         itr != indices.end();
         ++itr) {
      const unsigned i = *itr;
      ASSERT(i < allCoords.size()/2);
      result.push_back(allCoords[2*i]);
      result.push_back(allCoords[2*i + 1]);
    }
    return result;
  }
};

// 3D
template<typename RealType>
struct DimensionTraits<3, RealType> {
  typedef typename polytope::ReducedPLC<3, RealType> ConvexHull;
  typedef polytope::KeyTraits::Key CoordHash;
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
    return Point(ri[0], ri[1], ri[2], 
                 rlow[0], rlow[1], rlow[2],
                 dx, i);
  }
  static Point faceCentroid(const polytope::Tessellation<3, RealType>& mesh,
                            const unsigned iface,
                            const RealType* rlow,
                            const RealType& dx) {
    ASSERT(iface < mesh.faces.size());
    const unsigned nnodes = mesh.faces[iface].size();
    ASSERT(nnodes >= 3);
    RealType pface[3] = {0.0, 0.0, 0.0};
    unsigned ni;
    for (typename vector<unsigned>::const_iterator itr = mesh.faces[iface].begin();
         itr != mesh.faces[iface].end();
         ++itr) {
      ni = mesh.faces[iface][*itr];
      pface[0] += mesh.nodes[3*ni];
      pface[1] += mesh.nodes[3*ni + 1];
      pface[2] += mesh.nodes[3*ni + 2];
    }
    pface[0] /= nnodes;
    pface[1] /= nnodes;
    pface[2] /= nnodes;
    return constructPoint(pface, rlow, dx, iface);
  }
  static RealType maxLength(const RealType* low, const RealType* high) {
    return std::max(std::max(high[0] - low[0], high[1] - low[1]), high[2] - low[2]);
  }
  static vector<RealType> extractCoords(const vector<RealType>& allCoords,
                                        const vector<unsigned>& indices) {
    vector<RealType> result;
    result.reserve(3*indices.size());
    for (vector<unsigned>::const_iterator itr = indices.begin();
         itr != indices.end();
         ++itr) {
      const unsigned i = *itr;
      ASSERT(i < allCoords.size()/3);
      result.push_back(allCoords[3*i]);
      result.push_back(allCoords[3*i + 1]);
      result.push_back(allCoords[3*i + 2]);
    }
    return result;
  }
};

// Traits to handle mapping RealType -> MPI data type.
template<typename RealType> struct DataTypeTraits;

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

//------------------------------------------------------------------------------
// Sort a vector of stuff by the given keys.
//------------------------------------------------------------------------------
template<typename Value, typename Key>
void
sortByKeys(vector<Value>& values, const vector<Key>& keys) {
  ASSERT(values.size() == keys.size());
  vector<pair<Key, Value> > stuff;
  stuff.reserve(values.size());
  for (unsigned i = 0; i != values.size(); ++i) stuff.push_back(make_pair(keys[i], values[i]));
  ASSERT(stuff.size() == values.size());
  sort(stuff.begin(), stuff.end(), ComparePairByFirstElement<Key, Value>());
  for (unsigned i = 0; i != values.size(); ++i) values[i] = stuff[i].second;
}

//------------------------------------------------------------------------------
// Find the set of cells that touch each mesh node.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
vector<set<unsigned> > 
computeNodeCells(const polytope::Tessellation<Dimension, RealType>& mesh) {
  vector<set<unsigned> > result(mesh.nodes.size()/Dimension);
  for (unsigned i = 0; i != mesh.cells.size(); ++i) {
    for (vector<int>::const_iterator faceItr = mesh.cells[i].begin();
         faceItr != mesh.cells[i].end();
         ++faceItr) {
      const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
      for (vector<unsigned>::const_iterator nodeItr = mesh.faces[iface].begin();
           nodeItr != mesh.faces[iface].end();
           ++nodeItr) {
        ASSERT(*nodeItr < result.size());
        result[*nodeItr].insert(i);
      }
    }
  }
  return result;
}

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
  mPLCpointsPtr(0),
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
           const vector<RealType>& PLCpoints,
           const PLC<Dimension, RealType>& geometry,
           Tessellation<Dimension, RealType>& mesh) const {
  mType = plc;
  mPLCpointsPtr = &PLCpoints;
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
  typedef KeyTraits::Key Key;
  const double degeneracy = 1.0e-8;
  
  // Parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Compute the bounding box for normalizing our coordinates.
  RealType rlow[Dimension], rhigh[Dimension], genLow[Dimension];
  this->computeBoundingBox(points, rlow, rhigh);
  RealType fscale = 0;
  for (unsigned i = 0; i != Dimension; ++i) {
    fscale = max(fscale, rhigh[i] - rlow[i]);
    genLow[i] = RealType(0.0);
  }
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
  vector<unsigned> gen2domain(nlocal, rank);
  if (numProcs > 1) {

    // Compute the convex hull of each domain and distribute them to all processes.
    vector<ConvexHull> domainHulls; // (numProcs);
    vector<unsigned> domainCellOffset(1, 0);
    {
      const ConvexHull localHull = DimensionTraits<Dimension, RealType>::convexHull(generators, genLow, degeneracy);
      vector<char> localBuffer;
      serialize(localHull, localBuffer);
      for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
        vector<char> buffer = localBuffer;
        unsigned bufSize = localBuffer.size();
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
        buffer.resize(bufSize);
        MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, MPI_COMM_WORLD);
        vector<char>::const_iterator itr = buffer.begin();
        ConvexHull newHull;
        deserialize(newHull, itr, buffer.end());
        ASSERT(itr == buffer.end());
        domainHulls.push_back(newHull);
        domainCellOffset.push_back(domainCellOffset.back() + domainHulls[sendProc].points.size()/Dimension);
      }
    }
    ASSERT(domainHulls.size() == numProcs);
    ASSERT(domainCellOffset.size() == numProcs + 1);

    // cerr << "Domain cell offsets : ";
    // copy(domainCellOffset.begin(), domainCellOffset.end(), ostream_iterator<unsigned>(cerr, " "));
    // cerr << endl;
    // cerr << " mLow, mHigh : (" << mLow[0] << " " << mLow[1] << ") (" << mHigh[0] << " " << mHigh[1] << ")" << endl;

    // Create a tessellation of the hull vertices for all domains.
    vector<RealType> hullGenerators;
    for (unsigned i = 0; i != numProcs; ++i) {
      copy(domainHulls[i].points.begin(), domainHulls[i].points.end(), back_inserter(hullGenerators));
    }
    ASSERT(hullGenerators.size()/Dimension == domainCellOffset.back());
    Tessellation<Dimension, RealType> hullMesh;
    this->tessellationWrapper(hullGenerators, hullMesh);

    // // Blago!
    // vector<double> r2(hullMesh.cells.size(), rank);
    // map<string, double*> fields;
    // fields["domain"] = &r2[0];
    // cerr << "Writing hull mesh with " << hullMesh.cells.size() << endl;
    // polytope::SiloWriter<Dimension, RealType>::write(hullMesh, fields, "test_DistributedTessellator_hullMesh");
    // MPI_Barrier(MPI_COMM_WORLD);
    // // Blago!

    // Find the set of domains we need to communicate with according to two criteria:
    // 1.  Any domain hull that intersects our own.
    // 2.  Any domain hull that has elements adjacent to one of ours in the hullMesh.
    set<unsigned> neighborSet;

    // First any hulls that intersect ours.
    for (unsigned otherProc = 0; otherProc != numProcs; ++otherProc) {
      if (otherProc != rank and
          convexIntersect(domainHulls[otherProc], domainHulls[rank])) neighborSet.insert(otherProc);
    }

    // We need the set of cells that share nodes.
    const vector<set<unsigned> > hullNodeCells = computeNodeCells(hullMesh);

    // Now any hulls that have elements adjacent to ours in the hull mesh.
    for (unsigned icell = domainCellOffset[rank]; icell != domainCellOffset[rank + 1]; ++icell) {
      for (typename vector<int>::const_iterator faceItr = hullMesh.cells[icell].begin();
           faceItr != hullMesh.cells[icell].end(); ++faceItr) {
        const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
        ASSERT(iface < hullMesh.faceCells.size());
        for (vector<unsigned>::const_iterator nodeItr = hullMesh.faces[iface].begin();
             nodeItr != hullMesh.faces[iface].end();
             ++nodeItr) {
          const unsigned inode = *nodeItr;
          ASSERT(inode < hullNodeCells.size());
          for (set<unsigned>::const_iterator itr = hullNodeCells[inode].begin();
               itr != hullNodeCells[inode].end();
               ++itr) {
            const unsigned otherProc = bisectSearch(domainCellOffset, *itr);
            if (otherProc != rank) neighborSet.insert(otherProc);
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
      // cerr << "Communicating with : ";
      // copy(mesh.neighborDomains.begin(), mesh.neighborDomains.end(), ostream_iterator<unsigned>(cerr, " "));
      // cerr << endl;
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
    sendRequests.reserve(2*mesh.neighborDomains.size());
    for (vector<unsigned>::const_iterator otherItr = mesh.neighborDomains.begin();
         otherItr != mesh.neighborDomains.end();
         ++otherItr) {
      const unsigned otherProc = *otherItr;
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&localBufferSize, 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &sendRequests.back());
      if (localBufferSize > 0) {
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&localBuffer.front(), localBufferSize, MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &sendRequests.back());
      }
    }
    ASSERT(sendRequests.size() <= 2*mesh.neighborDomains.size());

    // Get the info from each of our neighbors and append it to the result.
    // Simultaneously build the mapping of local generator ID to domain.
    for (vector<unsigned>::const_iterator otherItr = mesh.neighborDomains.begin();
         otherItr != mesh.neighborDomains.end();
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
        gen2domain.resize(generators.size()/Dimension, otherProc);
      }
    }

    // Make sure all our sends are completed.
    // cerr << rank << " : Waiting for generator sends to complete." << endl;
    vector<MPI_Status> sendStatus(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &sendRequests.front(), &sendStatus.front());
    // cerr << rank << " : DONE." << endl;
  }
  ASSERT(gen2domain.size() == generators.size()/Dimension);

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

  // // Blago!
  // ASSERT(gen2domain.size() == mesh.cells.size());
  // vector<double> r2(mesh.cells.size());
  // for (unsigned i = 0; i != mesh.cells.size(); ++i) r2[i] = gen2domain[i];
  // map<string, double*> nodeFields, edgeFields, faceFields, cellFields;
  // cellFields["domain"] = &r2[0];
  // cerr << "Writing mesh with " << mesh.cells.size() << endl;
  // polytope::SiloWriter<Dimension, RealType>::write(mesh, nodeFields, edgeFields, faceFields, cellFields, "test_DistributedTessellator_fullMesh");
  // MPI_Barrier(MPI_COMM_WORLD);
  // // Blago!

  // If requested, build the communication info for the shared nodes & faces
  // with our neighbor domains.
  if (mBuildCommunicationInfo and numProcs > 1) {

    // Build the reverse lookup from procID to index in the neighborDomain set.
    map<unsigned, unsigned> proc2offset;
    for (unsigned i = 0; i != mesh.neighborDomains.size(); ++i) proc2offset[mesh.neighborDomains[i]] = i;
    ASSERT(proc2offset.size() == mesh.neighborDomains.size());

    // Look for the faces we share with other processors.
    mesh.sharedFaces.resize(mesh.neighborDomains.size());
    for (int icell = 0; icell != nlocal; ++icell) {

      // Look for shared faces.
      const vector<int>& faces = mesh.cells[icell];
      for (vector<int>::const_iterator faceItr = faces.begin();
           faceItr != faces.end();
           ++faceItr) {
        const int iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
        ASSERT(mesh.faceCells[iface].size() == 1 or mesh.faceCells[iface].size() == 2);
        if (mesh.faceCells[iface].size() == 2) {
          const int jcell1 = mesh.faceCells[iface][0] < 0 ? ~mesh.faceCells[iface][0] : mesh.faceCells[iface][0];
          const int jcell2 = mesh.faceCells[iface][1] < 0 ? ~mesh.faceCells[iface][1] : mesh.faceCells[iface][1];
          ASSERT(jcell1 == icell or jcell2 == icell);
          const int jcell = jcell1 == icell ? jcell2 : jcell1;
          ASSERT(jcell < gen2domain.size());
          const unsigned otherProc = gen2domain[jcell];
          if (otherProc != rank) {
            ASSERT(proc2offset.find(otherProc) != proc2offset.end());
            const unsigned joff = proc2offset[otherProc];
            mesh.sharedFaces[joff].push_back(iface);
          }
        }
      }
    }

    // Look for shared nodes.
    mesh.sharedNodes.resize(mesh.neighborDomains.size());
    const vector<set<unsigned> > nodeCells = computeNodeCells(mesh);
    for (unsigned inode = 0; inode != nodeCells.size(); ++inode) {
      const set<unsigned>& cells = nodeCells[inode];
      for (typename set<unsigned>::const_iterator cellItr = cells.begin();
           cellItr != cells.end();
           ++cellItr) {
        ASSERT(*cellItr < gen2domain.size());
        const unsigned otherProc = gen2domain[*cellItr];
        if (otherProc != rank) {
          ASSERT(proc2offset.find(otherProc) != proc2offset.end());
          const unsigned joff = proc2offset[otherProc];
          mesh.sharedNodes[joff].push_back(inode);
        }
      }
    }

    // // Blago!
    // for (unsigned procID = 0; procID != numProcs; ++procID) {
    //   if (procID == rank) {
    //     for (unsigned iproc = 0; iproc != mesh.neighborDomains.size(); ++iproc) {
    //       cerr << rank << " <-> " << mesh.neighborDomains[iproc] << " : " << mesh.nodes.size()/Dimension << endl;
    //       cerr << "  Nodes : ";
    //       copy(mesh.sharedNodes[iproc].begin(), mesh.sharedNodes[iproc].end(), ostream_iterator<unsigned>(cerr, " "));
    //       cerr << endl;
    //     }
    //   }
    //   MPI_Barrier(MPI_COMM_WORLD);
    // }
    // // Blago!

    // Remove any duplicate shared nodes.  (Faces should already be unique).
    for (unsigned i = 0; i != mesh.sharedNodes.size(); ++i) {
      sort(mesh.sharedNodes[i].begin(), mesh.sharedNodes[i].end());
      mesh.sharedNodes[i].erase(unique(mesh.sharedNodes[i].begin(), mesh.sharedNodes[i].end()), 
                                mesh.sharedNodes[i].end());
    }

    // Compute the bounding box for the mesh coordinates.
    this->computeBoundingBox(mesh.nodes, rlow, rhigh);
    const RealType dx = DimensionTraits<Dimension, RealType>::maxLength(rlow, rhigh)/KeyTraits::maxKey1d;

    // Sort the shared elements by Morton ordering.  This should make the ordering consistent
    // on all domains without communication.
    unsigned numNeighbors = mesh.neighborDomains.size();
    for (unsigned idomain = 0; idomain != numNeighbors; ++idomain) {

      // Nodes.
      vector<Point> nodePoints;
      nodePoints.reserve(mesh.sharedNodes[idomain].size());
      for (vector<unsigned>::const_iterator itr = mesh.sharedNodes[idomain].begin();
           itr != mesh.sharedNodes[idomain].end();
           ++itr) nodePoints.push_back(DimensionTraits<Dimension, RealType>::constructPoint(&(mesh.nodes[2 * (*itr)]),
                                                                                            &rlow[0],
                                                                                            dx,
                                                                                            *itr));
      sort(nodePoints.begin(), nodePoints.end());
      for (unsigned i = 0; i != mesh.sharedNodes[idomain].size(); ++i) mesh.sharedNodes[idomain][i] = nodePoints[i].index;

      // vector<Key> nodeKeys = mortonOrderIndices(nodePoints);
      // sortByKeys(mesh.sharedNodes[idomain], nodeKeys);

      // Faces.
      vector<Point> facePoints;
      for (vector<unsigned>::const_iterator itr = mesh.sharedFaces[idomain].begin();
           itr != mesh.sharedFaces[idomain].end();
           ++itr) facePoints.push_back(DimensionTraits<Dimension, RealType>::faceCentroid(mesh,
                                                                                          *itr,
                                                                                          &rlow[0],
                                                                                          dx));
      sort(facePoints.begin(), facePoints.end());
      for (unsigned i = 0; i != mesh.sharedFaces[idomain].size(); ++i) mesh.sharedFaces[idomain][i] = facePoints[i].index;

      // vector<Key> faceKeys = mortonOrderIndices(facePoints);
      // sortByKeys(mesh.sharedFaces[idomain], faceKeys);
    }
  }

  // Remove the elements of the tessellation corresponding to
  // other domains generators, and renumber the resulting elements.
  vector<unsigned> cellMask(mesh.cells.size(), 1);
  fill(cellMask.begin() + nlocal, cellMask.end(), 0);
  deleteCells(mesh, cellMask);

  // Remove any neighbors we don't actually share any info with.
  unsigned numNeighbors = mesh.neighborDomains.size();
  for (int i = numNeighbors - 1; i != -1; --i) {
    if (mesh.sharedNodes[i].size() == 0 and mesh.sharedFaces[i].size() == 0) {
      // cerr << "Removing neighbor " << i << " of " << numNeighbors << endl;
      mesh.neighborDomains.erase(mesh.neighborDomains.begin() + i);
      mesh.sharedNodes.erase(mesh.sharedNodes.begin() + i);
      mesh.sharedFaces.erase(mesh.sharedFaces.begin() + i);
    }
  }
  numNeighbors = mesh.neighborDomains.size();

  // // Blago!
  // for (unsigned procID = 0; procID != numProcs; ++procID) {
  //   if (procID == rank) {
  //     for (unsigned iproc = 0; iproc != mesh.neighborDomains.size(); ++iproc) {
  //       cerr << rank << " <-> " << mesh.neighborDomains[iproc] << " : " << mesh.nodes.size()/Dimension << endl;
  //       cerr << "  Nodes : ";
  //       copy(mesh.sharedNodes[iproc].begin(), mesh.sharedNodes[iproc].end(), ostream_iterator<unsigned>(cerr, " "));
  //       cerr << endl;
  //     }
  //   }
  //   MPI_Barrier(MPI_COMM_WORLD);
  // }
  // // Blago!

  // In parallel we need to make sure the shared nodes are bit perfect the same.
  if (numProcs > 0) {

    // Figure out which domain owns the shared nodes.
    const unsigned nNodes = mesh.nodes.size()/Dimension;
    vector<unsigned> ownNode(nNodes, rank);
    for (unsigned idomain = 0; idomain != numNeighbors; ++idomain) {
      for (vector<unsigned>::const_iterator itr = mesh.sharedNodes[idomain].begin();
           itr != mesh.sharedNodes[idomain].end();
           ++itr) {
        ASSERT(*itr < nNodes);
        ownNode[*itr] = std::min(ownNode[*itr], mesh.neighborDomains[idomain]);
      }
    }

    // // Blago!
    // for (unsigned p = 0; p != numProcs; ++p) {
    //   if (p == rank) {
    //     for (unsigned i = 0; i != mesh.nodes.size()/Dimension; ++i) {
    //       if (abs(mesh.nodes[2*i] - 0.411765) + abs(mesh.nodes[2*i+1] - 0.235294) < 1.0e-5) {
    //         cerr << "Domain " << rank << " thinks node owned by " << ownNode[i] << endl;
    //       }
    //     }
    //   }
    //   MPI_Barrier(MPI_COMM_WORLD);
    // }
    // // Blago!

    // Post the sends for any nodes we own, and note which nodes we expect to receive.
    vector<unsigned> bufSizes(numNeighbors, 0);
    list<vector<double> > sendCoords;
    list<vector<unsigned> > allRecvNodes;
    vector<MPI_Request> sendRequests;
    sendRequests.reserve(2*numNeighbors);
    for (unsigned idomain = 0; idomain != numNeighbors; ++idomain) {
      vector<unsigned> sendNodes;
      allRecvNodes.push_back(vector<unsigned>());
      vector<unsigned>& recvNodes = allRecvNodes.back();
      for (vector<unsigned>::const_iterator itr = mesh.sharedNodes[idomain].begin();
           itr != mesh.sharedNodes[idomain].end();
           ++itr) {
        ASSERT(ownNode[*itr] <= rank);
        if (ownNode[*itr] == rank) {
          sendNodes.push_back(*itr);
        } else if (ownNode[*itr] == mesh.neighborDomains[idomain]) {
          recvNodes.push_back(*itr);
        }
      }

      // Send any nodes we have for this neighbor.
      bufSizes[idomain] = Dimension*sendNodes.size();
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&bufSizes[idomain], 1, MPI_UNSIGNED,
                mesh.neighborDomains[idomain], 9, MPI_COMM_WORLD, &sendRequests.back());
      if (sendNodes.size() > 0) {
        sendCoords.push_back(DimensionTraits<Dimension, RealType>::extractCoords(mesh.nodes, sendNodes));
        vector<RealType>& coords = sendCoords.back();
        ASSERT(coords.size() == Dimension*sendNodes.size());
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&coords.front(), Dimension*sendNodes.size(), DataTypeTraits<RealType>::MpiDataType(),
                  mesh.neighborDomains[idomain], 10, MPI_COMM_WORLD, &sendRequests.back());
      }
    }
    ASSERT(sendRequests.size() <= 2*numNeighbors);

    // Iterate over the neighbors again and look for any receive information.
    list<vector<unsigned> >::const_iterator recvNodesItr = allRecvNodes.begin();
    for (unsigned idomain = 0; idomain != numNeighbors; ++idomain, ++recvNodesItr) {
      ASSERT(recvNodesItr != allRecvNodes.end());
      const vector<unsigned>& recvNodes = *recvNodesItr;

      unsigned otherSize;
      MPI_Status recvStatus;
      MPI_Recv(&otherSize, 1, MPI_UNSIGNED, mesh.neighborDomains[idomain], 9, MPI_COMM_WORLD, &recvStatus);
      // ASSERT2(otherSize == Dimension*recvNodes.size(),
      //         "Bad message size (" << mesh.neighborDomains[idomain] << " " << otherSize << ") (" << rank << " " << Dimension*recvNodes.size() << ")");
      if (recvNodes.size() > 0) {
        vector<RealType> recvCoords(otherSize); // Dimension*recvNodes.size());
        MPI_Status recvStatus;
        MPI_Recv(&recvCoords.front(), Dimension*recvNodes.size(), DataTypeTraits<RealType>::MpiDataType(),
                 mesh.neighborDomains[idomain], 10, MPI_COMM_WORLD, &recvStatus);

        if (otherSize != Dimension*recvNodes.size()) {
          cout << "Bad message size from " << mesh.neighborDomains[idomain] << " being sent to " << rank << endl;
          cout << "Other coordinates : ";
          copy(recvCoords.begin(), recvCoords.end(), ostream_iterator<RealType>(cout, " "));
          cout << endl;
          cout << "My coordinates : ";
          for (unsigned j = 0; j != recvNodes.size(); ++j) {
            const unsigned i = recvNodes[j];
            for (unsigned k = 0; k != Dimension; ++k) cout << mesh.nodes[Dimension*i + k] << " ";
          }
          cout << endl;
          ASSERT(false);
        }

        // Unpack the coordinates to the receive nodes.
        for (unsigned j = 0; j != recvNodes.size(); ++j) {
          const unsigned i = recvNodes[j];
          for (unsigned k = 0; k != Dimension; ++k) mesh.nodes[Dimension*i + k] = recvCoords[Dimension*j + k];
        }
      }
    }

    // Wait until all our send are complete.
    vector<MPI_Status> sendStatus(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &sendRequests.front(), &sendStatus.front());
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
    mSerialTessellator->tessellate(points, *mPLCpointsPtr, *mPLCptr, mesh);
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
