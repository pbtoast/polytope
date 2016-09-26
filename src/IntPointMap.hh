//------------------------------------------------------------------------------
// A struct to help collapsing vertices to unique values, allowing a tolerance
// for comparison of coordinate values.
//
// This is necessary 'cause the STL types that take a comparator don't like
// fuzzy comparisons that can result in a==b==c, a < c.
//
// Current implementation assumes you will match to at most one ID for a given
// Point value.
//------------------------------------------------------------------------------
#ifndef __polytope_IntPointMap__
#define __polytope_IntPointMap__

#include <vector>
#include <utility>
#include <algorithm>
#include "DimensionTraits.hh"
#include "polytope_internal.hh"

namespace polytope {
namespace internal {

template<int Dimension, typename IntType>
struct IntPointMap {

  typedef typename DimensionTraits<Dimension, IntType>::IntPoint IntPoint;
  typedef std::vector<std::pair<IntType, std::vector<size_t> > > Container;
  std::vector<Container> mCoords2IDs;
  IntType mtol;
  size_t mCurrentID;

  //----------------------------------------------------------------------------
  // Constructor.
  //----------------------------------------------------------------------------
  IntPointMap(const IntType tol):
    mCoords2IDs(Dimension),
    mtol(tol),
    mCurrentID(0) {}

  //----------------------------------------------------------------------------
  // Find the index corresponding to the given Point.
  //----------------------------------------------------------------------------
  size_t index(const IntPoint& p) {
    // const bool barf = (p == IntPoint(-71582789, -214748366));
    // const bool barf = (p == IntPoint(71582786,-71582789));
    std::vector<typename Container::iterator> lowers(Dimension), uppers(Dimension);
    bool empty = false;
    // if (barf) std::cerr << "----------------------------------------" << std::endl;
    for (unsigned idim = 0; idim != Dimension; ++idim) {
      lowers[idim] = std::lower_bound(mCoords2IDs[idim].begin(), mCoords2IDs[idim].end(),
                                      make_pair(p[idim] - mtol, std::vector<size_t>()),
                                      ComparePairByFirstElement<IntType, std::vector<size_t> >());
      uppers[idim] = std::upper_bound(mCoords2IDs[idim].begin(), mCoords2IDs[idim].end(), 
                                      make_pair(p[idim] + mtol, std::vector<size_t>()),
                                      ComparePairByFirstElement<IntType, std::vector<size_t> >());
      empty |= (lowers[idim] == uppers[idim]);
      // if (barf) std::cerr << " BARF: " << idim << " " << empty << " : " << std::distance(lowers[idim], uppers[idim]) << std::endl;
    }
    // if (barf) std::cerr << "----------------------------------------" << std::endl;

    // Can we already rule this a new point?
    if (empty) {
      for (unsigned idim = 0; idim != Dimension; ++idim) {
        typename Container::iterator itr = std::lower_bound(lowers[idim], uppers[idim], 
                                                            make_pair(p[idim], std::vector<size_t>()),
                                                            ComparePairByFirstElement<IntType, std::vector<size_t> >());
        if (itr == mCoords2IDs[idim].end() or itr->first != p[idim]) {
          mCoords2IDs[idim].insert(itr, std::make_pair(p[idim], std::vector<size_t>(1, mCurrentID)));
        } else {
          itr->second.push_back(mCurrentID);
        }
      }
      // if (barf) std::cerr << " BARF: emtpy: " << mCurrentID << std::endl;
      // if (barf) this->printState();
      return mCurrentID++;
    }

    // Otherwise there are points in each coordinate scan.  We need to check the intersections
    // of these sets.
    // if (barf) std::cerr << " BARF: not empty" << std::endl;
    std::vector<size_t> matchingIDs;
    for (unsigned idim = 0; idim != Dimension; ++idim) {
      std::vector<size_t> dimIDs;
      // if (barf) std::cerr << " ---> " << idim;
      for (typename Container::const_iterator itr = lowers[idim]; itr != uppers[idim]; ++itr) {
        // if (barf) for (typename std::vector<size_t>::const_iterator ii = itr->second.begin(); ii != itr->second.end(); ++ii) std::cerr << " " << (*ii);
        std::copy(itr->second.begin(), itr->second.end(), std::back_inserter(dimIDs));
      }
      std::sort(dimIDs.begin(), dimIDs.end());
      if (idim == 0) {
        matchingIDs = dimIDs;
      } else {
        std::vector<size_t> intersection;
        std::set_intersection(matchingIDs.begin(), matchingIDs.end(),
                              dimIDs.begin(), dimIDs.end(),
                              std::back_inserter(intersection));
        if (intersection.size() == 0) {
          // If the intersection is empty, we have a new distinct point.
          // Start by updating the internal state with the new ID.
          for (unsigned idim = 0; idim != Dimension; ++idim) {
            typename Container::iterator itr = std::lower_bound(lowers[idim], uppers[idim], 
                                                                make_pair(p[idim], std::vector<size_t>()),
                                                                ComparePairByFirstElement<IntType, std::vector<size_t> >());
            if (itr == mCoords2IDs[idim].end() or itr->first != p[idim]) {
              mCoords2IDs[idim].insert(itr, std::make_pair(p[idim], std::vector<size_t>(1, mCurrentID)));
            } else {
              itr->second.push_back(mCurrentID);
            }
          }
          // if (barf) std::cerr << std::endl << " BARF: intersection empty, so returning " << mCurrentID << std::endl;
          // if (barf) this->printState();
          return mCurrentID++;
        }
        matchingIDs = intersection;
      }
      // if (barf) std::cerr << std::endl;
    }

    // Now we get to the assumption that in the end we'll only have one final match.  If there is more than one value matching
    // within our degeneracy this breaks other current assumptions in polytop, so we'll need to stop and do something more clever.
    POLY_ASSERT2(matchingIDs.size() == 1,
                 "More than one degenerate coordinate match at point " << p << " (total " << matchingIDs.size() << ")");
    // if (barf) std::cerr << " BARF: found match " << matchingIDs[0] << std::endl;
    // if (barf) this->printState();
    return matchingIDs[0];
  }

  //----------------------------------------------------------------------------
  // Find the index corresponding to the given Point.
  //----------------------------------------------------------------------------
  bool have(const IntPoint& p) const {
    std::vector<typename Container::const_iterator> lowers(Dimension), uppers(Dimension);
    bool empty = false;
    for (unsigned idim = 0; idim != Dimension; ++idim) {
      lowers[idim] = std::lower_bound(mCoords2IDs[idim].begin(), mCoords2IDs[idim].end(),
                                      make_pair(p[idim] - mtol, std::vector<size_t>()),
                                      ComparePairByFirstElement<IntType, std::vector<size_t> >());
      uppers[idim] = std::upper_bound(mCoords2IDs[idim].begin(), mCoords2IDs[idim].end(), 
                                      make_pair(p[idim] + mtol, std::vector<size_t>()),
                                      ComparePairByFirstElement<IntType, std::vector<size_t> >());
      empty |= (lowers[idim] == uppers[idim]);
    }
    if (empty) return false;

    // There are points in each coordinate scan.  We need to check the intersections
    // of these sets.
    std::vector<size_t> matchingIDs;
    for (unsigned idim = 0; idim != Dimension; ++idim) {
      std::vector<size_t> dimIDs;
      for (typename Container::const_iterator itr = lowers[idim]; itr != uppers[idim]; ++itr) {
        std::copy(itr->second.begin(), itr->second.end(), std::back_inserter(dimIDs));
      }
      std::sort(dimIDs.begin(), dimIDs.end());
      if (idim == 0) {
        matchingIDs = dimIDs;
      } else {
        std::vector<size_t> intersection;
        std::set_intersection(matchingIDs.begin(), matchingIDs.end(),
                              dimIDs.begin(), dimIDs.end(),
                              std::back_inserter(intersection));
        if (intersection.size() == 0) {
          // If the intersection is empty, we have a new distinct point.
          return false;
        }
      }
    }

    // We must have the point.
    return true;
  }

  //----------------------------------------------------------------------------
  // Print the current known keys and indices.
  //----------------------------------------------------------------------------
  void printState() const {
    std::cout << "IntPointMat coords -> IDs:" << std::endl;
    for (unsigned idim = 0; idim != Dimension; ++idim) {
      std::cout << "Dim " << idim << std::endl;
      for (unsigned i = 0; i != mCoords2IDs[idim].size(); ++i) {
        std::cout << "  " << mCoords2IDs[idim][i].first << " :";
        const std::vector<size_t>& ids = mCoords2IDs[idim][i].second;
        for (unsigned j = 0; j != ids.size(); ++j) std::cout << " " << ids[j];
        std::cout << std::endl;
      }
    }
  }
};

}
}

#else

namespace polytope {
namespace internal {
template<int Dimension, typename IntType> struct IntPointMap;
}
}

#endif
