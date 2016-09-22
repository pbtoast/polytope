//------------------------------------------------------------------------------
// removeElements
//
// Removes the specified elements by index from a std::vector.
// This is needed because there doesn't seem to be a good solution for this in
// the STL.  The std::vector::erase method involves N^2 operations when you're
// removing many elements, while the std::remove and std::remove_if do not work
// for removing elements by index/iterator.
//
// This is borrowed from a similar earlier method in Spheral.  This version also
// returns a mapping of old->new indices.
//
// Created by JMO, Mon Dec 14 13:30:19 PST 2009
//------------------------------------------------------------------------------
#ifndef __polytope_removeElements__
#define __polytope_removeElements__

#include <vector>

namespace polytope {
namespace internal {

template<typename Value, typename index_t>
inline
std::vector<int>
removeElements(std::vector<Value>& vec,
               const std::vector<index_t>& elements) {

  const index_t originalSize = vec.size();
  const index_t newSize = originalSize - elements.size();
  std::vector<int> old2new(originalSize, -1);

  // Is there anything to do?
  if (elements.size() > 0) {

    // Pre-conditions.
    POLY_BEGIN_CONTRACT_SCOPE;
    {
      // We require the input IDs be sorted and unique.
      for (typename std::vector<index_t>::const_iterator itr = elements.begin();
           itr + 1 < elements.end();
           ++itr) {
        POLY_ASSERT(*itr < *(itr + 1));
      }
      if (elements.size() > 0) 
        POLY_ASSERT(elements[0] >= 0 && elements.back() < originalSize);
    }
    POLY_END_CONTRACT_SCOPE;

    // Fill in the return indexing up to the first place we start remvoving elements.
    for (index_t i = 0; i != elements[0]; ++i) old2new[i] = i;

    // Remove the elements.
    // We prefer not to use the vector::erase here 'cause if we're removing
    // many elements the copy and move behaviour of erase can make this
    // an N^2 thing.  Yuck!
    typename std::vector<index_t>::const_iterator delItr = elements.begin();
    typename std::vector<index_t>::const_iterator endItr = elements.end();
    index_t i = *delItr;
    index_t j = i + 1;
    ++delItr;
    while (j != originalSize and delItr != endItr) {
      if (j == *delItr) {
        ++delItr;
        ++j;
      } else {
        old2new[j] = i;
        vec[i] = vec[j];
        ++i;
        ++j;
      }
    }
    if (j != originalSize) {
      std::copy(vec.begin() + j, vec.end(), vec.begin() + i);
      for (size_t k = j; k != originalSize; ++k) old2new[k] = i + k - j;
    }

    // Resize vec to it's new size.
    vec.erase(vec.begin() + newSize, vec.end());

  } else {
    for (size_t i = 0; i != originalSize; ++i) old2new[i] = i;
  }

  // Post-conditions.
  POLY_ASSERT(vec.size() == newSize);
  POLY_ASSERT(old2new.size() == originalSize);
  return old2new;
}

}
}

#endif
