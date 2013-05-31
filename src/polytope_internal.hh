// polytope_internal.hh
//
// Put common includes for polytope here that you don't necessarily 
// want exposed in the public interface.
#include <map>

// An POLY_ASSERT macro, if one isn't already defined.
#ifndef POLY_ASSERT

namespace polytope {

// Forward declare our helper abort method.
void internal_abort();

#ifndef NDEBUG
#define POLY_ASSERT(x) \
  if (!(x)) \
  { \
    std::cout << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl; \
    polytope::internal_abort(); \
  }
#define POLY_ASSERT2(x, msg) \
  if (!(x)) \
  { \
    std::cout << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl << msg << std::endl; \
    polytope::internal_abort(); \
  }
#define POLY_BEGIN_CONTRACT_SCOPE if (false) { while(false)
#define POLY_END_CONTRACT_SCOPE } while(false)
#else
#define POLY_ASSERT(x)
#define POLY_ASSERT2(x, msg)
#define POLY_BEGIN_CONTRACT_SCOPE 
#define POLY_END_CONTRACT_SCOPE 
#endif

namespace internal {


//------------------------------------------------------------------------------
// Pair comparator for first index only
//------------------------------------------------------------------------------
template <typename T1, typename T2>
bool
pairCompareFirst( std::pair<T1,T2> x, std::pair<T1,T2> y ) {
  return (x.first < y.first);
}

//------------------------------------------------------------------------------
// Edge comparator
//------------------------------------------------------------------------------
template <typename T1, typename T2>
bool
pairCompare(const std::pair<T1,T2> pair1, const std::pair<T1,T2> pair2) {
   return (pair1.first  <  pair2.first      ? true :
           pair1.first  == pair2.first and
           pair1.second <  pair2.second     ? true : false);
}

//------------------------------------------------------------------------------
// Update the map of thingies to unique indices.
//------------------------------------------------------------------------------
template<typename Key,
	 typename Comparator>
int
addKeyToMap(const Key& key, std::map<Key, int, Comparator>& key2id) {
  const typename std::map<Key, int>::const_iterator itr = key2id.find(key);
  int result;
  if (itr == key2id.end()) {
    result = key2id.size();
    key2id[key] = result;
  } else {
    result = itr->second;
  }
  return result;
}

//------------------------------------------------------------------------------
// An implementation of the map specialized to help constructing counters.
// This thing just overloads the index operator to start the count at zero
// for new key values.
//------------------------------------------------------------------------------
template<typename Key, 
         typename Comparator = std::less<Key> >
class CounterMap: public std::map<Key, unsigned> {
public:
  CounterMap(): std::map<Key, unsigned>() {}
  virtual ~CounterMap() {}
  unsigned& operator[](const Key& key) {
    typename std::map<Key, unsigned>::iterator itr = this->find(key);
    if (itr == this->end()) {
      std::map<Key, unsigned>::operator[](key) = 0U;
      itr = this->find(key);
    }
    POLY_ASSERT(itr != this->end());
    return itr->second;
  }
};

//------------------------------------------------------------------------------
// Hash two node indices uniquely to represent an edge.
//------------------------------------------------------------------------------
inline
std::pair<int, int>
hashEdge(const int i, const int j) {
  POLY_ASSERT(i != j);
  return i < j ? std::make_pair(i, j) : std::make_pair(j, i);
}

//------------------------------------------------------------------------------
// Return the positive index.
//------------------------------------------------------------------------------
inline
int
positiveID(const int x) {
  return x >= 0 ? x : ~x;
}

}

}

// Classes within the library.
#include "ErrorHandler.hh"

#endif
