//----------------------------------------------------------------------------//
// A bunch of utility methods to help with serializing/unserializing objects
// in polytope.
//----------------------------------------------------------------------------//
#ifndef __polytope_serialize__
#define __polytope_serialize__

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <stdint.h>

namespace polytope {

//------------------------------------------------------------------------------
// Serialize.  Due to limitations on partial specialization for functions, we 
// define a Functor object to handle specializing for each data type.
//------------------------------------------------------------------------------
// Generic primitive types.
template<typename T>
struct Serializer {

  static void serializeImpl(const T& value, 
                            std::vector<char>& buffer) {
    const unsigned n = sizeof(T);
    const char* data = reinterpret_cast<const char*>(&value);
    std::copy(data, data + n, std::back_inserter(buffer));
  }

  static void deserializeImpl(T& value, 
                              std::vector<char>::const_iterator& bufItr, 
                              const std::vector<char>::const_iterator& endItr) {
    const unsigned n = sizeof(T);
    char* data = reinterpret_cast<char*>(&value);
    std::copy(bufItr, bufItr + n, data);
    bufItr += n;
    ASSERT(bufItr <= endItr);
  }
};

// std::vector of known type.
template<typename T>
struct Serializer<std::vector<T> > {

  static void serializeImpl(const std::vector<T>& value, 
                            std::vector<char>& buffer) {
    const unsigned n = value.size();
    Serializer<unsigned>::serializeImpl(n, buffer);
    for (unsigned i = 0; i != n; ++i) Serializer<T>::serializeImpl(value[i], buffer);
  }

  static void deserializeImpl(std::vector<T>& value, 
                              std::vector<char>::const_iterator& bufItr, 
                              const std::vector<char>::const_iterator& endItr) {
    unsigned n;
    Serializer<unsigned>::deserializeImpl(n, bufItr, endItr);
    value.resize(n);
    for (unsigned i = 0; i != n; ++i) Serializer<T>::deserializeImpl(value[i], bufItr, endItr);
  }
};

//------------------------------------------------------------------------------
// Dispatch to the functors!
//------------------------------------------------------------------------------
template<typename T>
void serialize(const T& value, 
               std::vector<char>& buffer) {
  Serializer<T>::serializeImpl(value, buffer);
}

template<typename T>
void deserialize(T& value,
                 std::vector<char>::const_iterator& bufItr,
                 const std::vector<char>::const_iterator& endItr) {
  Serializer<T>::deserializeImpl(value, bufItr, endItr);
}

}

#endif
