#ifndef __PBGWRAPS_POLYTOPE_CXXTYPES__
#define __PBGWRAPS_POLYTOPE_CXXTYPES__

#include <stdexcept>
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <algorithm>

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef std::set<int>                set_of_int;
typedef std::set<unsigned>           set_of_unsigned;
typedef std::set<double>             set_of_double;

typedef std::vector<int>                vector_of_int;
typedef std::vector<unsigned>           vector_of_unsigned;
typedef std::vector<double>             vector_of_double;

typedef std::vector<std::vector<unsigned> >    vector_of_vector_of_unsigned;
typedef std::vector<std::vector<int> >         vector_of_vector_of_int;
typedef std::vector<std::vector<double> >      vector_of_vector_of_double;
typedef std::vector<std::set<unsigned> >       vector_of_set_of_unsigned;

typedef std::vector<std::vector<std::vector<int> > >  vector_of_vector_of_vector_of_int;

namespace polytope {

//------------------------------------------------------------------------------
// Helpful trait classes.
//------------------------------------------------------------------------------
template<typename Container>
struct DataTypeTraits {
  static unsigned numElements(const Container& x) { return x.size(); }
};

//------------------------------------------------------------------------------
// Index into a Container.
//------------------------------------------------------------------------------
template<typename Container>
inline
typename Container::value_type
indexContainer(Container& container,
               int index) {
  const int n = DataTypeTraits<Container>::numElements(container);
  index = (index < 0 ? (n + 1 + index) : index);
  if (index < n) {
    try {
      return container.at(index);
    } catch (std::out_of_range) {
      PyErr_SetString(PyExc_IndexError, "Container index out of range");
      return typename Container::value_type();
    }
  } else {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return typename Container::value_type();
  }
}

template<typename Container>
inline
typename Container::value_type*
indexContainerAsPointer(Container& container,
                        int index) {
  const int n = DataTypeTraits<Container>::numElements(container);
  index = (index < 0 ? (n + 1 + index) : index);
  if (index < n) {
    try {
      return &(container.at(index));
    } catch (std::out_of_range) {
      PyErr_SetString(PyExc_IndexError, "Container index out of range");
      return NULL;
    }
  } else {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return NULL;
  }
}

template<typename Container>
inline
typename Container::value_type&
indexContainerAsReference(Container& container,
                          int index) {
  const int n = DataTypeTraits<Container>::numElements(container);
  index = (index < 0 ? (n + 1 + index) : index);
  if (index < n) {
    try {
      return container.at(index);
    } catch (std::out_of_range) {
      PyErr_SetString(PyExc_IndexError, "Container index out of range");
      return typename Container::value_type();
    }
  } else {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return NULL;
  }
}

//------------------------------------------------------------------------------
// Extract slices from a container.
//------------------------------------------------------------------------------
template<typename Container>
inline
Container
sliceContainer(Container& container,
               int index1,
               int index2) {
   const int n = DataTypeTraits<Container>::numElements(container);
   index1 = (index1 < 0 ? (n + 1 + index1) : index1);
   index2 = (index2 < 0 ? (n + 1 + index2) : index2);
   Container result;
   for (int i = index1; i < index2; ++i) result.push_back(indexContainer(container, i));
   return result;
}

//------------------------------------------------------------------------------
// Assign to a postion in a container.
//------------------------------------------------------------------------------
template<typename Container>
inline
int
assignToContainerIndex(Container& container, 
                       int index,
                       const typename Container::value_type& value) {
  const int n = DataTypeTraits<Container>::numElements(container);
  index = (index < 0 ? (n + 1 + index) : index);
  if (index < 0 or index >= container.size()) {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return 1;
  } else {
    container[index] = value;
    return 0;
  }
}

template<typename Container>
inline
int
assignToContainerIndexPtr(Container& container, 
                          int index,
                          const typename Container::value_type& value) {
  const int n = DataTypeTraits<Container>::numElements(container);
  index = (index < 0 ? (n + 1 + index) : index);
  if (index < 0 or index >= container.size()) {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return 1;
  } else {
    *(container[index]) = *(value);
    return 0;
  }
}

//------------------------------------------------------------------------------
// Assign to a slice in a container.
//------------------------------------------------------------------------------
template<typename Container>
inline
int
assignToSlice(Container& container, 
              int index1,
              int index2,
              const Container& values) {
   const int n = DataTypeTraits<Container>::numElements(container);
   index1 = (index1 < 0 ? (n + 1 + index1) : index1);
   index2 = (index2 < 0 ? (n + 1 + index2) : index2);
   const int nv = values.size();
   if (index2 - index1 != nv) {
     PyErr_SetString(PyExc_IndexError, "Container slices different sizes.");
     return -1;
   }
   for (int i = index1; i < index2; ++i) {
      if (assignToContainerIndex(container, i, values[i - index1]) != 0) return -1;
   }
   return 0;
}

//------------------------------------------------------------------------------
// Append a value to container of pointers.
//------------------------------------------------------------------------------
template<typename Container>
inline
int
appendToContainerOfPointers(Container& container, 
                            typename Container::value_type value) {
  container.push_back(value);
  return 0;
}

//------------------------------------------------------------------------------
// In-place concatenation.
//------------------------------------------------------------------------------
template<typename Container>
inline
Container&
concatContainersInPlace(Container& lhs,
                        const Container& rhs) {
  const unsigned n = lhs.size() + rhs.size();
  std::copy(rhs.begin(), rhs.end(), std::back_inserter(lhs));
  assert(lhs.size() == n);
  return lhs;
}

//------------------------------------------------------------------------------
// In-place repeat.
//------------------------------------------------------------------------------
template<typename Container>
inline
Container&
repeatContainerInPlace(Container& self,
                       const unsigned count) {
  const unsigned n0 = self.size();
  if (count == 0) {
    self = Container();
  } else {
    self.reserve(n0 * count);
    for (unsigned i = 0; i != count - 1; ++i) {
       std::copy(self.begin(), self.begin() + n0, std::back_inserter(self));
    }
  }
  assert(self.size() == n0 * count);
  return self;
}

//------------------------------------------------------------------------------
// Support concatenation of containers as addition.
//------------------------------------------------------------------------------
template<typename Container>
inline
Container
concatContainers(const Container& lhs,
                 const Container& rhs) {
  Container result(lhs);
  return concatContainersInPlace(result, rhs);
}

//------------------------------------------------------------------------------
// Repeat
//------------------------------------------------------------------------------
template<typename Container>
inline
Container
repeatContainer(Container& self,
                const unsigned count) {
  Container result(self);
  return repeatContainerInPlace(result, count);
}

//------------------------------------------------------------------------------
// Test if a container contains the given value.
//------------------------------------------------------------------------------
template<typename Container>
struct
ContainsValueFunctor {
  inline static int impl(Container& container,
                         const typename Container::value_type& value) {
    return ((std::find(container.begin(), container.end(), value) == container.end()) ?
            0 :
            1);
  }
};

// Specialize for set syntax.
template<typename T>
struct
ContainsValueFunctor<std::set<T> > {
  typedef std::set<T> Container;
  inline static int impl(Container& container,
                         const typename Container::value_type& value) {
    return (container.find(value) == container.end() ?
            0 :
            1);
  }
};

template<typename Container>
inline
int
containsValue(Container& container,
              const typename Container::value_type& value) {
  return ContainsValueFunctor<Container>::impl(container, value);
}

}

#endif
