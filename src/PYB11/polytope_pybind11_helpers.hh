//------------------------------------------------------------------------------
// Some useful methods for binding polytope and pybind11.
//------------------------------------------------------------------------------
#include "polytope.hh"

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/functional.h"

#include <vector>
#include <map>

namespace py = pybind11;

namespace polytope {
namespace pybind11_helpers {

template<typename Key, typename Value>
std::map<Key, std::vector<Value>>
copyDictToMap(py::dict& pydict) {
  std::map<Key, std::vector<Value>> result;
  for (const auto item: pydict) {
    const auto& key = item.first.cast<Key>();
    const auto& vallist = item.second.cast<py::list>();
    result[key] = vector<Value>();
    for (const auto vali: vallist) result[key].push_back(vali.cast<Value>());
  }
  return result;
}

}
}
