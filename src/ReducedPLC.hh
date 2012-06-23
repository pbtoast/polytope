#ifndef POLYTOPE_REDUCEDPLC_HH
#define POLYTOPE_REDUCEDPLC_HH

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include "polytope_serialize.hh"
#include "PLC.hh"

namespace polytope {

//! \class ReducedPLC - A Piecewise Linear Complex in 3D, or a Planar Straight Line 
//! Graph (PSLG) in 2D.
//! The reduced PLC is a PLC which
//!  a) contains its own generating points,
//!  b) is reduced to just the generating points that are used in the PLC.
template<int Dimension, typename RealType>
class ReducedPLC: public PLC<Dimension, RealType> {
public:

  //! This array of size (Dimension*numPoints) contains components of 
  //! the points that are used in the facets of the PLC.
  std::vector<RealType> points;

  //! Default constructor.
  ReducedPLC(): 
    PLC<Dimension, RealType>(),
    points() {}

  //! Construct from a normal PLC, copying the necessary generator coordinates
  //! to our internal data.
  ReducedPLC(const PLC<Dimension, RealType>& plc,
             const std::vector<RealType>& allpoints): 
    PLC<Dimension, RealType>(plc),
    points()
  {
    std::set<int> indices;
    int i;
    for (i = 0; i != plc.facets.size(); ++i) {
      std::copy(plc.facets[i].begin(), plc.facets[i].end(), std::inserter(indices, indices.end()));
    }
    std::map<int, int> old2new;
    int j = 0;
    for (typename std::set<int>::const_iterator itr = indices.begin();
         itr != indices.end();
         ++itr ) {
      i = *itr;
      old2new[i] = j++;
      std::copy(allpoints.begin() + Dimension*i, allpoints.begin() + Dimension*(i + 1), std::back_inserter(points));
    }
    for (i = 0; i != this->facets.size(); ++i) {
      for (j = 0; j != this->facets[i].size(); ++j) {
        this->facets[i][j] = old2new[this->facets[i][j]];
      }
    }
  }

};

//------------------------------------------------------------------------------
// Serialize a ReducedPLC.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
struct Serializer<ReducedPLC<Dimension, RealType> > {

  static void serializeImpl(const ReducedPLC<Dimension, RealType>& value,
                            std::vector<char>& buffer) {
    const unsigned nf = value.facets.size();
    serialize(nf, buffer);
    for (unsigned i = 0; i != nf; ++i) serialize(value.facets[i], buffer);
    serialize(value.holes, buffer);
    serialize(value.points, buffer);
  }

  static void deserializeImpl(ReducedPLC<Dimension, RealType>& value,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator& endItr) {
    unsigned nf;
    deserialize(nf, bufItr, endItr);
    value.facets.resize(nf);
    for (unsigned i = 0; i != nf; ++i) deserialize(value.facets[i], bufItr, endItr);
    deserialize(value.holes, bufItr, endItr);
    deserialize(value.points, bufItr, endItr);
  }
};

}

#endif
