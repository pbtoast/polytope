#ifndef POLYTOPE_PLC_SERIALIZE_HH
#define POLYTOPE_PLC_SERIALIZE_HH

#include <vector>

#include "PLC.hh"
#include "ReducedPLC.hh"
#include "polytope_serialize.hh"

namespace polytope
{

//------------------------------------------------------------------------------
// Serialize a PLC.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
struct Serializer<PLC<Dimension, RealType> >
{

  static void serializeImpl(const PLC<Dimension, RealType>& value,
                            std::vector<char>& buffer) {
    serialize(value.facets, buffer);
    serialize(value.holes, buffer);
  }

  static void deserializeImpl(PLC<Dimension, RealType>& value,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator& endItr) {
    deserialize(value.facets, bufItr, endItr);
    deserialize(value.holes, bufItr, endItr);
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
