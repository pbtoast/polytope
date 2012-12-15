//----------------------------------------------------------------------------//
// SerialDistributedTessellator
// 
// Provides a parallel tessellation.
//
// This class assumes the user provides:
// 1.  A serial Tessellator.
// 2.  The generators in parallel, distributed in any manner the user likes
//     so long as the generators are not degenerate: i.e., don't repeast the
//     same generator on different domains.
//
// This is a dumb but reliable version of DistributedTessellator.  This version
// exchanges all generators to everyone, builds the full tessellation on each
// processor, and then just keeps the local portion needed.  Obviously not 
// scalable, but useful as a check for the real DistributedTessellator.
//----------------------------------------------------------------------------//
#ifndef __Polytope_SerialDistributedTessellator__
#define __Polytope_SerialDistributedTessellator__

#include "DistributedTessellator.hh"

namespace polytope {

template<int Dimension, typename RealType>
class SerialDistributedTessellator: public DistributedTessellator<Dimension, RealType> {

  //--------------------------- Public Interface ---------------------------//
public:

  //! Constructor.
  //! \param serialTessellator A serial implementation of Tessellator.
  //! \param assumeControl If set to true, the DistributedTessellator will 
  //!                      assume control of the serial tessellator.
  SerialDistributedTessellator(Tessellator<Dimension, RealType>* serialTessellator,
                               bool assumeControl = true,
                               bool buildCommunicationInfo = false);
  virtual ~SerialDistributedTessellator();

protected:
  // Internal method that implements the parallel algorithm -- called by
  // each of the three ways we support doing tessellations.
  virtual void computeDistributedTessellation(const std::vector<RealType>& points,
                                              Tessellation<Dimension, RealType>& mesh) const;

private:
  // Forbidden methods.
  SerialDistributedTessellator();
  SerialDistributedTessellator(const SerialDistributedTessellator&);
  SerialDistributedTessellator& operator=(const SerialDistributedTessellator&);
};

}

#endif
