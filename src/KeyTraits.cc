//-----------------------------------------------------------------------------//
// KeyTraits
//
// Encapsulate how we think about keys for the space filling curves.
// This is specialized assuming we're working with 64 bit unsigned long long.
//----------------------------------------------------------------------------//
#include <limits>
#include "KeyTraits.hh"

namespace polytope {

const uint32_t KeyTraits::numbits = 32U;
const uint32_t KeyTraits::numbits1d = 32U;
const KeyTraits::Key KeyTraits::zero = 0;
const KeyTraits::Key KeyTraits::one = 1;
const KeyTraits::Key KeyTraits::two = 2;
const KeyTraits::Key KeyTraits::minKey1d = std::numeric_limits<Key>::min()/2;
const KeyTraits::Key KeyTraits::maxKey1d = std::numeric_limits<Key>::max()/2;
const KeyTraits::Key KeyTraits::maxKey = std::numeric_limits<Key>::max()/2;

}
