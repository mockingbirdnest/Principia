#include "quantities/concepts.hpp"

#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {

using namespace principia::quantities::_concepts;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

TEST(Traits, IsQuantityV) {
  static_assert(quantity<int>);
  static_assert(quantity<double>);
  static_assert(quantity<Area>);
  static_assert(quantity<Frequency const>);
  // Not sure if the following is what we want, but at least let's nail it in a
  // test.
  static_assert(!quantity<Entropy&>);
  static_assert(!quantity<float const&>);
}

}  // namespace quantities
}  // namespace principia
