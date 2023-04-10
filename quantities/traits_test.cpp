#include "quantities/traits.hpp"

#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace quantities {

using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_traits;

TEST(Traits, IsQuantityV) {
  static_assert(is_quantity_v<int>);
  static_assert(is_quantity_v<double>);
  static_assert(is_quantity_v<Area>);
  static_assert(is_quantity_v<Frequency const>);
  // Not sure if the following is what we want, but at least let's nail it in a
  // test.
  static_assert(!is_quantity_v<Entropy&>);
  static_assert(!is_quantity_v<float const&>);
}

}  // namespace quantities
}  // namespace principia
