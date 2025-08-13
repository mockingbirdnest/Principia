#include "quantities/concepts.hpp"

#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace quantities {

using namespace principia::quantities::_concepts;
using namespace principia::quantities::_named_quantities;

TEST(Concepts, Quantity) {
  static_assert(convertible_to_quantity<int>);
  static_assert(convertible_to_quantity<double>);
  static_assert(convertible_to_quantity<Area>);
  static_assert(convertible_to_quantity<Frequency const>);
  static_assert(convertible_to_quantity<Entropy&>);
  static_assert(convertible_to_quantity<float const&>);

  static_assert(!convertible_to_quantity<int*>);
}

}  // namespace quantities
}  // namespace principia
