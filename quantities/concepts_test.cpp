#include "quantities/concepts.hpp"

#include <string>

#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {

using namespace principia::quantities::_concepts;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

TEST(Concepts, IsQuantityV) {
  static_assert(convertible_to_quantity<int>);
  static_assert(convertible_to_quantity<double>);
  static_assert(convertible_to_quantity<Area>);
  static_assert(convertible_to_quantity<Frequency const>);
  static_assert(convertible_to_quantity<Entropy&>);
  static_assert(convertible_to_quantity<float const&>);

  static_assert(!convertible_to_quantity<int*>);
}

TEST(Concepts, Algebra) {
  static_assert(!affine<std::string>);
  static_assert(affine<std::byte*>);
  static_assert(affine<Length>);
  static_assert(affine<int>);
  static_assert(affine<double>);

  static_assert(!additive_group<std::string>);
  static_assert(!additive_group<std::byte*>);
  static_assert(additive_group<Length>);
  static_assert(additive_group<int>);
  static_assert(additive_group<double>);

  static_assert(homogeneous_ring<Length>);
  static_assert(homogeneous_ring<int>);
  static_assert(homogeneous_ring<double>);

  static_assert(!ring<Length>);
  static_assert(ring<int>);
  static_assert(ring<double>);

  static_assert(homogeneous_field<Length>);
  static_assert(!homogeneous_field<int>);
  static_assert(homogeneous_field<double>);

  static_assert(!field<Length>);
  static_assert(!field<int>);
  static_assert(field<double>);
}

TEST(Concepts, LinearAlgebra) {
  static_assert(homogeneous_vector_space<Length, Length>);
  static_assert(!vector_space<Length, Length>);
  static_assert(vector_space<double, double>);

  static_assert(!real_vector_space<int>);
}

}  // namespace quantities
}  // namespace principia
