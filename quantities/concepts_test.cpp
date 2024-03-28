#include "quantities/concepts.hpp"

#include <string>

#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {

using namespace principia::numerics::_fixed_arrays;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

TEST(Traits, IsQuantityV) {
  static_assert(convertible_to_quantity<int>);
  static_assert(convertible_to_quantity<double>);
  static_assert(convertible_to_quantity<Area>);
  static_assert(convertible_to_quantity<Frequency const>);
  static_assert(convertible_to_quantity<Entropy&>);
  static_assert(convertible_to_quantity<float const&>);

  static_assert(!convertible_to_quantity<int*>);
}

TEST(Traits, Algebra) {
  static_assert(!affine<std::string>);
  static_assert(affine<std::byte*>);
  static_assert(affine<FixedMatrix<double, 2, 3>>);
  static_assert(affine<FixedMatrix<double, 3, 3>>);
  static_assert(affine<FixedMatrix<Length, 3, 3>>);
  static_assert(affine<Length>);
  static_assert(affine<int>);
  static_assert(affine<double>);

  static_assert(!additive_group<std::string>);
  static_assert(!additive_group<std::byte*>);
  static_assert(additive_group<FixedMatrix<double, 2, 3>>);
  static_assert(additive_group<FixedMatrix<double, 3, 3>>);
  static_assert(additive_group<FixedMatrix<Length, 3, 3>>);
  static_assert(additive_group<Length>);
  static_assert(additive_group<int>);
  static_assert(additive_group<double>);

  static_assert(!homogeneous_ring<FixedMatrix<double, 2, 3>>);
  static_assert(homogeneous_ring<FixedMatrix<double, 3, 3>>);
  static_assert(homogeneous_ring<FixedMatrix<Length, 3, 3>>);
  static_assert(homogeneous_ring<Length>);
  static_assert(homogeneous_ring<int>);
  static_assert(homogeneous_ring<double>);

  static_assert(ring<FixedMatrix<double, 3, 3>>);
  static_assert(!ring<FixedMatrix<Length, 3, 3>>);
  static_assert(!ring<Length>);
  static_assert(ring<int>);
  static_assert(ring<double>);

  static_assert(!homogeneous_field<FixedMatrix<double, 3, 3>>);
  static_assert(!homogeneous_field<FixedMatrix<Length, 3, 3>>);
  static_assert(homogeneous_field<Length>);
  static_assert(!homogeneous_field<int>);
  static_assert(homogeneous_field<double>);

  static_assert(!field<Length>);
  static_assert(!field<int>);
  static_assert(field<double>);
}

TEST(Traits, LinearAlgebra) {
  static_assert(!homogeneous_vector_space<FixedMatrix<double, 3, 3>,
                                          FixedMatrix<double, 3, 3>>);
  static_assert(homogeneous_vector_space<Length, Length>);
  static_assert(!vector_space<Length, Length>);
  static_assert(vector_space<double, double>);

  static_assert(homogeneous_vector_space<FixedMatrix<Length, 3, 3>, Time>);
  static_assert(real_vector_space<FixedMatrix<Length, 3, 3>>);
  static_assert(!vector_space<FixedMatrix<Length, 3, 3>, Time>);
  static_assert(!real_vector_space<int>);
}

}  // namespace quantities
}  // namespace principia
