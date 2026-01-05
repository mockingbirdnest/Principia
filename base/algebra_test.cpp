#include "base/algebra.hpp"

#include <string>

#include "base/multiprecision.hpp"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace base {

using namespace principia::base::_algebra;
using namespace principia::base::_multiprecision;
using namespace principia::quantities::_quantities;

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
  static_assert(ring<cpp_int>);

  static_assert(homogeneous_field<Length>);
  static_assert(!homogeneous_field<int>);
  static_assert(homogeneous_field<double>);

  static_assert(!field<Length>);
  static_assert(!field<int>);
  static_assert(field<double>);
  static_assert(!field<cpp_int>);
  static_assert(field<cpp_rational>);
  static_assert(field<cpp_bin_float_50>);
}

TEST(Concepts, LinearAlgebra) {
  static_assert(homogeneous_vector_space<Length, Length>);
  static_assert(!vector_space<Length, Length>);
  static_assert(vector_space<double, double>);

  static_assert(!real_vector_space<int>);
}

}  // namespace base
}  // namespace principia
