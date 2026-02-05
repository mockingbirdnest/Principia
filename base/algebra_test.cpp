#include "base/algebra.hpp"

#include <string>

#include "base/multiprecision.hpp"
#include "geometry/instant.hpp"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "testing_utilities/check_well_formedness.hpp"

namespace principia {
namespace base {

using namespace principia::base::_algebra;
using namespace principia::base::_multiprecision;
using namespace principia::geometry::_instant;
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

template<affine T>
constexpr std::string_view description() {
  return "affine";
}
template<additive_group T>
constexpr std::string_view description() {
  return "additive group";
}
template<homogeneous_ring T>
constexpr std::string_view description() {
  return "homogeneous ring";
}
template<ring T>
constexpr std::string_view description() {
  return "ring";
}
template<homogeneous_field T>
constexpr std::string_view description() {
  return "homogeneous field";
}
template<field T>
constexpr std::string_view description() {
  return "field";
}
template<affine T>
constexpr std::string_view vector_description() {
  return "affine";
}
template<module<int> T>
constexpr std::string_view vector_description() {
  return "ℤ-module";
}
template<real_affine_space T>
constexpr std::string_view vector_description() {
  return "real affine space";
}
template<real_vector_space T>
constexpr std::string_view vector_description() {
  return "real vector space";
}

template<additive_group T>
constexpr std::string_view AnAbelianGroupByAnyOtherName() {
  return "additive group";
}
template<module<int> T>
constexpr std::string_view AnAbelianGroupByAnyOtherName() {
  return "ℤ-module";
}

constexpr std::string_view meow = AnAbelianGroupByAnyOtherName<int>();

TEST(Concepts, Subsumption) {
  static_assert(description<Instant>() == "affine");
  static_assert(description<int>() == "ring");
  static_assert(description<Length>() == "homogeneous field");
  static_assert(description<double>() == "field");
  static_assert(vector_description<int>() == "ℤ-module");
  static_assert(vector_description<Instant>() == "real affine space");
  static_assert(vector_description<double>() == "real vector space");
}

}  // namespace base
}  // namespace principia
