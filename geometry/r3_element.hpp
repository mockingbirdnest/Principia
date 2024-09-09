#pragma once

#include <pmmintrin.h>

#include <iostream>
#include <string>

#include "base/not_null.hpp"
#include "base/tags.hpp"
#include "quantities/concepts.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace _r3_element {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::base::_tags;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

template<typename Scalar>
struct SphericalCoordinates;

// An `R3Element<Scalar>` is an element of Scalar³. `Scalar` should be a vector
// space over ℝ, represented by `double`. `R3Element` is the underlying data
// type for more advanced strongly typed structures suchas `Multivector`.
template<typename Scalar>
struct alignas(16) R3Element final {
 public:
  constexpr R3Element();
  constexpr explicit R3Element(uninitialized_t);
  R3Element(Scalar const& x, Scalar const& y, Scalar const& z);
  R3Element(__m128d xy, __m128d zt);

  Scalar& operator[](int index);
  Scalar const& operator[](int index) const;

  R3Element& operator+=(R3Element const& right);
  R3Element& operator-=(R3Element const& right);
  R3Element& operator*=(double right);
  R3Element& operator/=(double right);

  Scalar Norm() const;
  Square<Scalar> Norm²() const;

  // Uses the x-y plane as the equator, the x-axis as the reference direction on
  // the equator, and the z-axis as the north pole.
  SphericalCoordinates<Scalar> ToSpherical() const;

  // Returns a vector coplanar to `*this` and `r3_element`, orthogonal to
  // `r3_element`, and on the same side of `r3_element` as `*this`.
  // Uses the modified Gram-Schmidt algorithm.  Fails if `r3_element` is zero.
  template<typename S>
  R3Element OrthogonalizationAgainst(R3Element<S> const& r3_element) const;

  void WriteToMessage(not_null<serialization::R3Element*> message) const;
  static R3Element ReadFromMessage(serialization::R3Element const& message);

  union {
    struct {
      Scalar x;
      Scalar y;
      Scalar z;
    };
    struct {
      __m128d xy;
      __m128d zt;
    };
  };
};

template<typename Scalar>
struct SphericalCoordinates final {
  // Default, but prevents aggregate initialization of `SphericalCoordinates` to
  // obviate confusion over the order of `latitude` and `longitude`.
  SphericalCoordinates();

  // Uses the x-y plane as the equator, the x-axis as the reference direction on
  // the equator, and the z-axis as the north pole.
  R3Element<Scalar> ToCartesian();

  Scalar radius;
  Angle latitude;
  Angle longitude;
};

template<typename Scalar>
SphericalCoordinates<Scalar> RadiusLatitudeLongitude(Scalar const& radius,
                                                     Angle const& latitude,
                                                     Angle const& longitude);

template<typename Scalar>
std::ostream& operator<<(
    std::ostream& out,
    SphericalCoordinates<Scalar> const& spherical_coordinates);

template<typename Scalar>
R3Element<Scalar> operator+(R3Element<Scalar> const& right);
template<typename Scalar>
R3Element<Scalar> operator-(R3Element<Scalar> const& right);

template<typename Scalar>
R3Element<Scalar> operator+(R3Element<Scalar> const& left,
                            R3Element<Scalar> const& right);
template<typename Scalar>
R3Element<Scalar> operator-(R3Element<Scalar> const& left,
                            R3Element<Scalar> const& right);

// Dimensionful multiplication `LScalar * R3Element<RScalar>` is the tensor
// product LScalar ⊗ Scalar³. Since LScalar ⊗ Scalar³ ≅ (LScalar ⊗ Scalar)³,
// the result is an R3Element<Product<LScalar, RScalar>>.
template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<LScalar>
R3Element<Product<LScalar, RScalar>> operator*(LScalar const& left,
                                               R3Element<RScalar> const& right);

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3Element<Product<LScalar, RScalar>> operator*(R3Element<LScalar> const& left,
                                               RScalar const& right);

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3Element<Quotient<LScalar, RScalar>> operator/(R3Element<LScalar> const& left,
                                                RScalar const& right);

// FMA for ±vector * scalar ± vector.
template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3Element<Product<LScalar, RScalar>> FusedMultiplyAdd(
    R3Element<LScalar> const& a,
    RScalar const& b,
    R3Element<Product<LScalar, RScalar>> const& c);
template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3Element<Product<LScalar, RScalar>> FusedMultiplySubtract(
    R3Element<LScalar> const& a,
    RScalar const& b,
    R3Element<Product<LScalar, RScalar>> const& c);
template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3Element<Product<LScalar, RScalar>> FusedNegatedMultiplyAdd(
    R3Element<LScalar> const& a,
    RScalar const& b,
    R3Element<Product<LScalar, RScalar>> const& c);
template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3Element<Product<LScalar, RScalar>> FusedNegatedMultiplySubtract(
    R3Element<LScalar> const& a,
    RScalar const& b,
    R3Element<Product<LScalar, RScalar>> const& c);

// FMA for ±scalar * vector ± vector.
template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<LScalar>
R3Element<Product<LScalar, RScalar>> FusedMultiplyAdd(
    LScalar const& a,
    R3Element<RScalar> const& b,
    R3Element<Product<LScalar, RScalar>> const& c);
template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<LScalar>
R3Element<Product<LScalar, RScalar>> FusedMultiplySubtract(
    LScalar const& a,
    R3Element<RScalar> const& b,
    R3Element<Product<LScalar, RScalar>> const& c);
template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<LScalar>
R3Element<Product<LScalar, RScalar>> FusedNegatedMultiplyAdd(
    LScalar const& a,
    R3Element<RScalar> const& b,
    R3Element<Product<LScalar, RScalar>> const& c);
template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<LScalar>
R3Element<Product<LScalar, RScalar>> FusedNegatedMultiplySubtract(
    LScalar const& a,
    R3Element<RScalar> const& b,
    R3Element<Product<LScalar, RScalar>> const& c);

template<typename Scalar>
constexpr bool operator==(R3Element<Scalar> const& left,
                          R3Element<Scalar> const& right);
template<typename Scalar>
constexpr bool operator!=(R3Element<Scalar> const& left,
                          R3Element<Scalar> const& right);

template<typename Scalar>
R3Element<double> Normalize(R3Element<Scalar> const& r3_element);

template<typename Scalar>
R3Element<double> NormalizeOrZero(R3Element<Scalar> const& r3_element);

template<typename Scalar>
std::string DebugString(R3Element<Scalar> const& r3_element);

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         R3Element<Scalar> const& r3_element);

template<typename LScalar, typename RScalar>
R3Element<Product<LScalar, RScalar>> Cross(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right);

template<typename LScalar, typename RScalar>
Product<LScalar, RScalar> Dot(R3Element<LScalar> const& left,
                              R3Element<RScalar> const& right);

// Returns the `i`th basis vector, whose `i`th coordinate is 1, and whose
// other coordinates are 0.  `i` must be in [0, 2].
R3Element<double> BasisVector(int i);

}  // namespace internal

using internal::BasisVector;
using internal::Cross;
using internal::Dot;
using internal::FusedNegatedMultiplyAdd;
using internal::FusedNegatedMultiplySubtract;
using internal::FusedMultiplyAdd;
using internal::FusedMultiplySubtract;
using internal::Normalize;
using internal::NormalizeOrZero;
using internal::R3Element;
using internal::RadiusLatitudeLongitude;
using internal::SphericalCoordinates;

}  // namespace _r3_element
}  // namespace geometry
}  // namespace principia

#include "geometry/r3_element_body.hpp"
