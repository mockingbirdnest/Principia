
#pragma once

#include <pmmintrin.h>

#include <iostream>
#include <string>

#include "base/not_null.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace internal_r3_element {

using base::not_null;
using quantities::Angle;
using quantities::is_quantity;
using quantities::Product;
using quantities::Quantity;
using quantities::Quotient;
using quantities::Square;

template<typename Scalar>
struct SphericalCoordinates;

// An |R3Element<Scalar>| is an element of Scalar³. |Scalar| should be a vector
// space over ℝ, represented by |double|. |R3Element| is the underlying data
// type for more advanced strongly typed structures suchas |Multivector|.
template<typename Scalar>
struct alignas(16) R3Element final {
 public:
  constexpr R3Element();
  R3Element(Scalar const& x, Scalar const& y, Scalar const& z);
  R3Element(__m128d xy, __m128d zt);

  Scalar&       operator[](int index);
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

  // Returns a vector coplanar to |*this| and |r3_element|, orthogonal to
  // |r3_element|, and on the same side of |r3_element| as |*this|.
  // Uses the modified Gram-Schmidt algorithm.  Fails if |r3_element| is zero.
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
  // Default, but prevents aggregate initialization of |SphericalCoordinates| to
  // obviate confusion over the order of |latitude| and |longitude|.
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
R3Element<Scalar> operator+(R3Element<Scalar> const& right);
template<typename Scalar>
R3Element<Scalar> operator-(R3Element<Scalar> const& right);

template<typename Scalar>
R3Element<Scalar> operator+(R3Element<Scalar> const& left,
                            R3Element<Scalar> const& right);
template<typename Scalar>
R3Element<Scalar> operator-(R3Element<Scalar> const& left,
                            R3Element<Scalar> const& right);

// Dimensionful multiplication |LScalar * R3Element<RScalar>| is the tensor
// product LScalar ⊗ Scalar³. Since LScalar ⊗ Scalar³ ≅ (LScalar ⊗ Scalar)³,
// the result is an R3Element<Product<LScalar, RScalar>>.
// The special case where one of the scalars is |double| is handled separately
// above in order to allow implicit conversions to |double|.
template<typename LScalar, typename RScalar,
         typename = std::enable_if_t<is_quantity<LScalar>::value>>
R3Element<Product<LScalar, RScalar>>
operator*(LScalar const& left, R3Element<RScalar> const& right);

template<typename LScalar, typename RScalar,
         typename = std::enable_if_t<is_quantity<RScalar>::value>>
R3Element<Product<LScalar, RScalar>>
operator*(R3Element<LScalar> const& left, RScalar const& right);

template<typename LScalar, typename RScalar,
         typename = std::enable_if_t<is_quantity<RScalar>::value>>
R3Element<Quotient<LScalar, RScalar>>
operator/(R3Element<LScalar> const& left, RScalar const& right);

template<typename Scalar>
bool operator==(R3Element<Scalar> const& left,
                R3Element<Scalar> const& right);
template<typename Scalar>
bool operator!=(R3Element<Scalar> const& left,
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

// Returns the |i|th basis vector, whose |i|th coordinate is 1, and whose
// other coordinates are 0.  |i| must be in [0, 2].
R3Element<double> BasisVector(int i);

}  // namespace internal_r3_element

using internal_r3_element::BasisVector;
using internal_r3_element::Cross;
using internal_r3_element::Dot;
using internal_r3_element::Normalize;
using internal_r3_element::NormalizeOrZero;
using internal_r3_element::R3Element;
using internal_r3_element::RadiusLatitudeLongitude;
using internal_r3_element::SphericalCoordinates;

}  // namespace geometry
}  // namespace principia

#include "geometry/r3_element_body.hpp"
