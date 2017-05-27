
#pragma once

// We use ostream for logging purposes.
#include <iostream>  // NOLINT(readability/streams)
#include <string>

#include "base/not_constructible.hpp"
#include "base/not_null.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace internal_r3_element {

using base::not_constructible;
using base::not_null;
using quantities::Angle;
using quantities::Product;
using quantities::Quantity;
using quantities::Quotient;

template<typename Scalar>
struct SphericalCoordinates;

template<typename ScalarX,
         typename ScalarY = ScalarX,
         typename ScalarZ = ScalarX>
struct R3Element;

template<typename T>
struct enable_if_normed : not_constructible {};

template<typename T>
struct enable_if_normed<R3Element<T, T, T>> : not_constructible {
  using Scalar = T;
};

// An |R3Element| is an element of ℝ³. |ScalarX|, |ScalarY| and |ScalarZ| must
// be 1-dimensional vector spaces over ℝ, typically represented by |Quantity|
// or |double|.  If the three scalars are homogeneous (in the dimensional
// analysis sense, not in the projective sense), then the resulting vector space
// is also normed.
// |R3Element| is the underlying data type for more advanced strongly typed
// structures such as |Multivector|.
template<typename ScalarX, typename ScalarY, typename ScalarZ>
struct R3Element final {
 public:
  R3Element();
  R3Element(ScalarX const& x, ScalarY const& y, ScalarZ const& z);

  template<typename Scalar =
               enable_if_normed<R3Element<ScalarX, ScalarY, ScalarZ>>::Scalar>
  Scalar& operator[](int index);
  template<typename Scalar =
               enable_if_normed<R3Element<ScalarX, ScalarY, ScalarZ>>::Scalar>
  Scalar const& operator[](int index) const;

  R3Element& operator+=(R3Element const& right);
  R3Element& operator-=(R3Element const& right);
  R3Element& operator*=(double right);
  R3Element& operator/=(double right);

  template<typename Scalar =
               enable_if_normed<R3Element<ScalarX, ScalarY, ScalarZ>>::Scalar>
  Scalar Norm() const;

  // Returns a vector coplanar to |*this| and |r3_element|, orthogonal to
  // |r3_element|, and on the same side of |r3_element| as |*this|.
  // Uses the modified Gram-Schmidt algorithm.  Fails if |r3_element| is zero.
  template<typename S,
           typename Scalar =
               enable_if_normed<R3Element<ScalarX, ScalarY, ScalarZ>>::Scalar>
  R3Element<Scalar> OrthogonalizationAgainst(
      R3Element<S> const& r3_element) const;

  // Uses the x-y plane as the equator, the x-axis as the reference direction on
  // the equator, and the z-axis as the north pole.
  template<typename Scalar =
               enable_if_normed<R3Element<ScalarX, ScalarY, ScalarZ>>::Scalar>
  SphericalCoordinates<Scalar> ToSpherical() const;

  void WriteToMessage(not_null<serialization::R3Element*> message) const;
  static R3Element ReadFromMessage(serialization::R3Element const& message);

  ScalarX x;
  ScalarY y;
  ScalarZ z;
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

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator+(
    R3Element<ScalarX, ScalarY, ScalarZ> const& right);
template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator-(
    R3Element<ScalarX, ScalarY, ScalarZ> const& right);

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator+(
    R3Element<ScalarX, ScalarY, ScalarZ> const& left,
    R3Element<ScalarX, ScalarY, ScalarZ> const& right);
template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator-(
    R3Element<ScalarX, ScalarY, ScalarZ> const& left,
    R3Element<ScalarX, ScalarY, ScalarZ> const& right);

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator*(
    double left,
    R3Element<ScalarX, ScalarY, ScalarZ> const& right);
template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator*(
    R3Element<ScalarX, ScalarY, ScalarZ> const& left,
    double right);
template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator/(
    R3Element<ScalarX, ScalarY, ScalarZ> const& left,
    double right);

// Dimensionful multiplication |LScalar * R3Element<RScalar>| is the tensor
// product LScalar ⊗ Scalar³. Since LScalar ⊗ Scalar³ ≅ (LScalar ⊗ Scalar)³,
// the result is an R3Element<Product<LScalar, RScalar>>.
// The special case where one of the scalars is |double| is handled separately
// above in order to allow implicit conversions to |double|.
template<typename LDimension,
         typename RScalarX, typename RScalarY, typename RScalarZ>
R3Element<Product<Quantity<LDimension>, RScalarX>,
          Product<Quantity<LDimension>, RScalarY>,
          Product<Quantity<LDimension>, RScalarZ>>
operator*(Quantity<LDimension> const& left,
          R3Element<RScalarX, RScalarY, RScalarZ> const& right);
template<typename LScalarX, typename LScalarY, typename LScalarZ,
         typename RDimension>
R3Element<Product<LScalarX, Quantity<RDimension>>,
          Product<LScalarY, Quantity<RDimension>>,
          Product<LScalarZ, Quantity<RDimension>>>
operator*(R3Element<LScalarX, LScalarY, LScalarZ> const& left,
          Quantity<RDimension> const& right);
template<typename LScalarX, typename LScalarY, typename LScalarZ,
         typename RDimension>
R3Element<Quotient<LScalarX, Quantity<RDimension>>,
          Quotient<LScalarY, Quantity<RDimension>>,
          Quotient<LScalarZ, Quantity<RDimension>>>
operator/(R3Element<LScalarX, LScalarY, LScalarZ> const& left,
          Quantity<RDimension> const& right);

template<typename ScalarX, typename ScalarY, typename ScalarZ>
bool operator==(R3Element<ScalarX, ScalarY, ScalarZ> const& left,
                R3Element<ScalarX, ScalarY, ScalarZ> const& right);
template<typename ScalarX, typename ScalarY, typename ScalarZ>
bool operator!=(R3Element<ScalarX, ScalarY, ScalarZ> const& left,
                R3Element<ScalarX, ScalarY, ScalarZ> const& right);

template<typename Scalar>
R3Element<double> Normalize(R3Element<Scalar> const& r3_element);

template<typename Scalar>
R3Element<double> NormalizeOrZero(R3Element<Scalar> const& r3_element);

template<typename ScalarX, typename ScalarY, typename ScalarZ>
std::string DebugString(R3Element<ScalarX, ScalarY, ScalarZ> const& r3_element);

template<typename ScalarX, typename ScalarY, typename ScalarZ>
std::ostream& operator<<(
    std::ostream& out,
    R3Element<ScalarX, ScalarY, ScalarZ> const& r3_element);

// Beware, the template arguments are in nonstandard order to make defaulting
// work.
template<typename LScalarX,
         typename RScalarX,
         typename LScalarY = LScalarX,
         typename RScalarY = RScalarX,
         typename LScalarZ = LScalarX,
         typename RScalarZ = RScalarX>
R3Element<Product<LScalarX, RScalarX>,
          Product<LScalarY, RScalarY>,
          Product<LScalarZ, RScalarZ>>
Cross(R3Element<LScalarX, LScalarY, LScalarZ> const& left,
      R3Element<RScalarX, RScalarY, RScalarZ> const& right);

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
