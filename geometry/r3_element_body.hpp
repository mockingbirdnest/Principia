
#pragma once

#include "geometry/r3_element.hpp"

#include <string>

#include "base/macros.hpp"
#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/serialization.hpp"

namespace principia {
namespace geometry {
namespace internal_r3_element {

using quantities::ArcSin;
using quantities::ArcTan;
using quantities::Cos;
using quantities::DebugString;
using quantities::DoubleOrQuantitySerializer;
using quantities::Quantity;
using quantities::Sin;
using quantities::SIUnit;

// We want zero initialization here, so the default constructor won't do.
template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ>::R3Element() : x(), y(), z() {}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ>::R3Element(ScalarX const& x,
                                                ScalarY const& y,
                                                ScalarZ const& z)
    : x(x), y(y), z(z) {}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
typename enable_if_normed<R3Element<ScalarX, ScalarY, ScalarZ>>::Scalar&
R3Element<ScalarX, ScalarY, ScalarZ>::operator[](int const index) {
  switch (index) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      DLOG(FATAL) << FUNCTION_SIGNATURE << ": index = " << index;
      base::noreturn();
  }
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
typename enable_if_normed<R3Element<ScalarX, ScalarY, ScalarZ>>::Scalar const&
R3Element<ScalarX, ScalarY, ScalarZ>::operator[](int const index) const {
  switch (index) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      DLOG(FATAL) << FUNCTION_SIGNATURE << ": index = " << index;
      base::noreturn();
  }
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ>& R3Element<ScalarX, ScalarY, ScalarZ>::
operator+=(R3Element<ScalarX, ScalarY, ScalarZ> const& right) {
  x += right.x;
  y += right.y;
  z += right.z;
  return *this;
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ>& R3Element<ScalarX, ScalarY, ScalarZ>::
operator-=(R3Element<ScalarX, ScalarY, ScalarZ> const& right) {
  x -= right.x;
  y -= right.y;
  z -= right.z;
  return *this;
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ>& R3Element<ScalarX, ScalarY, ScalarZ>::
operator*=(double const right) {
  x *= right;
  y *= right;
  z *= right;
  return *this;
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ>& R3Element<ScalarX, ScalarY, ScalarZ>::
operator/=(double const right) {
  x /= right;
  y /= right;
  z /= right;
  return *this;
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
typename enable_if_normed<R3Element<ScalarX, ScalarY, ScalarZ>>::Scalar
R3Element<ScalarX, ScalarY, ScalarZ>::Norm() const {
  return quantities::Sqrt(Dot(*this, *this));
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
template<typename S>
R3Element<
    typename enable_if_normed<R3Element<ScalarX, ScalarY, ScalarZ>>::Scalar>
R3Element<ScalarX, ScalarY, ScalarZ>::OrthogonalizationAgainst(
    R3Element<S> const& r3_element) const {
  R3Element<double> const r3_element_normalized = Normalize(r3_element);
  return *this - Dot(*this, r3_element_normalized) * r3_element_normalized;
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
SphericalCoordinates<
    typename enable_if_normed<R3Element<ScalarX, ScalarY, ScalarZ>>::Scalar>
R3Element<ScalarX, ScalarY, ScalarZ>::ToSpherical() const {
  SphericalCoordinates<
      typename enable_if_normed<R3Element<ScalarX, ScalarY, ScalarZ>>::Scalar>
      result;
  result.radius = Norm();
  result.latitude = ArcSin(z / result.radius);
  result.longitude = ArcTan(y, x);
  return result;
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
void R3Element<ScalarX, ScalarY, ScalarZ>::WriteToMessage(
    not_null<serialization::R3Element*> const message) const {
  using SerializerX =
      DoubleOrQuantitySerializer<ScalarX, serialization::R3Element::Coordinate>;
  using SerializerY =
      DoubleOrQuantitySerializer<ScalarY, serialization::R3Element::Coordinate>;
  using SerializerZ =
      DoubleOrQuantitySerializer<ScalarZ, serialization::R3Element::Coordinate>;
  SerializerX::WriteToMessage(x, message->mutable_x());
  SerializerY::WriteToMessage(y, message->mutable_y());
  SerializerZ::WriteToMessage(z, message->mutable_z());
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ>::ReadFromMessage(
    serialization::R3Element const& message) {
  using SerializerX =
      DoubleOrQuantitySerializer<ScalarX, serialization::R3Element::Coordinate>;
  using SerializerY =
      DoubleOrQuantitySerializer<ScalarY, serialization::R3Element::Coordinate>;
  using SerializerZ =
      DoubleOrQuantitySerializer<ScalarZ, serialization::R3Element::Coordinate>;
  return {SerializerX::ReadFromMessage(message.x()),
          SerializerY::ReadFromMessage(message.y()),
          SerializerZ::ReadFromMessage(message.z())};
}

template<typename Scalar>
SphericalCoordinates<Scalar>::SphericalCoordinates() {}

template<typename Scalar>
R3Element<Scalar> SphericalCoordinates<Scalar>::ToCartesian() {
  double const cos_latitude = Cos(latitude);
  return {radius * Cos(longitude) * cos_latitude,
          radius * Sin(longitude) * cos_latitude,
          radius * Sin(latitude)};
}

template<typename Scalar>
SphericalCoordinates<Scalar> RadiusLatitudeLongitude(Scalar const& radius,
                                                     Angle const& latitude,
                                                     Angle const& longitude) {
  SphericalCoordinates<Scalar> result;
  result.latitude = latitude;
  result.longitude = longitude;
  result.radius = radius;
  return result;
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator+(
    R3Element<ScalarX, ScalarY, ScalarZ> const& right) {
  return R3Element<ScalarX, ScalarY, ScalarZ>(+right.x, +right.y, +right.z);
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator-(
    R3Element<ScalarX, ScalarY, ScalarZ> const& right) {
  return R3Element<ScalarX, ScalarY, ScalarZ>(-right.x, -right.y, -right.z);
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator+(
    R3Element<ScalarX, ScalarY, ScalarZ> const& left,
    R3Element<ScalarX, ScalarY, ScalarZ> const& right) {
  return R3Element<ScalarX, ScalarY, ScalarZ>(left.x + right.x,
                                              left.y + right.y,
                                              left.z + right.z);
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator-(
    R3Element<ScalarX, ScalarY, ScalarZ> const& left,
    R3Element<ScalarX, ScalarY, ScalarZ> const& right) {
  return R3Element<ScalarX, ScalarY, ScalarZ>(left.x - right.x,
                                              left.y - right.y,
                                              left.z - right.z);
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator*(
    double const left,
    R3Element<ScalarX, ScalarY, ScalarZ> const& right) {
  return R3Element<ScalarX, ScalarY, ScalarZ>(left * right.x,
                                              left * right.y,
                                              left * right.z);
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator*(
    R3Element<ScalarX, ScalarY, ScalarZ> const& left,
    double const right) {
  return R3Element<ScalarX, ScalarY, ScalarZ>(left.x * right,
                                              left.y * right,
                                              left.z * right);
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
R3Element<ScalarX, ScalarY, ScalarZ> operator/(
    R3Element<ScalarX, ScalarY, ScalarZ> const& left,
    double const right) {
  return R3Element<ScalarX, ScalarY, ScalarZ>(left.x / right,
                                              left.y / right,
                                              left.z / right);
}

template<typename LDimension,
         typename RScalarX, typename RScalarY, typename RScalarZ>
R3Element<Product<Quantity<LDimension>, RScalarX>,
          Product<Quantity<LDimension>, RScalarY>,
          Product<Quantity<LDimension>, RScalarZ>>
operator*(Quantity<LDimension> const& left,
          R3Element<RScalarX, RScalarY, RScalarZ> const& right) {
  return R3Element<Product<Quantity<LDimension>, RScalarX>,
                   Product<Quantity<LDimension>, RScalarY>,
                   Product<Quantity<LDimension>, RScalarZ>>(
      left * right.x,
      left * right.y,
      left * right.z);
}

template<typename LScalarX, typename LScalarY, typename LScalarZ,
         typename RDimension>
R3Element<Product<LScalarX, Quantity<RDimension>>,
          Product<LScalarY, Quantity<RDimension>>,
          Product<LScalarZ, Quantity<RDimension>>>
operator*(R3Element<LScalarX, LScalarY, LScalarZ> const& left,
          Quantity<RDimension> const& right) {
  return R3Element<Product<LScalarX, Quantity<RDimension>>,
                   Product<LScalarY, Quantity<RDimension>>,
                   Product<LScalarZ, Quantity<RDimension>>>(
      left.x * right,
      left.y * right,
      left.z * right);
}

template<typename LScalarX, typename LScalarY, typename LScalarZ,
         typename RDimension>
R3Element<Quotient<LScalarX, Quantity<RDimension>>,
          Quotient<LScalarY, Quantity<RDimension>>,
          Quotient<LScalarZ, Quantity<RDimension>>>
operator/(R3Element<LScalarX, LScalarY, LScalarZ> const& left,
          Quantity<RDimension> const& right) {
  return R3Element<Quotient<LScalarX, Quantity<RDimension>>,
                   Quotient<LScalarY, Quantity<RDimension>>,
                   Quotient<LScalarZ, Quantity<RDimension>>>(
      left.x / right,
      left.y / right,
      left.z / right);
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
bool operator==(R3Element<ScalarX, ScalarY, ScalarZ> const& left,
                R3Element<ScalarX, ScalarY, ScalarZ> const& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
bool operator!=(R3Element<ScalarX, ScalarY, ScalarZ> const& left,
                R3Element<ScalarX, ScalarY, ScalarZ> const& right) {
  return left.x != right.x || left.y != right.y || left.z != right.z;
}

template<typename Scalar>
R3Element<double> Normalize(R3Element<Scalar> const& r3_element) {
  Scalar const norm = r3_element.Norm();
#ifdef _DEBUG
  CHECK_NE(Scalar(), norm);
#endif
  return r3_element / norm;
}

template<typename Scalar>
R3Element<double> NormalizeOrZero(R3Element<Scalar> const& r3_element) {
  Scalar const norm = r3_element.Norm();
  if (norm == Scalar()) {
    static R3Element<double> const zeroes = {0, 0, 0};
    return zeroes;
  } else {
    return r3_element / norm;
  }
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
std::string DebugString(
    R3Element<ScalarX, ScalarY, ScalarZ> const& r3_element) {
  std::string result = "{";
  result += DebugString(r3_element.x);
  result += ", ";
  result += DebugString(r3_element.y);
  result += ", ";
  result += DebugString(r3_element.z);
  result +="}";
  return result;
}

template<typename ScalarX, typename ScalarY, typename ScalarZ>
std::ostream& operator<<(
    std::ostream& out,
    R3Element<ScalarX, ScalarY, ScalarZ> const& r3_element) {
  out << DebugString(r3_element);
  return out;
}

template<typename LScalarX,
         typename RScalarX,
         typename LScalarY,
         typename RScalarY,
         typename LScalarZ,
         typename RScalarZ>
R3Element<Product<LScalarX, RScalarX>,
          Product<LScalarY, RScalarY>,
          Product<LScalarZ, RScalarZ>>
Cross(R3Element<LScalarX, LScalarY, LScalarZ> const& left,
      R3Element<RScalarX, RScalarY, RScalarZ> const& right) {
  return R3Element<Product<LScalarX, RScalarX>,
                   Product<LScalarY, RScalarY>,
                   Product<LScalarZ, RScalarZ>>(
      left.y * right.z - left.z * right.y,
      left.z * right.x - left.x * right.z,
      left.x * right.y - left.y * right.x);
}

template<typename LScalar, typename RScalar>
Product<LScalar, RScalar> Dot(R3Element<LScalar> const& left,
                              R3Element<RScalar> const& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

inline R3Element<double> BasisVector(int const i) {
  DCHECK_GE(i, 0) << i;
  DCHECK_LT(i, 3) << i;
  return {static_cast<double>(i == 0),
          static_cast<double>(i == 1),
          static_cast<double>(i == 2)};
}

}  // namespace internal_r3_element
}  // namespace geometry
}  // namespace principia
