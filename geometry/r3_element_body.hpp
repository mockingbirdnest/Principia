
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
template<typename Scalar>
R3Element<Scalar>::R3Element() : x(), y(), z() {}

template<typename Scalar>
R3Element<Scalar>::R3Element(Scalar const& x,
                             Scalar const& y,
                             Scalar const& z) : x(x), y(y), z(z) {}

template<typename Scalar>
Scalar& R3Element<Scalar>::operator[](int const index) {
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

template<typename Scalar>
Scalar const& R3Element<Scalar>::operator[](
    int const index) const {
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

template<typename Scalar>
R3Element<Scalar>& R3Element<Scalar>::operator+=(
    R3Element<Scalar> const& right) {
  x += right.x;
  y += right.y;
  z += right.z;
  return *this;
}

template<typename Scalar>
R3Element<Scalar>& R3Element<Scalar>::operator-=(
    R3Element<Scalar> const& right) {
  x -= right.x;
  y -= right.y;
  z -= right.z;
  return *this;
}

template<typename Scalar>
R3Element<Scalar>& R3Element<Scalar>::operator*=(double const right) {
  x *= right;
  y *= right;
  z *= right;
  return *this;
}

template<typename Scalar>
R3Element<Scalar>& R3Element<Scalar>::operator/=(double const right) {
  x /= right;
  y /= right;
  z /= right;
  return *this;
}

template<typename Scalar>
Scalar R3Element<Scalar>::Norm() const {
  return quantities::Sqrt(Dot(*this, *this));
}

template<typename Scalar>
SphericalCoordinates<Scalar> R3Element<Scalar>::ToSpherical() const {
  SphericalCoordinates<Scalar> result;
  result.radius = Norm();
  result.latitude = ArcSin(z / result.radius);
  result.longitude = ArcTan(y, x);
  return result;
}

template<typename Scalar>
template<typename S>
R3Element<Scalar> R3Element<Scalar>::OrthogonalizationAgainst(
    R3Element<S> const& r3_element) const {
  R3Element<double> const r3_element_normalized = Normalize(r3_element);
  return *this - Dot(*this, r3_element_normalized) * r3_element_normalized;
}

template<typename Scalar>
void R3Element<Scalar>::WriteToMessage(
    not_null<serialization::R3Element*> const message) const {
  using Serializer =
      DoubleOrQuantitySerializer<Scalar, serialization::R3Element::Coordinate>;
  Serializer::WriteToMessage(x, message->mutable_x());
  Serializer::WriteToMessage(y, message->mutable_y());
  Serializer::WriteToMessage(z, message->mutable_z());
}

template<typename Scalar>
R3Element<Scalar> R3Element<Scalar>::ReadFromMessage(
    serialization::R3Element const& message) {
  using Serializer =
      DoubleOrQuantitySerializer<Scalar, serialization::R3Element::Coordinate>;
  return {Serializer::ReadFromMessage(message.x()),
          Serializer::ReadFromMessage(message.y()),
          Serializer::ReadFromMessage(message.z())};
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

template<typename Scalar>
R3Element<Scalar> operator+(R3Element<Scalar> const& right) {
  return R3Element<Scalar>(+right.x, +right.y, +right.z);
}

template<typename Scalar>
R3Element<Scalar> operator-(R3Element<Scalar> const& right) {
  return R3Element<Scalar>(-right.x, -right.y, -right.z);
}

template<typename Scalar>
R3Element<Scalar> operator+(R3Element<Scalar> const& left,
                            R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left.x + right.x,
                           left.y + right.y,
                           left.z + right.z);
}

template<typename Scalar>
R3Element<Scalar> operator-(R3Element<Scalar> const& left,
                            R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left.x - right.x,
                           left.y - right.y,
                           left.z - right.z);
}

template<typename Scalar>
R3Element<Scalar> operator*(double const left,
                            R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left * right.x,
                           left * right.y,
                           left * right.z);
}

template<typename Scalar>
R3Element<Scalar> operator*(R3Element<Scalar> const& left,
                            double const right) {
  return R3Element<Scalar>(left.x * right,
                           left.y * right,
                           left.z * right);
}

template<typename Scalar>
R3Element<Scalar> operator/(R3Element<Scalar> const& left,
                            double const right) {
  return R3Element<Scalar>(left.x / right,
                           left.y / right,
                           left.z / right);
}

template<typename LDimension, typename RScalar>
R3Element<Product<Quantity<LDimension>, RScalar>>
operator*(Quantity<LDimension> const& left, R3Element<RScalar> const& right) {
  return R3Element<Product<Quantity<LDimension>, RScalar>>(
      left * right.x,
      left * right.y,
      left * right.z);
}

template<typename LScalar, typename RDimension>
R3Element<Product<LScalar, Quantity<RDimension>>>
operator*(R3Element<LScalar> const& left, Quantity<RDimension> const& right) {
  return R3Element<Product<LScalar, Quantity<RDimension>>>(
      left.x * right,
      left.y * right,
      left.z * right);
}

template<typename LScalar, typename RDimension>
R3Element<Quotient<LScalar, Quantity<RDimension>>>
operator/(R3Element<LScalar> const& left,
          Quantity<RDimension> const& right) {
  return R3Element<Quotient<LScalar, Quantity<RDimension>>>(
      left.x / right,
      left.y / right,
      left.z / right);
}

template<typename Scalar>
bool operator==(R3Element<Scalar> const& left,
                R3Element<Scalar> const& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

template<typename Scalar>
bool operator!=(R3Element<Scalar> const& left,
                R3Element<Scalar> const& right) {
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

template<typename Scalar>
std::string DebugString(R3Element<Scalar> const& r3_element) {
  std::string result = "{";
  result += DebugString(r3_element.x);
  result += ", ";
  result += DebugString(r3_element.y);
  result += ", ";
  result += DebugString(r3_element.z);
  result +="}";
  return result;
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         R3Element<Scalar> const& r3_element) {
  out << DebugString(r3_element);
  return out;
}

template<typename LScalar, typename RScalar>
R3Element<Product<LScalar, RScalar>> Cross(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right) {
  return R3Element<Product<LScalar, RScalar>>(
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
