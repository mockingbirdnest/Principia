#pragma once

#include <string>

#include "base/macros.hpp"
#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/serialization.hpp"
#include "r3_element.hpp"

namespace principia {

using quantities::Cos;
using quantities::Sin;
using quantities::DoubleOrQuantitySerializer;
using quantities::Quantity;
using quantities::SIUnit;

namespace geometry {

// We want zero initialization here, so the default constructor won't do.
template<typename Scalar>
inline R3Element<Scalar>::R3Element() : x(), y(), z() {}

template<typename Scalar>
inline R3Element<Scalar>::R3Element(Scalar const& x,
                                    Scalar const& y,
                                    Scalar const& z) : x(x), y(y), z(z) {}

template<typename Scalar>
inline Scalar& R3Element<Scalar>::operator[](int const index) {
  switch (index) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      LOG(FATAL) << FUNCTION_SIGNATURE << ": index = " << index;
      base::noreturn();
  }
}

template<typename Scalar>
inline Scalar const& R3Element<Scalar>::operator[](int const index) const {
  switch (index) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      LOG(FATAL) << FUNCTION_SIGNATURE << ": index = " << index;
      base::noreturn();
  }
}

template<typename Scalar>
inline R3Element<Scalar>& R3Element<Scalar>::operator+=(
    R3Element<Scalar> const& right) {
  x += right.x;
  y += right.y;
  z += right.z;
  return *this;
}

template<typename Scalar>
inline R3Element<Scalar>& R3Element<Scalar>::operator-=(
    R3Element<Scalar> const& right) {
  x -= right.x;
  y -= right.y;
  z -= right.z;
  return *this;
}

template<typename Scalar>
inline R3Element<Scalar>& R3Element<Scalar>::operator*=(double const right) {
  x *= right;
  y *= right;
  z *= right;
  return *this;
}

template<typename Scalar>
inline R3Element<Scalar>& R3Element<Scalar>::operator/=(double const right) {
  x /= right;
  y /= right;
  z /= right;
  return *this;
}

template<typename Scalar>
inline Scalar R3Element<Scalar>::Norm() const {
  return quantities::Sqrt(Dot(*this, *this));
}

template<typename Scalar>
template<typename S>
void R3Element<Scalar>::Orthogonalize(
    not_null<R3Element<S>*> const r3_element) const {
  Scalar const this_norm = this->Norm();
  CHECK_NE(0 * SIUnit<Scalar>(), this_norm);
  R3Element<double> const this_normalized = *this / this_norm;
  *r3_element -= Dot(*r3_element, this_normalized) * this_normalized;
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
SphericalCoordinates<Scalar>::SphericalCoordinates() {};

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
R3Element<Scalar> ToCartesian(SphericalCoordinates<Scalar> const& coordinates) {
  double cos_? = Cos(coordinates.latitude);
  double sin_? = Sin(coordinates.latitude);
  return {coordinates.radius * Cos(coordinates.longitude) * cos_?,
          coordinates.radius * Sin(coordinates.longitude) * cos_?,
          coordinates.radius * sin_?};
}

template<typename Scalar>
inline R3Element<Scalar> operator+(R3Element<Scalar> const& right) {
  return R3Element<Scalar>(+right.x, +right.y, +right.z);
}

template<typename Scalar>
inline R3Element<Scalar> operator-(R3Element<Scalar> const& right) {
  return R3Element<Scalar>(-right.x, -right.y, -right.z);
}

template<typename Scalar>
inline R3Element<Scalar> operator+(
    R3Element<Scalar> const& left,
    R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left.x + right.x,
                           left.y + right.y,
                           left.z + right.z);
}

template<typename Scalar>
inline R3Element<Scalar> operator-(
    R3Element<Scalar> const& left,
    R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left.x - right.x,
                           left.y - right.y,
                           left.z - right.z);
}

template<typename Scalar>
inline R3Element<Scalar> operator*(double const left,
                                   R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left * right.x,
                           left * right.y,
                           left * right.z);
}

template<typename Scalar>
inline R3Element<Scalar> operator*(R3Element<Scalar> const& left,
                                   double const right) {
  return R3Element<Scalar>(left.x * right,
                           left.y * right,
                           left.z * right);
}

template<typename Scalar>
inline R3Element<Scalar> operator/(R3Element<Scalar> const& left,
                                   double const right) {
  return R3Element<Scalar>(left.x / right,
                           left.y / right,
                           left.z / right);
}

template<typename LDimension, typename RScalar>
inline R3Element<quantities::Product<quantities::Quantity<LDimension>, RScalar>>
operator*(quantities::Quantity<LDimension> const& left,
          R3Element<RScalar> const& right) {
  return R3Element<quantities::Product<quantities::Quantity<LDimension>,
                                       RScalar>>(
      left * right.x,
      left * right.y,
      left * right.z);
}

template<typename LScalar, typename RDimension>
inline R3Element<quantities::Product<LScalar, quantities::Quantity<RDimension>>>
operator*(R3Element<LScalar> const& left,
          quantities::Quantity<RDimension> const& right) {
  return R3Element<quantities::Product<LScalar,
                                       quantities::Quantity<RDimension>>>(
      left.x * right,
      left.y * right,
      left.z * right);
}

template<typename LScalar, typename RDimension>
inline R3Element<quantities::Quotient<LScalar,
                                      quantities::Quantity<RDimension>>>
operator/(R3Element<LScalar> const& left,
          quantities::Quantity<RDimension> const& right) {
  return R3Element<quantities::Quotient<LScalar,
                                        quantities::Quantity<RDimension>>>(
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
  return r3_element / r3_element.Norm();
}

template<typename Scalar>
std::string DebugString(R3Element<Scalar> const& r3_element) {
  std::string result = "{";
  result += quantities::DebugString(r3_element.x);
  result += ", ";
  result += quantities::DebugString(r3_element.y);
  result += ", ";
  result += quantities::DebugString(r3_element.z);
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
inline R3Element<quantities::Product<LScalar, RScalar>> Cross(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right) {
  return R3Element<quantities::Product<LScalar, RScalar>>(
      left.y * right.z - left.z * right.y,
      left.z * right.x - left.x * right.z,
      left.x * right.y - left.y * right.x);
}

template<typename LScalar, typename RScalar>
inline quantities::Product<LScalar, RScalar> Dot(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

}  // namespace geometry
}  // namespace principia
