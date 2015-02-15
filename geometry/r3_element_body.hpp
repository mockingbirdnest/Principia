#pragma once

#include <string>

#include "base/macros.hpp"
#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {

using quantities::Quantity;
using quantities::SIUnit;

namespace geometry {

template<typename T>
class R3ElementSerializer {};

template<typename Dimensions>
class R3ElementSerializer<Quantity<Dimensions>> {
 public:
  using T = Quantity<Dimensions>;
  static void WriteToMessage(
      T const& t,
      not_null<serialization::R3Element::Coordinate*> const message) {
    t.WriteToMessage(message->mutable_quantity());
  }

  static T ReadFromMessage(
      serialization::R3Element::Coordinate const& message) {
    CHECK(message.has_quantity());
    return T::ReadFromMessage(message.quantity());
  }
};

template<>
class R3ElementSerializer<double> {
 public:
  static void WriteToMessage(
      double const& d,
      not_null<serialization::R3Element::Coordinate*> const message) {
    message->set_double_(d);
  }

  static double ReadFromMessage(
      serialization::R3Element::Coordinate const& message) {
    CHECK(message.has_double_());
    return message.double_();
  }
};

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
      LOG(FATAL) << __FUNCSIG__ << ": index = " << index;
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
      LOG(FATAL) << __FUNCSIG__ << ": index = " << index;
      base::noreturn();
  }
}

template<typename Scalar>
inline R3Element<Scalar>& R3Element<Scalar>::operator+=(
    R3Element<Scalar> const& right) {
  return *this = *this + right;
}

template<typename Scalar>
inline R3Element<Scalar>& R3Element<Scalar>::operator-=(
    R3Element<Scalar> const& right) {
  return *this = *this - right;
}

template<typename Scalar>
inline R3Element<Scalar>& R3Element<Scalar>::operator*=(double const right) {
  return *this = *this * right;
}

template<typename Scalar>
inline R3Element<Scalar>& R3Element<Scalar>::operator/=(double const right) {
  return *this = *this / right;
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
  R3ElementSerializer<Scalar>::WriteToMessage(x, message->mutable_x());
  R3ElementSerializer<Scalar>::WriteToMessage(y, message->mutable_y());
  R3ElementSerializer<Scalar>::WriteToMessage(z, message->mutable_z());
}

template<typename Scalar>
R3Element<Scalar> R3Element<Scalar>::ReadFromMessage(
    serialization::R3Element const& message) {
  return {R3ElementSerializer<Scalar>::ReadFromMessage(message.x()),
          R3ElementSerializer<Scalar>::ReadFromMessage(message.y()),
          R3ElementSerializer<Scalar>::ReadFromMessage(message.z())};
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
