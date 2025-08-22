#pragma once

#include "geometry/quaternion.hpp"

#include <string>

#include "geometry/sign.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace _quaternion {
namespace internal {

using namespace principia::geometry::_sign;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_quantities;

inline Quaternion::Quaternion(double const real_part)
    : real_part_(real_part) {}

inline Quaternion::Quaternion(double const real_part,
                              R3Element<double> const& imaginary_part)
    : real_part_(real_part),
      imaginary_part_(imaginary_part) {}

inline Quaternion& Quaternion::operator+=(Quaternion const& right) {
  real_part_ += right.real_part_;
  imaginary_part_ += right.imaginary_part_;
  return *this;
}

inline Quaternion& Quaternion::operator-=(Quaternion const& right)  {
  real_part_ -= right.real_part_;
  imaginary_part_ -= right.imaginary_part_;
  return *this;
}

inline Quaternion& Quaternion::operator*=(Quaternion const& right) {
  // TODO(phl): Can this be optimized?
  return *this = *this * right;
}

inline Quaternion& Quaternion::operator/=(Quaternion const& right) {
  // TODO(phl): Can this be optimized?
  return *this = *this / right;
}

inline Quaternion& Quaternion::operator*=(double const right) {
  real_part_ *= right;
  imaginary_part_ *= right;
  return *this;
}

inline Quaternion& Quaternion::operator/=(double const right) {
  real_part_ /= right;
  imaginary_part_ /= right;
  return *this;
}

inline double Quaternion::real_part() const {
  return real_part_;
}

inline R3Element<double> const&
Quaternion::imaginary_part() const {
  return imaginary_part_;
}

inline double Quaternion::Norm() const {
  return Sqrt(Norm²());
}

inline double Quaternion::Norm²() const {
  return real_part_ * real_part_ + imaginary_part_.Norm²();
}

inline Quaternion Quaternion::Conjugate() const {
  return Quaternion(real_part_, -imaginary_part_);
}

inline Quaternion Quaternion::Inverse() const {
  return Conjugate() /
      (real_part_ * real_part_ + Dot(imaginary_part_, imaginary_part_));
}

inline void Quaternion::WriteToMessage(
    not_null<serialization::Quaternion*> const message) const {
  message->set_real_part(real_part_);
  imaginary_part_.WriteToMessage(message->mutable_imaginary_part());
}

inline Quaternion Quaternion::ReadFromMessage(
    serialization::Quaternion const& message) {
  return Quaternion(message.real_part(),
                    R3Element<double>::ReadFromMessage(
                        message.imaginary_part()));
}

inline Quaternion operator+(Quaternion const& right) {
  return right;
}

inline Quaternion operator-(Quaternion const& right) {
  return Quaternion(-right.real_part(), -right.imaginary_part());
}

inline Quaternion operator+(Quaternion const& left, Quaternion const& right) {
  return Quaternion(left.real_part() + right.real_part(),
                    left.imaginary_part() + right.imaginary_part());
}

inline Quaternion operator-(Quaternion const& left, Quaternion const& right) {
  return Quaternion(left.real_part() - right.real_part(),
                    left.imaginary_part() - right.imaginary_part());
}

inline Quaternion operator*(Quaternion const& left, Quaternion const& right) {
  return Quaternion(left.real_part() * right.real_part() -
                        Dot(left.imaginary_part(), right.imaginary_part()),
                    left.real_part() * right.imaginary_part() +
                        right.real_part() * left.imaginary_part() +
                        Cross(left.imaginary_part(), right.imaginary_part()));
}

inline Quaternion operator/(Quaternion const& left, Quaternion const& right) {
  return left * right.Inverse();
}

inline Quaternion operator*(double const left, Quaternion const& right) {
  return Quaternion(left * right.real_part(),
                    left * right.imaginary_part());
}

inline Quaternion operator*(Quaternion const& left, double const right) {
  return Quaternion(left.real_part() * right,
                    left.imaginary_part() * right);
}

inline Quaternion operator/(Quaternion const& left, double const right) {
  return Quaternion(left.real_part() / right,
                    left.imaginary_part() / right);
}

inline Quaternion Normalize(Quaternion const& quaternion) {
  return quaternion / quaternion.Norm();
}

inline std::ostream& operator<<(std::ostream& out,
                                Quaternion const& quaternion) {
  double const w = quaternion.real_part();
  double const x = quaternion.imaginary_part().x;
  double const y = quaternion.imaginary_part().y;
  double const z = quaternion.imaginary_part().z;

  auto const unsigned_debug_string = [](double const d) {
    std::string const s = DebugString(Abs(d));
    CHECK_EQ('+', s[0]);
    return s.substr(1);
  };

  return out << w << " "
             << Sign(x) << " " << unsigned_debug_string(x) << " i "
             << Sign(y) << " " << unsigned_debug_string(y) << " j "
             << Sign(z) << " " << unsigned_debug_string(z) << " k";
}

}  // namespace internal
}  // namespace _quaternion
}  // namespace geometry
}  // namespace principia
