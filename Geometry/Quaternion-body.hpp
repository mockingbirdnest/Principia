#pragma once

#include "Geometry/Quaternion.hpp"
#include "Geometry/R3Element.hpp"
#include "Quantities/Dimensionless.hpp"

namespace principia {
namespace geometry {

inline Quaternion::Quaternion(quantities::Dimensionless const& real_part)
    : real_part_(real_part) {}

inline Quaternion::Quaternion(
    quantities::Dimensionless const& real_part,
    R3Element<quantities::Dimensionless> const& imaginary_part)
    : real_part_(real_part),
      imaginary_part_(imaginary_part) {}

inline quantities::Dimensionless const& Quaternion::real_part() const {
  return real_part_;
}

inline R3Element<quantities::Dimensionless> const& 
Quaternion::imaginary_part() const {
  return imaginary_part_;
}

inline Quaternion Quaternion::Conjugate() const {
  return Quaternion(real_part_, -imaginary_part_);
}

inline Quaternion Quaternion::Inverse() const {
  return Conjugate() /
      (real_part_ * real_part_ + Dot(imaginary_part_, imaginary_part_));
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

inline Quaternion operator*(quantities::Dimensionless const& left,
                            Quaternion const& right) {
  return Quaternion(left * right.real_part(),
                    left * right.imaginary_part());
}

inline Quaternion operator*(Quaternion const& left,
                            quantities::Dimensionless const& right) {
  return Quaternion(left.real_part() * right,
                    left.imaginary_part() * right);
}

inline Quaternion operator/(Quaternion const& left,
                            quantities::Dimensionless const& right) {
  return Quaternion(left.real_part() / right,
                    left.imaginary_part() / right);
}

inline void operator+=(Quaternion& left, Quaternion const& right) {
  left.real_part_ += right.real_part_;
  left.imaginary_part_ += right.imaginary_part_;
}

inline void operator-=(Quaternion& left, Quaternion const& right)  {
  left.real_part_ -= right.real_part_;
  left.imaginary_part_ -= right.imaginary_part_;
}

inline void operator*=(Quaternion& left, Quaternion const& right) {
  // TODO(phl): Can this be optimized?
  left = left * right;
}

inline void operator/=(Quaternion& left, Quaternion const& right) {
  // TODO(phl): Can this be optimized?
  left = left / right;
}

inline void operator*=(Quaternion& left, 
                       quantities::Dimensionless const& right) {
  left.real_part_ *= right;
  left.imaginary_part_ *= right;
}

inline void operator/=(Quaternion& left, 
                       quantities::Dimensionless const& right) {
  left.real_part_ /= right;
  left.imaginary_part_ /= right;
}

}  // namespace geometry
}  // namespace principia
