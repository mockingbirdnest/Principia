#pragma once

#include "geometry/r3_element.hpp"

namespace principia {
namespace geometry {

class Quaternion {
 public:
  Quaternion();
  explicit Quaternion(double const real_part);
  Quaternion(double const real_part,
             R3Element<double> const& imaginary_part);

  double const& real_part() const;
  R3Element<double> const& imaginary_part() const;

  Quaternion Conjugate() const;
  Quaternion Inverse() const;

  Quaternion& operator+=(Quaternion const& right);
  Quaternion& operator-=(Quaternion const& right);
  Quaternion& operator*=(Quaternion const& right);
  Quaternion& operator/=(Quaternion const& right);

  Quaternion& operator*=(double const right);
  Quaternion& operator/=(double const right);

 private:
  double real_part_;
  R3Element<double> imaginary_part_;
};

Quaternion operator+(Quaternion const& right);
Quaternion operator-(Quaternion const& right);

Quaternion operator+(Quaternion const& left, Quaternion const& right);
Quaternion operator-(Quaternion const& left, Quaternion const& right);
Quaternion operator*(Quaternion const& left, Quaternion const& right);
Quaternion operator/(Quaternion const& left, Quaternion const& right);

Quaternion operator*(double const left,
                     Quaternion const& right);
Quaternion operator*(Quaternion const& left,
                     double const right);
Quaternion operator/(Quaternion const& left,
                     double const right);

}  // namespace geometry
}  // namespace principia

#include "geometry/quaternion_body.hpp"
