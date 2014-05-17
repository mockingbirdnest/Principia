#pragma once

#include "geometry/r3_element.hpp"
#include "quantities/dimensionless.hpp"

namespace principia {
namespace geometry {

class Quaternion {
 public:
  Quaternion() = default;
  explicit Quaternion(quantities::Dimensionless const& real_part);
  Quaternion(quantities::Dimensionless const& real_part,
             R3Element<quantities::Dimensionless> const& imaginary_part);

  quantities::Dimensionless const& real_part() const;
  R3Element<quantities::Dimensionless> const& imaginary_part() const;

  Quaternion Conjugate() const;
  Quaternion Inverse() const;

 private:
  quantities::Dimensionless real_part_;
  R3Element<quantities::Dimensionless> imaginary_part_;

  friend void operator+=(Quaternion& left,
                         Quaternion const& right);
  friend void operator-=(Quaternion& left,
                         Quaternion const& right);
  friend void operator*=(Quaternion& left,
                         quantities::Dimensionless const& right);
  friend void operator/=(Quaternion& left,
                         quantities::Dimensionless const& right);
};

Quaternion operator+(Quaternion const& right);
Quaternion operator-(Quaternion const& right);

Quaternion operator+(Quaternion const& left, Quaternion const& right);
Quaternion operator-(Quaternion const& left, Quaternion const& right);
Quaternion operator*(Quaternion const& left, Quaternion const& right);
Quaternion operator/(Quaternion const& left, Quaternion const& right);

Quaternion operator*(quantities::Dimensionless const& left,
                     Quaternion const& right);
Quaternion operator*(Quaternion const& left,
                     quantities::Dimensionless const& right);
Quaternion operator/(Quaternion const& left,
                     quantities::Dimensionless const& right);

void operator+=(Quaternion& left, Quaternion const& right);
void operator-=(Quaternion& left, Quaternion const& right);
void operator*=(Quaternion& left, Quaternion const& right);
void operator/=(Quaternion& left, Quaternion const& right);

void operator*=(Quaternion& left, quantities::Dimensionless const& right);
void operator/=(Quaternion& left, quantities::Dimensionless const& right);

}  // namespace geometry
}  // namespace principia

#include "geometry/quaternion_body.hpp"
