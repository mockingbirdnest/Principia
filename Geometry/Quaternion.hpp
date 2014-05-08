#pragma once

#include "Geometry/R3Element.hpp"
#include "Quantities/Dimensionless.hpp"

namespace principia {
namespace geometry {

class Quaternion {
 public:
  Quaternion();
  Quaternion(quantities::Dimensionless const& real_part,
             R3Element<quantities::Dimensionless> const& imaginary_part);

  quantities::Dimensionless const& real_part() const;
  R3Element<quantities::Dimensionless> const& imaginary_part() const;

 private:
  quantities::Dimensionless const& real_part_;
  R3Element<quantities::Dimensionless> const& imaginary_part_;
};

Quaternion operator+(Quaternion const& right);
Quaternion operator-(Quaternion const& right);

Quaternion operator+(Quaternion const& left,
                            Quaternion const& right);
Quaternion operator-(Quaternion const& left,
                            Quaternion const& right);

Quaternion operator*(quantities::Dimensionless const& left,
                            Quaternion const& right);
Quaternion operator*(Quaternion const& left,
                            quantities::Dimensionless const& right);
Quaternion operator/(Quaternion const& left,
                            quantities::Dimensionless const& right);

void operator+=(Quaternion& left,
                Quaternion const& right);
void operator-=(Quaternion& left,
                Quaternion const& right);

void operator*=(Quaternion& left,
                quantities::Dimensionless const& right);
void operator/=(Quaternion& left,
                quantities::Dimensionless const& right);

}  // namespace geometry
}  // namespace principia

#include "Geometry/Quaternion-body.hpp"
