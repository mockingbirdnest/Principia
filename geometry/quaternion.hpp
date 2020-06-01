
#pragma once

#include "geometry/r3_element.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace internal_quaternion {

using base::not_null;

// An element of the skew field of quaternions ℍ (where ℝ is modeled by
// |double|).
class Quaternion final {
 public:
  constexpr Quaternion() = default;
  explicit Quaternion(double real_part);
  Quaternion(double real_part, R3Element<double> const& imaginary_part);

  double real_part() const;
  R3Element<double> const& imaginary_part() const;

  double Norm() const;
  double Norm²() const;

  Quaternion Conjugate() const;
  Quaternion Inverse() const;

  Quaternion& operator+=(Quaternion const& right);
  Quaternion& operator-=(Quaternion const& right);
  Quaternion& operator*=(Quaternion const& right);
  Quaternion& operator/=(Quaternion const& right);

  Quaternion& operator*=(double right);
  Quaternion& operator/=(double right);

  void WriteToMessage(not_null<serialization::Quaternion*> message) const;
  static Quaternion ReadFromMessage(serialization::Quaternion const& message);

 private:
  double real_part_ = 0;
  R3Element<double> imaginary_part_;
};

bool operator==(Quaternion const& left, Quaternion const& right);
bool operator!=(Quaternion const& left, Quaternion const& right);

Quaternion operator+(Quaternion const& right);
Quaternion operator-(Quaternion const& right);

Quaternion operator+(Quaternion const& left, Quaternion const& right);
Quaternion operator-(Quaternion const& left, Quaternion const& right);
Quaternion operator*(Quaternion const& left, Quaternion const& right);
Quaternion operator/(Quaternion const& left, Quaternion const& right);

Quaternion operator*(double left, Quaternion const& right);
Quaternion operator*(Quaternion const& left, double right);
Quaternion operator/(Quaternion const& left, double right);

Quaternion Normalize(Quaternion const& quaternion);

std::ostream& operator<<(std::ostream& out, Quaternion const& quaternion);

}  // namespace internal_quaternion

using internal_quaternion::Normalize;
using internal_quaternion::Quaternion;

}  // namespace geometry
}  // namespace principia

#include "geometry/quaternion_body.hpp"
