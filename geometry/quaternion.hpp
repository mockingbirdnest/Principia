#pragma once

#include "base/not_null.hpp"
#include "geometry/r3_element.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace _quaternion {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_r3_element;

// An element of the skew field of quaternions ℍ (where ℝ is modeled by
// `double`).
class Quaternion final {
 public:
  constexpr Quaternion() = default;
  explicit Quaternion(double real_part);
  Quaternion(double real_part, R3Element<double> const& imaginary_part);

  friend bool operator==(Quaternion const& left,
                         Quaternion const& right) = default;
  friend bool operator!=(Quaternion const& left,
                         Quaternion const& right) = default;

  Quaternion& operator+=(Quaternion const& right);
  Quaternion& operator-=(Quaternion const& right);
  Quaternion& operator*=(Quaternion const& right);
  Quaternion& operator/=(Quaternion const& right);

  Quaternion& operator*=(double right);
  Quaternion& operator/=(double right);

  double real_part() const;
  R3Element<double> const& imaginary_part() const;

  double Norm() const;
  double Norm²() const;

  Quaternion Conjugate() const;
  Quaternion Inverse() const;

  void WriteToMessage(not_null<serialization::Quaternion*> message) const;
  static Quaternion ReadFromMessage(serialization::Quaternion const& message);

 private:
  double real_part_ = 0;
  R3Element<double> imaginary_part_;
};

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

}  // namespace internal

using internal::Normalize;
using internal::Quaternion;

}  // namespace _quaternion
}  // namespace geometry
}  // namespace principia

#include "geometry/quaternion_body.hpp"
