
#pragma once

#include <string>

#include "base/not_null.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace internal_sign {

using base::not_null;
using quantities::Quantity;

// An element of the multiplicative group ({+1, -1}, *). Useful for instance to
// represent the determinant of an orthogonal map.
class Sign final {
 public:
  explicit Sign(double x);
  template<typename Dimensions>
  explicit Sign(Quantity<Dimensions> const& x);
  // The deleted constructor forbids construction from an integer via
  // integer-to-double conversion.  Integers have no signed 0, so this could
  // lead to confusing behaviour.
  template<typename T>
  explicit Sign(T x) = delete;

  template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
  static constexpr Sign OfNonZero(T x);

  constexpr bool is_negative() const;
  constexpr bool is_positive() const;

  constexpr bool operator==(Sign const& other) const;
  constexpr bool operator!=(Sign const& other) const;

  void WriteToMessage(not_null<serialization::Sign*> message) const;
  static Sign ReadFromMessage(serialization::Sign const& message);

 private:
  constexpr explicit Sign(bool negative);

  bool negative_;
  friend Sign operator*(Sign const& left, Sign const& right);
  template<typename T>
  friend T operator*(Sign const& left, T const& right);
};

Sign operator*(Sign const& left, Sign const& right);

// This operator is applicable to any type that has a unary minus operator.
template<typename T>
T operator*(Sign const& left, T const& right);

std::string DebugString(Sign const& sign);

std::ostream& operator<<(std::ostream& out, Sign const& sign);

}  // namespace internal_sign

using internal_sign::Sign;

}  // namespace geometry
}  // namespace principia

#include "geometry/sign_body.hpp"
