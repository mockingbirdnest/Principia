
#pragma once

#include <string>

#include "base/not_null.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace internal_sign {

using base::not_null;

// An element of the multiplicative group ({+1, -1}, *). Useful for instance to
// represent the determinant of an orthogonal map.
class Sign final {
 public:
  template<typename Scalar> explicit Sign(Scalar const& s);

  bool Negative() const;
  bool Positive() const;

  bool operator==(Sign const& other) const;
  bool operator!=(Sign const& other) const;

  void WriteToMessage(not_null<serialization::Sign*> const message) const;
  static Sign ReadFromMessage(serialization::Sign const& message);

 private:
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
