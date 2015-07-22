#pragma once

#include "base/not_null.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using base::not_null;

namespace geometry {

// An element of the multiplicative group ({+1, -1}, *). Useful for instance to
// represent the determinant of an orthogonal map.
class Sign {
 public:
  template<typename Scalar> explicit Sign(Scalar const& s);
  ~Sign() = default;

  bool Negative() const;
  bool Positive() const;

  void WriteToMessage(not_null<serialization::Sign*> const message) const;
  static Sign ReadFromMessage(serialization::Sign const& message);

  bool operator==(Sign const& right) const;
  bool operator!=(Sign const& right) const;

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

}  // namespace geometry
}  // namespace principia

#include "geometry/sign_body.hpp"
