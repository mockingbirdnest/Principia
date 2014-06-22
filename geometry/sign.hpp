#pragma once

namespace principia {
namespace geometry {

// An element of the multiplicative group ({+1, -1}, *). Useful for instance to
// represent the determinant of an orthogonal map.
class Sign {
 public:
  template<typename Scalar> explicit Sign(Scalar const& s);
  ~Sign() = default;

  bool Negative() const;
  bool Positive() const;

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
