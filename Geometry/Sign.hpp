#pragma once

namespace principia {
namespace geometry {

// A simple abstraction for something that can take the values -1 and 1.  Useful
// for instance to represent the determinant of an orthogonal map.
class Sign {
 public:
  template <typename Scalar> explicit Sign(const Scalar& s);
  ~Sign() = default;

  inline bool Negative() const;
  inline bool Positive() const;

 private:
  const bool negative_;
  friend Sign operator*(const Sign& left, const Sign& right);
};

inline Sign operator*(const Sign& left, const Sign& right);

}  // namespace geometry
}  // namespace principia

#include "Geometry/Sign-body.hpp"