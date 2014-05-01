#pragma once

namespace Principia {
namespace Geometry {

class Sign {
 public:
  explicit Sign(const bool positive);
  ~Sign();

  inline bool Negative() const;
  inline bool Positive() const;

 private:
  const bool negative_;
  friend Sign operator*(const Sign& left, const Sign& right);
};

inline Sign operator*(const Sign& left, const Sign& right);

}  // namespace Geometry
}  // namespace Principia

#include "Sign-body.hpp"