#pragma once

#include <pmmintrin.h>

#include <cstdint>

namespace principia {
namespace numerics {
namespace _m128d {
namespace internal {

class M128D {
 public:
  M128D(double value);
  explicit M128D(std::int64_t value);

  explicit operator double() const;

  friend M128D operator+(M128D left, M128D right);
  friend M128D operator-(M128D left, M128D right);
  friend M128D operator*(M128D left, M128D right);
  friend M128D operator/(M128D left, M128D right);

  friend M128D operator&(M128D left, M128D right);
  friend M128D operator|(M128D left, M128D right);
  friend M128D operator^(M128D left, M128D right);

  friend M128D FusedMultiplyAdd(M128D a, M128D b, M128D c);
  friend M128D FusedMultiplySubtract(M128D a, M128D b, M128D c);
  friend M128D FusedNegatedMultiplyAdd(M128D a, M128D b, M128D c);
  friend M128D FusedNegatedMultiplySubtract(M128D a, M128D b, M128D c);

 private:
  explicit M128D(__m128d value);

  __m128d value_;
};

M128D operator+(M128D left, M128D right);
M128D operator-(M128D left, M128D right);
M128D operator*(M128D left, M128D right);
M128D operator/(M128D left, M128D right);

M128D operator&(M128D left, M128D right);
M128D operator|(M128D left, M128D right);
M128D operator^(M128D left, M128D right);

// ⟦ab + c⟧.
M128D FusedMultiplyAdd(M128D a, M128D b, M128D c);

// ⟦ab - c⟧.
M128D FusedMultiplySubtract(M128D a, M128D b, M128D c);

// ⟦-ab + c⟧.
M128D FusedNegatedMultiplyAdd(M128D a, M128D b, M128D c);

// ⟦-ab - c⟧.
M128D FusedNegatedMultiplySubtract(M128D a, M128D b, M128D c);

}  // namespace internal

using internal::M128D;

}  // namespace _m128d
}  // namespace numerics
}  // namespace principia

#include "numerics/m128d_body.hpp"
