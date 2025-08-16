#pragma once

#include <pmmintrin.h>

#include <concepts>
#include <cstdint>
#include <ostream>

namespace principia {
namespace quantities {
namespace _m128d {
namespace internal {

class M128D {
 public:
  M128D() = default;

  // No constructors taking integral types, because the integral-to-floating
  // conversions are costly.  The client should write `M128D(2.0)`, not
  // `M128D(2)`.
  template<std::floating_point T>
  explicit M128D(T value);
  explicit M128D(__m128d value);

  // No conversion to integral types because the floating-to-integral conversion
  // is costly and should be explicit in client code.
  explicit operator double() const;
  explicit operator __m128d() const;

  // These functions are the equivalent of `reinterpret_cast`, they just copy
  // the bits without any integral/floating conversion.
  template<std::integral T>
  static M128D MakeFromBits(T value);
  template<std::integral T>
  T Bits() const;

  M128D& operator+=(M128D right);
  M128D& operator-=(M128D right);
  M128D& operator*=(M128D right);
  M128D& operator/=(M128D right);

  friend M128D operator+(M128D right);
  friend M128D operator-(M128D right);
  friend M128D operator+(M128D left, M128D right);
  friend M128D operator-(M128D left, M128D right);
  friend M128D operator*(M128D left, M128D right);
  friend M128D operator*(M128D left, double right);//TODO(phl)remove
  friend M128D operator*(double left, M128D right);
  friend M128D operator/(M128D left, M128D right);

  // The ℤ-module structure.  It is important to use `std::integral` here to
  // make sure that these operations are not callable with an implicitly-
  // converted `double`.
  template<std::integral T>
  M128D& operator*=(T right);
  template<std::integral T>
  friend M128D operator*(M128D left, T right);
  template<std::integral T>
  friend M128D operator*(T left, M128D right);

  friend M128D operator~(M128D right);
  friend M128D operator&(M128D left, M128D right);
  friend M128D operator|(M128D left, M128D right);
  friend M128D operator^(M128D left, M128D right);

  // Comparisons should be implemented as vanilla float comparisons, not as
  // calls to intrinsics like `_mm_comieq_sd`.  The former generates an
  // instruction like `ucomisd` and does a conditional branch, while the latter
  // produces `setcc` instruction to build an integer from the flags, and then
  // does a `test`.  For this reason, we need overloads that take `double` (we
  // are not going to build an `M128D` and immediately tear it apart).
  // This also helps with OSACA constexpr-ness.
  friend bool operator==(M128D left, M128D right);
  friend bool operator==(M128D left, double right);
  friend bool operator==(double left, M128D right);
  friend bool operator!=(M128D left, M128D right);
  friend bool operator!=(M128D left, double right);
  friend bool operator!=(double left, M128D right);
  friend bool operator<(M128D left, M128D right);
  friend bool operator<(M128D left, double right);
  friend bool operator<(double left, M128D right);
  friend bool operator<=(M128D left, M128D right);
  friend bool operator<=(M128D left, double right);
  friend bool operator<=(double left, M128D right);
  friend bool operator>=(M128D left, M128D right);
  friend bool operator>=(M128D left, double right);
  friend bool operator>=(double left, M128D right);
  friend bool operator>=(M128D left, double right);
  friend bool operator>(M128D left, M128D right);
  friend bool operator>(M128D left, double right);
  friend bool operator>(double left, M128D right);

  friend M128D Abs(M128D a);
  friend M128D Sign(M128D a);

  friend M128D FusedMultiplyAdd(M128D a, M128D b, M128D c);
  friend M128D FusedMultiplySubtract(M128D a, M128D b, M128D c);
  friend M128D FusedNegatedMultiplyAdd(M128D a, M128D b, M128D c);
  friend M128D FusedNegatedMultiplySubtract(M128D a, M128D b, M128D c);

 private:
  __m128d value_;

  static M128D const all_ones_;
  static M128D const negated_sign_bit_;
  static M128D const sign_bit_;
};

M128D operator+(M128D right);
M128D operator-(M128D right);
M128D operator+(M128D left, M128D right);
M128D operator-(M128D left, M128D right);
M128D operator*(M128D left, M128D right);
M128D operator/(M128D left, M128D right);

template<std::integral T>
M128D operator*(M128D left, T right);
template<std::integral T>
M128D operator*(T left, M128D right);

M128D operator~(M128D right);
M128D operator&(M128D left, M128D right);
M128D operator|(M128D left, M128D right);
M128D operator^(M128D left, M128D right);

bool operator==(M128D left, M128D right);
bool operator==(M128D left, double right);
bool operator==(double left, M128D right);
bool operator!=(M128D left, M128D right);
bool operator!=(M128D left, double right);
bool operator!=(double left, M128D right);
bool operator<(M128D left, M128D right);
bool operator<(M128D left, double right);
bool operator<(double left, M128D right);
bool operator<=(M128D left, M128D right);
bool operator<=(M128D left, double right);
bool operator<=(double left, M128D right);
bool operator>=(M128D left, M128D right);
bool operator>=(M128D left, double right);
bool operator>=(double left, M128D right);
bool operator>=(M128D left, double right);
bool operator>(M128D left, M128D right);
bool operator>(M128D left, double right);
bool operator>(double left, M128D right);

M128D Abs(M128D a);
M128D Sign(M128D a);

// ⟦ab + c⟧.
M128D FusedMultiplyAdd(M128D a, M128D b, M128D c);

// ⟦ab - c⟧.
M128D FusedMultiplySubtract(M128D a, M128D b, M128D c);

// ⟦-ab + c⟧.
M128D FusedNegatedMultiplyAdd(M128D a, M128D b, M128D c);

// ⟦-ab - c⟧.
M128D FusedNegatedMultiplySubtract(M128D a, M128D b, M128D c);

std::ostream& operator<<(std::ostream& os, M128D x);

}  // namespace internal

using internal::Abs;
using internal::FusedMultiplyAdd;
using internal::FusedMultiplySubtract;
using internal::FusedNegatedMultiplyAdd;
using internal::FusedNegatedMultiplySubtract;
using internal::M128D;
using internal::Sign;

}  // namespace _m128d
}  // namespace quantities
}  // namespace principia

#include "quantities/m128d_body.hpp"
