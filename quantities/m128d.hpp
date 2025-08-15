#pragma once

#include <pmmintrin.h>

#include <cstdint>
#include <ostream>

namespace principia {
namespace quantities {
namespace _m128d {
namespace internal {

class M128D {
 public:
  M128D() = default;
  explicit M128D(double value);
  explicit M128D(std::int64_t value);
  explicit M128D(std::uint64_t value);

  M128D(M128D const volatile& v);
  M128D(M128D const& v) = default;
  M128D& operator=(M128D const&) = default;

  explicit operator double() const;
  explicit operator std::int64_t() const;
  explicit operator std::uint64_t() const;
  explicit operator __m128d() const;

  M128D& operator+=(M128D right);
  M128D& operator-=(M128D right);
  M128D& operator*=(M128D right);
  template<std::integral T>
  M128D& operator*=(T right);
  M128D& operator/=(M128D right);

  friend M128D operator+(M128D right);
  friend M128D operator-(M128D right);
  friend M128D operator+(M128D left, M128D right);
  friend M128D operator-(M128D left, M128D right);
  friend M128D operator*(M128D left, M128D right);
  template<std::integral T>
  friend M128D operator*(M128D left, T right);
  template<std::integral T>
  friend M128D operator*(T left, M128D right);
  friend M128D operator*(M128D left, double right);
  friend M128D operator*(double left, M128D right);
  friend M128D operator/(M128D left, M128D right);

  friend M128D operator~(M128D right);
  friend M128D operator&(M128D left, M128D right);
  friend M128D operator|(M128D left, M128D right);
  friend M128D operator^(M128D left, M128D right);

  friend bool operator==(M128D left, M128D right);
  friend bool operator!=(M128D left, M128D right);
  friend bool operator<(M128D left, M128D right);
  friend bool operator<(M128D left, double right);
  friend bool operator<=(M128D left, M128D right);
  friend bool operator<=(M128D left, double right);
  friend bool operator>=(M128D left, M128D right);
  friend bool operator>=(M128D left, double right);
  friend bool operator>(M128D left, M128D right);

  friend M128D Abs(M128D a);
  friend M128D Sign(M128D a);

  friend M128D FusedMultiplyAdd(M128D a, M128D b, M128D c);
  friend M128D FusedMultiplySubtract(M128D a, M128D b, M128D c);
  friend M128D FusedNegatedMultiplyAdd(M128D a, M128D b, M128D c);
  friend M128D FusedNegatedMultiplySubtract(M128D a, M128D b, M128D c);

 private:
  explicit M128D(__m128d value);

  __m128d value_;
};

M128D operator+(M128D right);
M128D operator-(M128D right);
M128D operator+(M128D left, M128D right);
M128D operator-(M128D left, M128D right);
M128D operator*(M128D left, M128D right);
template<std::integral T>
M128D operator*(M128D left, T right);
template<std::integral T>
M128D operator*(T left, M128D right);
M128D operator/(M128D left, M128D right);

M128D operator~(M128D right);
M128D operator&(M128D left, M128D right);
M128D operator|(M128D left, M128D right);
M128D operator^(M128D left, M128D right);

bool operator==(M128D left, M128D right);
bool operator!=(M128D left, M128D right);
bool operator<(M128D left, M128D right);
bool operator<(M128D left, double right);
bool operator<=(M128D left, M128D right);
bool operator<=(M128D left, double right);
bool operator>=(M128D left, M128D right);
bool operator>=(M128D left, double right);
bool operator>(M128D left, M128D right);

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
