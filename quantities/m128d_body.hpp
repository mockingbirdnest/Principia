#pragma once

#include <immintrin.h>

#include "quantities/m128d.hpp"

namespace principia {
namespace quantities {
namespace _m128d {
namespace internal {

template<std::floating_point T>
M128D::M128D(T const value) : value_(_mm_set_sd(value)) {}

inline M128D::M128D(__m128d const value) : value_(value) {}

inline M128D::operator double() const {
  return _mm_cvtsd_f64(value_);
}

inline M128D::operator __m128d() const {
  return value_;
}

template<std::integral T>
M128D M128D::MakeFromBits(T const value) {
  return M128D(_mm_castsi128_pd(_mm_cvtsi64_si128(value)));
}

template<std::integral T>
T M128D::Bits() const {
  return _mm_cvtsi128_si64(_mm_castpd_si128(value_));
}

inline M128D& M128D::operator+=(M128D const right) {
  *this = *this + right;
  return *this;
}

inline M128D& M128D::operator-=(M128D const right) {
  *this = *this - right;
  return *this;
}

inline M128D& M128D::operator*=(M128D const right) {
  *this = *this * right;
  return *this;
}

inline M128D& M128D::operator/=(M128D const right) {
  *this = *this / right;
  return *this;
}

inline M128D const M128D::all_ones_ =
    M128D::MakeFromBits(0xffff'ffff'ffff'ffffull);
inline M128D const M128D::negated_sign_bit_ =
    M128D::MakeFromBits(0x7fff'ffff'ffff'ffffull);
inline M128D const M128D::sign_bit_ =
    M128D::MakeFromBits(0x8000'0000'0000'0000ull);

inline M128D operator+(M128D const right) {
  return right;
}

inline M128D operator-(M128D const right) {
  return M128D(_mm_xor_pd(right.value_, M128D::sign_bit_.value_));
}

inline M128D operator+(M128D const left, M128D const right) {
  return M128D(_mm_add_sd(left.value_, right.value_));
}

inline M128D operator-(M128D const left, M128D const right) {
  return M128D(_mm_sub_sd(left.value_, right.value_));
}

inline M128D operator*(M128D const left, M128D const right) {
  return M128D(_mm_mul_sd(left.value_, right.value_));
}

inline M128D operator/(M128D const left, M128D const right) {
  return M128D(_mm_div_sd(left.value_, right.value_));
}

template<std::integral T>
M128D& M128D::operator*=(T const right) {
  *this = *this * M128D(static_cast<double>(right));
  return *this;
}

template<std::integral T>
M128D operator*(M128D const left, T const right) {
  return left * M128D(static_cast<double>(right));
}

template<std::integral T>
M128D operator*(T const left, M128D const right) {
  return M128D(static_cast<double>(left)) * right;
}

inline M128D operator~(M128D const right) {
  return M128D(_mm_xor_pd(right.value_, M128D::all_ones_.value_));
}

inline M128D operator&(M128D const left, M128D const right) {
  return M128D(_mm_and_pd(left.value_, right.value_));
}

inline M128D operator|(M128D const left, M128D const right) {
  return M128D(_mm_or_pd(left.value_, right.value_));
}

inline M128D operator^(M128D const left, M128D const right) {
  return M128D(_mm_xor_pd(left.value_, right.value_));
}

inline bool operator==(M128D const left, M128D const right) {
  return static_cast<double>(left) == static_cast<double>(right);
}

template<std::floating_point T>
bool operator==(M128D const left, T const right) {
  return static_cast<double>(left) == right;
}

template<std::floating_point T>
bool operator==(T const left, M128D const right) {
  return left == static_cast<double>(right);
}

inline bool operator!=(M128D const left, M128D const right) {
  return static_cast<double>(left) != static_cast<double>(right);
}

template<std::floating_point T>
bool operator!=(M128D const left, T const right) {
  return static_cast<double>(left) != right;
}

template<std::floating_point T>
bool operator!=(T const left, M128D const right) {
  return left != static_cast<double>(right);
}

inline bool operator<(M128D const left, M128D const right) {
  return static_cast<double>(left) < static_cast<double>(right);
}

template<std::floating_point T>
bool operator<(M128D const left, T const right) {
  return static_cast<double>(left) < right;
}

template<std::floating_point T>
bool operator<(T const left, M128D const right) {
  return left < static_cast<double>(right);
}

inline bool operator<=(M128D const left, M128D const right) {
  return static_cast<double>(left) <= static_cast<double>(right);
}

template<std::floating_point T>
bool operator<=(M128D const left, T const right) {
  return static_cast<double>(left) <= right;
}

template<std::floating_point T>
bool operator<=(T const left, M128D const right) {
  return left <= static_cast<double>(right);
}

inline bool operator>=(M128D const left, M128D const right) {
  return static_cast<double>(left) >= static_cast<double>(right);
}

template<std::floating_point T>
bool operator>=(M128D const left, T const right) {
  return static_cast<double>(left) >= right;
}

template<std::floating_point T>
bool operator>=(T const left, M128D const right) {
  return left >= static_cast<double>(right);
}

inline bool operator>(M128D const left, M128D const right) {
  return static_cast<double>(left) > static_cast<double>(right);
}

template<std::floating_point T>
bool operator>(M128D const left, T const right) {
  return static_cast<double>(left) > right;
}

template<std::floating_point T>
bool operator>(T const left, M128D const right) {
  return left > static_cast<double>(right);
}

inline M128D Abs(M128D const a) {
  return a & M128D::negated_sign_bit_;
}

inline M128D Sign(M128D const a) {
  return a & M128D::sign_bit_;
}

inline M128D FusedMultiplyAdd(M128D const a, M128D const b, M128D const c) {
  return M128D(_mm_fmadd_sd(a.value_, b.value_, c.value_));
}

inline M128D FusedMultiplySubtract(M128D const a,
                                   M128D const b,
                                   M128D const c) {
  return M128D(_mm_fmsub_sd(a.value_, b.value_, c.value_));
}

inline M128D FusedNegatedMultiplyAdd(M128D const a,
                                     M128D const b,
                                     M128D const c) {
  return M128D(_mm_fnmadd_sd(a.value_, b.value_, c.value_));
}

inline M128D FusedNegatedMultiplySubtract(M128D const a,
                                          M128D const b,
                                          M128D const c) {
  return M128D(_mm_fnmsub_sd(a.value_, b.value_, c.value_));
}

inline std::ostream& operator<<(std::ostream& os, M128D const x) {
  return os << static_cast<double>(x);
}

}  // namespace internal
}  // namespace _m128d
}  // namespace quantities
}  // namespace principia
