#pragma once

#include <immintrin.h>

#include "numerics/m128d.hpp"

namespace principia {
namespace numerics {
namespace _m128d {
namespace internal {

M128D const all_ones(0xffff'ffff'ffff'ffffull);
M128D const negated_sign_bit(0x7fff'ffff'ffff'ffffull);
M128D const sign_bit(0x8000'0000'0000'0000ull);

inline M128D::M128D(double const lo) : value_(_mm_set_sd(lo)) {}

inline M128D::M128D(std::int64_t const value)
    : value_(_mm_castsi128_pd(_mm_cvtsi64_si128(value))) {}

inline M128D::M128D(std::uint64_t const value)
    : value_(_mm_castsi128_pd(_mm_cvtsi64_si128(value))) {}

inline M128D::operator double() const {
  return _mm_cvtsd_f64(value_);
}

inline M128D::operator std::int64_t() const {
  return _mm_cvtsi128_si64(_mm_castpd_si128(value_));
}

inline M128D::operator std::uint64_t() const {
  return _mm_cvtsi128_si64(_mm_castpd_si128(value_));
}

inline M128D::operator __m128d() const {
  return value_;
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

inline M128D::M128D(__m128d const value) : value_(value) {}

inline M128D operator+(M128D const right) {
  return right;
}

inline M128D operator-(M128D const right) {
  return M128D(_mm_xor_pd(right.value_, sign_bit.value_));
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

inline M128D operator~(M128D const right) {
  return M128D(_mm_xor_pd(right.value_, all_ones.value_));
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
  return _mm_comieq_sd(left.value_, right.value_);
}

inline bool operator!=(M128D const left, M128D const right) {
  return _mm_comineq_sd(left.value_, right.value_);
}

inline bool operator<(M128D const left, M128D const right) {
  return _mm_comilt_sd(left.value_, right.value_);
}

inline bool operator<=(M128D const left, M128D const right) {
  return _mm_comile_sd(left.value_, right.value_);
}

inline bool operator>=(M128D const left, M128D const right) {
  return _mm_comige_sd(left.value_, right.value_);
}

inline bool operator>(M128D const left, M128D const right) {
  return _mm_comigt_sd(left.value_, right.value_);
}

inline M128D Abs(M128D const a) {
  return a & negated_sign_bit;
}

inline M128D Sign(M128D const a) {
  return a & sign_bit;
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

}  // namespace internal
}  // namespace _m128d
}  // namespace numerics
}  // namespace principia
