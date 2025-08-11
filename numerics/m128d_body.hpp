#pragma once

#include "numerics/m128d.hpp"

namespace principia {
namespace numerics {
namespace _m128d {
namespace internal {

inline M128D::M128D(double const lo) : value_(_mm_set_sd(lo)) {}

inline M128D::M128D(std::int64_t const value)
    : value_(_mm_castsi128_pd(_mm_cvtsi64_si128(value))) {}

inline M128D::operator double() const {
  return _mm_cvtsd_f64(value_);
}

inline M128D::M128D(__m128d const value) : value_(value) {}

M128D operator+(M128D const left, M128D const right) {
  return M128D(_mm_add_sd(left.value_, right.value_));
}

M128D operator-(M128D const left, M128D const right) {
  return M128D(_mm_sub_sd(left.value_, right.value_));
}

M128D operator*(M128D const left, M128D const right) {
  return M128D(_mm_mul_sd(left.value_, right.value_));
}

M128D operator/(M128D const left, M128D const right) {
  return M128D(_mm_div_sd(left.value_, right.value_));
}

M128D operator&(M128D const left, M128D const right) {
  return M128D(_mm_and_pd(left.value_, right.value_));
}

M128D operator|(M128D const left, M128D const right) {
  return M128D(_mm_or_pd(left.value_, right.value_));
}

M128D operator^(M128D const left, M128D const right) {
  return M128D(_mm_and_pd(left.value_, right.value_));
}

M128D FusedMultiplyAdd(M128D const a, M128D const b, M128D const c) {
  return M128D(_mm_fmadd_sd(a.value_, b.value_, c.value_));
}

M128D FusedMultiplySubtract(M128D const a, M128D const b, M128D const c) {
  return M128D(_mm_fmsub_sd(a.value_, b.value_, c.value_));
}

M128D FusedNegatedMultiplyAdd(M128D const a, M128D const b, M128D const c) {
  return M128D(_mm_fnmadd_sd(a.value_, b.value_, c.value_));
}

M128D FusedNegatedMultiplySubtract(M128D const a,
                                   M128D const b,
                                   M128D const c) {
  return M128D(_mm_fnmsub_sd(a.value_, b.value_, c.value_));
}

}  // namespace internal
}  // namespace _m128d
}  // namespace numerics
}  // namespace principia
