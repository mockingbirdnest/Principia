#pragma once
#include "numerics/wide.hpp"

#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _wide {
namespace internal {

using namespace principia::quantities::_quantities;

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
inline Wide<Scalar> Wide<Scalar>::Identity() requires ring<Scalar> {
  return Widen(1.0);
}

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Scalar Wide<Scalar>::narrow() const {
  return narrow_value_;
}

#if PRINCIPIA_USE_AVX()
template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar>::Wide(__m256d value) : value_(value) {}
#else
template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar>::Wide(__m256d value) : value_(value) {}
#endif

template<quantity T>
Wide<T> Widen(T narrow) {
#if PRINCIPIA_USE_AVX()
  return Wide<T>(ToM256D(narrow));
#else
  return Wide<T>(ToM128D(narrow));
#endif
}

inline Wide<double> wide_sign_bit = Widen(-0.0);

template<additive_group T>
Wide<T> operator+(Wide<T> const right) {
  return right;
}

template<additive_group T>
Wide<T> operator-(Wide<T> const right) {
#if PRINCIPIA_USE_AVX()
  return Wide<T>(_mm256_xor_pd(right.value_, wide_sign_bit.value_));
#else
  return Wide<T>(_mm_xor_pd(right.value_, wide_sign_bit.value_));
#endif
}

template<affine T, affine U>
Wide<Sum<T, U>> operator+(Wide<T> const left, Wide<U> const right) {
#if PRINCIPIA_USE_AVX()
  return Wide<T>(_mm256_add_pd(left.value_, right.value_));
#else
  return Wide<T>(_mm_add_pd(left.value_, right.value_));
#endif
}

template<affine T, affine U>
Wide<Difference<T, U>> operator-(Wide<T> const left, Wide<U> const right) {
#if PRINCIPIA_USE_AVX()
  return Wide<Difference<T, U>>(_mm256_sub_pd(left.value_, right.value_));
#else
  return Wide<Difference<T, U>>(_mm_sub_pd(left.value_, right.value_));
#endif
}

template<typename T, typename U>
Wide<Product<T, U>> operator*(Wide<T> const left, Wide<U> const right) {
#if PRINCIPIA_USE_AVX()
  return Wide<Product<T, U>>(_mm256_mul_pd(left.value_, right.value_));
#else
  return Wide<Product<T, U>>(_mm_mul_pd(left.value_, right.value_));
#endif
}

template<typename T, typename U>
Wide<Quotient<T, U>> operator/(Wide<T> const left, Wide<U> const right) {
#if PRINCIPIA_USE_AVX()
  return Wide<Quotient<T, U>>(_mm256_div_pd(left.value_, right.value_));
#else
  return Wide<Quotient<T, U>>(_mm_div_pd(left.value_, right.value_));
#endif
}

template<additive_group T>
Wide<T>& operator+=(Wide<T>& left, Wide<T> right) {
  return left = left + right;
}

template<additive_group T>
Wide<T>& operator-=(Wide<T>& left, Wide<T> right) {
  return left = left - right;
}

template<r3_element_coordinates T, typename U>
Product<T, U> operator*(T left, Wide<U> right) {
  return Product<T, U>(left.coordinates() * right);
}

template<typename T, r3_element_coordinates U>
Product<T, U> operator*(Wide<T> left, U right) {
  return Product<T, U>(left * right.coordinates());
}

template<r3_element_coordinates T>
T& operator*=(T& left, Wide<double> right) {
  return left = T(left.coordinates() * right);
}

template<r3_element_coordinates T, typename U>
Quotient<T, U> operator/(T left, Wide<U> right) {
  return Quotient<T, U>(left.coordinates() / right);
}

template<r3_element_coordinates T>
T& operator/=(T& left, Wide<double> right) {
  return left = T(left.coordinates() / right);
}

template<typename T, typename U>
R3Element<Product<T, U>> operator*(R3Element<T> left, Wide<U> right) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_mul_pd(left.xyzt, right.value_));
#else
  return R3Element<Product<T, U>>(_mm_mul_pd(left.xy, right.value_),
                                  _mm_mul_sd(left.zt, right.value_));
#endif
}

template<typename T, typename U>
R3Element<Product<T, U>> operator*(Wide<T> left, R3Element<U> right) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_mul_pd(left.value_, right.xyzt));
#else
  return R3Element<Product<T, U>>(_mm_mul_pd(left.value_, right.xy),
                                  _mm_mul_sd(left.value_, right.zt));
#endif
}

template<typename T>
R3Element<T>& operator*=(R3Element<T>& left, Wide<double> right) {
#if PRINCIPIA_USE_AVX()
  left.xyzt = _mm256_mul_pd(left.xyzt, right.value_);
#else
  left.xy = _mm_mul_pd(left.xy, right.value_);
  left.zt = _mm_mul_sd(left.zt, right.value_);
#endif
  return left;
}

template<typename T, typename U>
R3Element<Quotient<T, U>> operator/(R3Element<T> left, Wide<U> right) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_div_pd(left.xyzt, right.value_));
#else
  return R3Element<Product<T, U>>(_mm_div_pd(left.xy, right.value_),
                                  _mm_div_sd(left.zt, right.value_));
#endif
}

template<typename T>
R3Element<T>& operator/=(R3Element<T>& left, Wide<double> right) {
#if PRINCIPIA_USE_AVX()
  left.xyzt = _mm256_div_pd(left.xyzt, right.value_);
#else
  left.xy = _mm_div_pd(left.xy, right.value_);
  left.zt = _mm_div_sd(left.zt, right.value_);
#endif
  return left;
}

}  // namespace internal
}  // namespace _wide
}  // namespace numerics
}  // namespace principia
