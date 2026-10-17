#pragma once

#include "numerics/wide.hpp"

#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _wide {
namespace internal {

using namespace principia::quantities::_quantities;

#if PRINCIPIA_USE_AVX()
#define MM(intrinsic) _mm256_ ## intrinsic
#else
#define MM(intrinsic) _mm_ ## intrinsic
#endif

#if PRINCIPIA_USE_AVX()
template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar>::Wide(Scalar const narrow) : Wide(ToM256D(narrow)) {}
#else
template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar>::Wide(Scalar const narrow) : Wide(ToM128D(narrow)) {}
#endif

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar> Wide<Scalar>::Identity() requires ring<Scalar> {
  return Wide(1.0);
}

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Scalar Wide<Scalar>::narrow() const {
#if PRINCIPIA_USE_AVX()
  return FromM256D<Scalar>(value_);
#else
  return FromM128D<Scalar>(value_);
#endif
}

template<additive_group T>
Wide<T>& operator+=(Wide<T>& left, Wide<T> const right) {
  return left = left + right;
}

template<additive_group T>
Wide<T>& operator-=(Wide<T>& left, Wide<T> const right) {
  return left = left - right;
}

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar>::Wide(UnderlyingType const value) : value_(value) {}

inline UnderlyingType wide_sign_bit = MM(set1_pd)(-0.0);

template<additive_group T>
Wide<T> operator+(Wide<T> const right) {
  return right;
}

template<additive_group T>
Wide<T> operator-(Wide<T> const right) {
  return Wide<T>(MM(xor_pd)(right.value_, wide_sign_bit));
}

template<affine T, affine U>
Wide<Sum<T, U>> operator+(Wide<T> const left, Wide<U> const right) {
  return Wide<T>(MM(add_pd)(left.value_, right.value_));
}

template<affine T, affine U>
Wide<Difference<T, U>> operator-(Wide<T> const left, Wide<U> const right) {
  return Wide<Difference<T, U>>(MM(sub_pd)(left.value_, right.value_));
}

template<typename T, typename U>
Wide<Product<T, U>> operator*(Wide<T> const left, Wide<U> const right) {
  return Wide<Product<T, U>>(MM(mul_pd)(left.value_, right.value_));
}

template<typename T, typename U>
Wide<Quotient<T, U>> operator/(Wide<T> const left, Wide<U> const right) {
  return Wide<Quotient<T, U>>(MM(div_pd)(left.value_, right.value_));
}

template<typename T, typename U>
R3Element<Product<T, U>> operator*(R3Element<T> const& left,
                                   Wide<U> const right) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_mul_pd(left.xyzt, right.value_));
#else
  return R3Element<Product<T, U>>(_mm_mul_pd(left.xy, right.value_),
                                  _mm_mul_sd(left.zt, right.value_));
#endif
}

template<typename T, typename U>
R3Element<Product<T, U>> operator*(Wide<T> const left,
                                   R3Element<U> const& right) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_mul_pd(left.value_, right.xyzt));
#else
  return R3Element<Product<T, U>>(_mm_mul_pd(left.value_, right.xy),
                                  _mm_mul_sd(left.value_, right.zt));
#endif
}

template<typename T, typename U>
R3Element<Quotient<T, U>> operator/(R3Element<T> const& left,
                                    Wide<U> const right) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_div_pd(left.xyzt, right.value_));
#else
  return R3Element<Product<T, U>>(_mm_div_pd(left.xy, right.value_),
                                  _mm_div_sd(left.zt, right.value_));
#endif
}

template<typename T>
R3Element<T>& operator*=(R3Element<T>& left, Wide<double> const right) {
#if PRINCIPIA_USE_AVX()
  left.xyzt = _mm256_mul_pd(left.xyzt, right.value_);
#else
  left.xy = _mm_mul_pd(left.xy, right.value_);
  left.zt = _mm_mul_sd(left.zt, right.value_);
#endif
  return left;
}

template<typename T>
R3Element<T>& operator/=(R3Element<T>& left, Wide<double> const right) {
#if PRINCIPIA_USE_AVX()
  left.xyzt = _mm256_div_pd(left.xyzt, right.value_);
#else
  left.xy = _mm_div_pd(left.xy, right.value_);
  left.zt = _mm_div_sd(left.zt, right.value_);
#endif
  return left;
}

template<r3_element_coordinates T, typename U>
Product<T, U> operator*(T const& left, Wide<U> const right) {
  return Product<T, U>(left.coordinates() * right);
}

template<typename T, r3_element_coordinates U>
Product<T, U> operator*(Wide<T> const left, U const& right) {
  return Product<T, U>(left * right.coordinates());
}

template<r3_element_coordinates T, typename U>
Quotient<T, U> operator/(T const& left, Wide<U> const right) {
  return Quotient<T, U>(left.coordinates() / right);
}

template<r3_element_coordinates T>
T& operator*=(T& left, Wide<double> const right) {
  return left = T(left.coordinates() * right);
}

template<r3_element_coordinates T>
T& operator/=(T& left, Wide<double> const right) {
  return left = T(left.coordinates() / right);
}

#undef MM

}  // namespace internal
}  // namespace _wide
}  // namespace numerics
}  // namespace principia
