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

inline UnderlyingType wide_sign_bit = MM(set1_pd)(-0.0);

#if PRINCIPIA_USE_AVX()
template<real_affine_space Scalar>
  requires (real_dimension<Scalar> == 1)  // NOLINT
Wide<Scalar>::Wide(Scalar const narrow)
  requires quantity<Scalar>
    : Wide(ToM256D(narrow)) {}

template<real_affine_space Scalar>
  requires (real_dimension<Scalar> == 1)  // NOLINT
Wide<Scalar>::Wide(Scalar const narrow)
  requires (!quantity<Scalar>) && quantity<Difference<Scalar>>  // NOLINT
    : Wide(ToM256D(narrow - Scalar{})) {}
#else
template<real_affine_space Scalar>
  requires (real_dimension<Scalar> == 1)  // NOLINT
Wide<Scalar>::Wide(Scalar const narrow)
  requires quantity<Scalar>
    : Wide(ToM128D(narrow)) {}

template<real_affine_space Scalar>
  requires (real_dimension<Scalar> == 1)  // NOLINT
Wide<Scalar>::Wide(Scalar const narrow)
  requires (!quantity<Scalar>) && quantity<Difference<Scalar>>  // NOLINT
    : Wide(ToM128D(narrow - Scalar{})) {}
#endif

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar> Wide<Scalar>::Identity() requires ring<Scalar> {
  return Wide(1.0);
}

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Scalar Wide<Scalar>::narrow() const {
  if constexpr (quantity<Scalar>) {
#if PRINCIPIA_USE_AVX()
  return FromM256D<Scalar>(value_);
#else
  return FromM128D<Scalar>(value_);
#endif
  } else {
#if PRINCIPIA_USE_AVX()
  return Scalar{} + FromM256D<Difference<Scalar>>(value_);
#else
  return Scalar{} + FromM128D<Difference<Scalar>>(value_);
#endif
  }
}

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar>& Wide<Scalar>::operator+=(Wide<Scalar> const right)
  requires additive_group<Scalar> {
  return *this = *this + right;
}

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar>& Wide<Scalar>::operator-=(Wide<Scalar> const right)
  requires additive_group<Scalar> {
  return *this = *this - right;
}

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
template<std::integral U>
Wide<Scalar>& Wide<Scalar>::operator*=(U const right)
  requires additive_group<Scalar> {
  return *this = *this * right;
}

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar>& Wide<Scalar>::operator*=(Wide<Scalar> const right)
  requires ring<Scalar> {
  return *this = *this * right;
}

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar>& Wide<Scalar>::operator/=(Wide<Scalar> const right)
  requires field<Scalar> {
  return *this = *this / right;
}

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar>::Wide(UnderlyingType const value) : value_(value) {}

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

template<typename T, std::integral U>
Wide<T> operator*(Wide<T> const left, U const right) {
  return left * Wide<double>(right);
}

template<typename T, std::integral U>
Wide<T> operator*(U const left, Wide<T> const right) {
  return Wide<double>(left) * right;
}

template<typename T, typename U>
Wide<Product<T, U>> FusedMultiplyAdd(Wide<T> const a,
                                     Wide<U> const b,
                                     Wide<Product<T, U>> const c) {
  return Wide<Product<T, U>>(MM(fmadd_pd)(a.value_, b.value_, c.value_));
}

template<typename T, typename U>
Wide<Product<T, U>> FusedMultiplySubtract(Wide<T> const a,
                                          Wide<U> const b,
                                          Wide<Product<T, U>> const c) {
  return Wide<Product<T, U>>(MM(fmsub_pd)(a.value_, b.value_, c.value_));
}

template<typename T, typename U>
Wide<Product<T, U>> FusedNegatedMultiplyAdd(Wide<T> const a,
                                            Wide<U> const b,
                                            Wide<Product<T, U>> const c) {
  return Wide<Product<T, U>>(MM(fnmadd_pd)(a.value_, b.value_, c.value_));
}

template<typename T, typename U>
Wide<Product<T, U>> FusedNegatedMultiplySubtract(Wide<T> const a,
                                                 Wide<U> const b,
                                                 Wide<Product<T, U>> const c) {
  return Wide<Product<T, U>>(MM(fnmsub_pd)(a.value_, b.value_, c.value_));
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
  return R3Element<Quotient<T, U>>(_mm256_div_pd(left.xyzt, right.value_));
#else
  return R3Element<Quotient<T, U>>(_mm_div_pd(left.xy, right.value_),
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

template<typename T, typename U>
R3Element<Product<T, U>> FusedMultiplyAdd(R3Element<T> const& a,
                                          Wide<U> const b,
                                          R3Element<Product<T, U>> const& c) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_fmadd_pd(a.xyzt, b.value_, c.xyzt));
#else
  return R3Element<Product<T, U>>(_mm_fmadd_pd(a.xy, b.value_, c.xy),
                                  _mm_fmadd_sd(a.zt, b.value_, c.zt));
#endif
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedMultiplyAdd(Wide<T> const a,
                                          R3Element<U> const& b,
                                          R3Element<Product<T, U>> const& c) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_fmadd_pd(a.value_, b.xyzt, c.xyzt));
#else
  return R3Element<Product<T, U>>(_mm_fmadd_pd(a.value_, b.xy, c.xy),
                                  _mm_fmadd_sd(a.value_, b.zt, c.zt));
#endif
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedMultiplySubtract(
    R3Element<T> const& a,
    Wide<U> b,
    R3Element<Product<T, U>> const& c) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_fmsub_pd(a.xyzt, b.value_, c.xyzt));
#else
  return R3Element<Product<T, U>>(_mm_fmsub_pd(a.xy, b.value_, c.xy),
                                  _mm_fmsub_sd(a.zt, b.value_, c.zt));
#endif
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedMultiplySubtract(
    Wide<T> a,
    R3Element<U> const& b,
    R3Element<Product<T, U>> const& c) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_fmsub_pd(a.value_, b.xyzt, c.xyzt));
#else
  return R3Element<Product<T, U>>(_mm_fmsub_pd(a.value_, b.xy, c.xy),
                                  _mm_fmsub_sd(a.value_, b.zt, c.zt));
#endif
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedNegatedMultiplyAdd(
    R3Element<T> const& a,
    Wide<U> b,
    R3Element<Product<T, U>> const& c) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_fnmadd_pd(a.xyzt, b.value_, c.xyzt));
#else
  return R3Element<Product<T, U>>(_mm_fnmadd_pd(a.xy, b.value_, c.xy),
                                  _mm_fnmadd_sd(a.zt, b.value_, c.zt));
#endif
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedNegatedMultiplyAdd(
    Wide<T> a,
    R3Element<U> const& b,
    R3Element<Product<T, U>> const& c) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_fnmadd_pd(a.value_, b.xyzt, c.xyzt));
#else
  return R3Element<Product<T, U>>(_mm_fnmadd_pd(a.value_, b.xy, c.xy),
                                  _mm_fnmadd_sd(a.value_, b.zt, c.zt));
#endif
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedNegatedMultiplySubtract(
    R3Element<T> const& a,
    Wide<U> b,
    R3Element<Product<T, U>> const& c) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_fnmsub_pd(a.xyzt, b.value_, c.xyzt));
#else
  return R3Element<Product<T, U>>(_mm_fnmsub_pd(a.xy, b.value_, c.xy),
                                  _mm_fnmsub_sd(a.zt, b.value_, c.zt));
#endif
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedNegatedMultiplySubtract(
    Wide<T> a,
    R3Element<U> const& b,
    R3Element<Product<T, U>> const& c) {
#if PRINCIPIA_USE_AVX()
  return R3Element<Product<T, U>>(_mm256_fnmsub_pd(a.value_, b.xyzt, c.xyzt));
#else
  return R3Element<Product<T, U>>(_mm_fnmsub_pd(a.value_, b.xy, c.xy),
                                  _mm_fnmsub_sd(a.value_, b.zt, c.zt));
#endif
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

template<r3_element_coordinates T, typename U>
Product<T, U> FusedMultiplyAdd(T const& a,
                               Wide<U> b,
                               Product<T, U> const& c) {
  return Product<T, U>(FusedMultiplyAdd(a.coordinates(), b, c.coordinates()));
}

template<typename T, r3_element_coordinates U>
Product<T, U> FusedMultiplyAdd(Wide<T> a,
                               U const& b,
                               Product<T, U> const& c) {
  return Product<T, U>(FusedMultiplyAdd(a, b.coordinates(), c.coordinates()));
}

template<r3_element_coordinates T, typename U>
Product<T, U> FusedMultiplySubtract(T const& a,
                                    Wide<U> b,
                                    Product<T, U> const& c) {
  return Product<T, U>(
      FusedMultiplySubtract(a.coordinates(), b, c.coordinates()));
}

template<typename T, r3_element_coordinates U>
Product<T, U> FusedMultiplySubtract(Wide<T> a,
                                    U const& b,
                                    Product<T, U> const& c) {
  return Product<T, U>(
      FusedMultiplySubtract(a, b.coordinates(), c.coordinates()));
}

template<r3_element_coordinates T, typename U>
Product<T, U> FusedNegatedMultiplyAdd(T const& a,
                                      Wide<U> b,
                                      Product<T, U> const& c) {
  return Product<T, U>(
      FusedNegatedMultiplyAdd(a.coordinates(), b, c.coordinates()));
}

template<typename T, r3_element_coordinates U>
Product<T, U> FusedNegatedMultiplyAdd(Wide<T> a,
                                      U const& b,
                                      Product<T, U> const& c) {
  return Product<T, U>(
      FusedNegatedMultiplyAdd(a, b.coordinates(), c.coordinates()));
}

template<r3_element_coordinates T, typename U>
Product<T, U> FusedNegatedMultiplySubtract(T const& a,
                                           Wide<U> b,
                                           Product<T, U> const& c) {
  return Product<T, U>(
      FusedNegatedMultiplySubtract(a.coordinates(), b, c.coordinates()));
}

template<typename T, r3_element_coordinates U>
Product<T, U> FusedNegatedMultiplySubtract(Wide<T> a,
                                           U const& b,
                                           Product<T, U> const& c) {
  return Product<T, U>(
      FusedNegatedMultiplySubtract(a, b.coordinates(), c.coordinates()));
}

#undef MM

}  // namespace internal
}  // namespace _wide
}  // namespace numerics
}  // namespace principia
