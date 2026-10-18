#pragma once

#include "numerics/wide.hpp"

#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _wide {
namespace internal {

using namespace principia::quantities::_quantities;

// A wrapper for result types that are template instances.
#define TYPE(...) __VA_ARGS__

// In the macros below, 'R' stands for 'R3Element' and 'W' for 'Wide'.
#if PRINCIPIA_USE_AVX()
#define MM_PD(intrinsic) _mm256_##intrinsic##_pd
#define INTRINSIC_RW(type, intrinsic, left, right) \
  R3Element<type>(MM_PD(intrinsic)(left.xyzt, right.value_))
#define INTRINSIC_RWR(type, intrinsic, a, b, c) \
  R3Element<type>(MM_PD(intrinsic)(a.xyzt, b.value_, c.xyzt))
#define INTRINSIC_WR(type, intrinsic, left, right) \
  R3Element<type>(MM_PD(intrinsic)(left.value_, right.xyzt))
#define INTRINSIC_WW(type, intrinsic, left, right) \
  Wide<type>(MM_PD(intrinsic)(left.value_, right.value_))
#define INTRINSIC_WRR(type, intrinsic, a, b, c) \
  R3Element<type>(MM_PD(intrinsic)(a.value_, b.xyzt, c.xyzt))
#define INTRINSIC_WWW(type, intrinsic, a, b, c) \
  Wide<type>(MM_PD(intrinsic)(a.value_, b.value_, c.value_));
#else
#define MM_PD(intrinsic) _mm_##intrinsic##_pd
#define MM_SD(intrinsic) _mm_##intrinsic##_sd
#define INTRINSIC_RW(type, intrinsic, left, right)         \
  R3Element<type>(MM_PD(intrinsic)(left.xy, right.value_), \
                  MM_SD(intrinsic)(left.zt, right.value_))
#define INTRINSIC_RWR(type, intrinsic, a, b, c)           \
  R3Element<type>(MM_PD(intrinsic)(a.xy, b.value_, c.xy), \
                  MM_SD(intrinsic)(a.zt, b.value_, c.zt))
#define INTRINSIC_WR(type, intrinsic, left, right)         \
  R3Element<type>(MM_PD(intrinsic)(left.value_, right.xy), \
                  MM_SD(intrinsic)(left.value_, right.zt))
#define INTRINSIC_WW(type, intrinsic, left, right) \
  Wide<type>(MM_PD(intrinsic)(left.value_, right.value_))
#define INTRINSIC_WRR(type, intrinsic, a, b, c)           \
  R3Element<type>(MM_PD(intrinsic)(a.value_, b.xy, c.xy), \
                  MM_SD(intrinsic)(a.value_, b.zt, c.zt))
#define INTRINSIC_WWW(type, intrinsic, a, b, c) \
  Wide<type>(MM_PD(intrinsic)(a.value_, b.value_, c.value_));
#endif

inline UnderlyingType wide_sign_bit = MM_PD(set1)(-0.0);

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
Wide<Scalar>::operator Scalar() const {
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
Wide<Scalar>& Wide<Scalar>::operator+=(Wide<Difference<Scalar>> const right) {
  return *this = *this + right;
}

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
Wide<Scalar>& Wide<Scalar>::operator-=(Wide<Difference<Scalar>> const right) {
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
template<ring R>
Wide<Scalar>& Wide<Scalar>::operator*=(Wide<R> const right)
  requires module<Scalar, R> {
  return *this = *this * right;
}

template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
template<field K>
Wide<Scalar>& Wide<Scalar>::operator/=(Wide<K> const right)
  requires vector_space<Scalar, K> {
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
  return Wide<T>(MM_PD(xor)(right.value_, wide_sign_bit));
}

template<affine T, affine U>
Wide<Sum<T, U>> operator+(Wide<T> const left, Wide<U> const right) {
  return INTRINSIC_WW(TYPE(Sum<T, U>), add, left, right);
}

template<affine T, affine U>
Wide<Difference<T, U>> operator-(Wide<T> const left, Wide<U> const right) {
  return INTRINSIC_WW(TYPE(Difference<T, U>), sub, left, right);
}

template<typename T, typename U>
Wide<Product<T, U>> operator*(Wide<T> const left, Wide<U> const right) {
  return INTRINSIC_WW(TYPE(Product<T, U>), mul, left, right);
}

template<typename T, typename U>
Wide<Quotient<T, U>> operator/(Wide<T> const left, Wide<U> const right) {
  return INTRINSIC_WW(TYPE(Quotient<T, U>), div, left, right);
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
  return INTRINSIC_WWW(TYPE(Product<T, U>), fmadd, a, b, c);
}

template<typename T, typename U>
Wide<Product<T, U>> FusedMultiplySubtract(Wide<T> const a,
                                          Wide<U> const b,
                                          Wide<Product<T, U>> const c) {
  return INTRINSIC_WWW(TYPE(Product<T, U>), fmsub, a, b, c);
}

template<typename T, typename U>
Wide<Product<T, U>> FusedNegatedMultiplyAdd(Wide<T> const a,
                                            Wide<U> const b,
                                            Wide<Product<T, U>> const c) {
  return INTRINSIC_WWW(TYPE(Product<T, U>), fnmadd, a, b, c);
}

template<typename T, typename U>
Wide<Product<T, U>> FusedNegatedMultiplySubtract(Wide<T> const a,
                                                 Wide<U> const b,
                                                 Wide<Product<T, U>> const c) {
  return INTRINSIC_WWW(TYPE(Product<T, U>), fnmsub, a, b, c);
}

template<typename T, typename U>
R3Element<Product<T, U>> operator*(R3Element<T> const& left,
                                   Wide<U> const right) {
  return INTRINSIC_RW(TYPE(Product<T, U>), mul, left, right);
}

template<typename T, typename U>
R3Element<Product<T, U>> operator*(Wide<T> const left,
                                   R3Element<U> const& right) {
  return INTRINSIC_WR(TYPE(Product<T, U>), mul, left, right);
}

template<typename T, typename U>
R3Element<Quotient<T, U>> operator/(R3Element<T> const& left,
                                    Wide<U> const right) {
  return INTRINSIC_RW(TYPE(Quotient<T, U>), div, left, right);
}

template<typename T>
R3Element<T>& operator*=(R3Element<T>& left, Wide<double> const right) {
  return left = INTRINSIC_RW(TYPE(Product<T, double>), mul, left, right);
}

template<typename T>
R3Element<T>& operator/=(R3Element<T>& left, Wide<double> const right) {
  return left = INTRINSIC_RW(TYPE(Quotient<T, double>), div, left, right);
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedMultiplyAdd(R3Element<T> const& a,
                                          Wide<U> const b,
                                          R3Element<Product<T, U>> const& c) {
  return INTRINSIC_RWR(TYPE(Product<T, U>), fmadd, a, b, c);
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedMultiplyAdd(Wide<T> const a,
                                          R3Element<U> const& b,
                                          R3Element<Product<T, U>> const& c) {
  return INTRINSIC_WRR(TYPE(Product<T, U>), fmadd, a, b, c);
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedMultiplySubtract(
    R3Element<T> const& a,
    Wide<U> const b,
    R3Element<Product<T, U>> const& c) {
  return INTRINSIC_RWR(TYPE(Product<T, U>), fmsub, a, b, c);
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedMultiplySubtract(
    Wide<T> const a,
    R3Element<U> const& b,
    R3Element<Product<T, U>> const& c) {
  return INTRINSIC_WRR(TYPE(Product<T, U>), fmsub, a, b, c);
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedNegatedMultiplyAdd(
    R3Element<T> const& a,
    Wide<U> const b,
    R3Element<Product<T, U>> const& c) {
  return INTRINSIC_RWR(TYPE(Product<T, U>), fnmadd, a, b, c);
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedNegatedMultiplyAdd(
    Wide<T> const a,
    R3Element<U> const& b,
    R3Element<Product<T, U>> const& c) {
  return INTRINSIC_WRR(TYPE(Product<T, U>), fnmadd, a, b, c);
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedNegatedMultiplySubtract(
    R3Element<T> const& a,
    Wide<U> const b,
    R3Element<Product<T, U>> const& c) {
  return INTRINSIC_RWR(TYPE(Product<T, U>), fnmsub, a, b, c);
}

template<typename T, typename U>
R3Element<Product<T, U>> FusedNegatedMultiplySubtract(
    Wide<T> const a,
    R3Element<U> const& b,
    R3Element<Product<T, U>> const& c) {
  return INTRINSIC_WRR(TYPE(Product<T, U>), fnmsub, a, b, c);
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

#undef TYPE
#undef MM_PD
#undef MM_SD
#undef INTRINSIC_RW
#undef INTRINSIC_RWR
#undef INTRINSIC_WR
#undef INTRINSIC_WW
#undef INTRINSIC_WRR
#undef INTRINSIC_WWW

}  // namespace internal
}  // namespace _wide
}  // namespace numerics
}  // namespace principia
