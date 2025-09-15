#pragma once

#include "geometry/r3_element.hpp"

#include <pmmintrin.h>

#include <string>
#include <type_traits>

#include "base/macros.hpp"  // 🧙 For PRINCIPIA_USE_SSE3_INTRINSICS.
#include "glog/logging.h"
#include "numerics/elementary_functions.hpp"
#include "numerics/fma.hpp"
#include "quantities/serialization.hpp"

namespace principia {
namespace geometry {
namespace _r3_element {
namespace internal {

using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_fma;
using namespace principia::quantities::_serialization;

// We want zero initialization here, so the default constructor won't do.
template<typename Scalar>
constexpr R3Element<Scalar>::R3Element() : x(), y(), z() {
  static_assert(std::is_standard_layout<R3Element>::value,
                "R3Element has a nonstandard layout");
}

template<typename Scalar>
constexpr R3Element<Scalar>::R3Element(uninitialized_t) {
  static_assert(std::is_standard_layout<R3Element>::value,
                "R3Element has a nonstandard layout");
}

template<typename Scalar>
R3Element<Scalar>::R3Element(Scalar const& x,
                             Scalar const& y,
                             Scalar const& z) : x(x), y(y), z(z) {
  static_assert(std::is_standard_layout<R3Element>::value,
                "R3Element has a nonstandard layout");
}

template<typename Scalar>
R3Element<Scalar>::R3Element(__m128d const xy, __m128d const zt)
    : xy(xy), zt(zt) {
  static_assert(std::is_standard_layout<R3Element>::value,
                "R3Element has a nonstandard layout");
}

template<typename Scalar>
Scalar& R3Element<Scalar>::operator[](int const index) {
  switch (index) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      DLOG(FATAL) << FUNCTION_SIGNATURE << ": index = " << index;
      std::abort();
  }
}

template<typename Scalar>
Scalar const& R3Element<Scalar>::operator[](
    int const index) const {
  switch (index) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      DLOG(FATAL) << FUNCTION_SIGNATURE << ": index = " << index;
      std::abort();
  }
}

template<typename Scalar>
R3Element<Scalar>& R3Element<Scalar>::operator+=(
    R3Element<Scalar> const& right) {
#if PRINCIPIA_USE_SSE3_INTRINSICS()
  xy = _mm_add_pd(xy, right.xy);
  zt = _mm_add_sd(zt, right.zt);
#else
  x += right.x;
  y += right.y;
  z += right.z;
#endif
  return *this;
}

template<typename Scalar>
R3Element<Scalar>& R3Element<Scalar>::operator-=(
    R3Element<Scalar> const& right) {
#if PRINCIPIA_USE_SSE3_INTRINSICS()
  xy = _mm_sub_pd(xy, right.xy);
  zt = _mm_sub_sd(zt, right.zt);
#else
  x -= right.x;
  y -= right.y;
  z -= right.z;
#endif
  return *this;
}

template<typename Scalar>
R3Element<Scalar>& R3Element<Scalar>::operator*=(double const right) {
#if PRINCIPIA_USE_SSE3_INTRINSICS()
  __m128d const right_128d = ToM128D(right);
  xy = _mm_mul_pd(xy, right_128d);
  zt = _mm_mul_sd(zt, right_128d);
#else
  x *= right;
  y *= right;
  z *= right;
#endif
  return *this;
}

template<typename Scalar>
R3Element<Scalar>& R3Element<Scalar>::operator/=(double const right) {
#if PRINCIPIA_USE_SSE3_INTRINSICS()
  __m128d const right_128d = ToM128D(right);
  xy = _mm_div_pd(xy, right_128d);
  zt = _mm_div_sd(zt, right_128d);
#else
  x /= right;
  y /= right;
  z /= right;
#endif
  return *this;
}

template<typename Scalar>
Scalar R3Element<Scalar>::Norm() const {
  return Sqrt(Norm²());
}

template<typename Scalar>
Square<Scalar> R3Element<Scalar>::Norm²() const {
  return x * x + y * y + z * z;
}

template<typename Scalar>
SphericalCoordinates<Scalar> R3Element<Scalar>::ToSpherical() const {
  SphericalCoordinates<Scalar> result;
  result.radius = Norm();
  result.latitude = ArcSin(z / result.radius);
  result.longitude = ArcTan(y, x);
  return result;
}

template<typename Scalar>
template<typename S>
R3Element<Scalar> R3Element<Scalar>::OrthogonalizationAgainst(
    R3Element<S> const& r3_element) const {
  R3Element<double> const r3_element_normalized = Normalize(r3_element);
  return *this - Dot(*this, r3_element_normalized) * r3_element_normalized;
}

template<typename Scalar>
void R3Element<Scalar>::WriteToMessage(
    not_null<serialization::R3Element*> const message) const {
  using Serializer =
      DoubleOrQuantitySerializer<Scalar, serialization::R3Element::Coordinate>;
  Serializer::WriteToMessage(x, message->mutable_x());
  Serializer::WriteToMessage(y, message->mutable_y());
  Serializer::WriteToMessage(z, message->mutable_z());
}

template<typename Scalar>
R3Element<Scalar> R3Element<Scalar>::ReadFromMessage(
    serialization::R3Element const& message) {
  using Serializer =
      DoubleOrQuantitySerializer<Scalar, serialization::R3Element::Coordinate>;
  return {Serializer::ReadFromMessage(message.x()),
          Serializer::ReadFromMessage(message.y()),
          Serializer::ReadFromMessage(message.z())};
}

template<typename Scalar>
SphericalCoordinates<Scalar>::SphericalCoordinates() {}

template<typename Scalar>
R3Element<Scalar> SphericalCoordinates<Scalar>::ToCartesian() {
  auto const [sin_latitude, cos_latitude] = SinCos(latitude);
  auto const [sin_longitude, cos_longitude] = SinCos(longitude);
  return {radius * cos_longitude * cos_latitude,
          radius * sin_longitude * cos_latitude,
          radius * sin_latitude};
}

template<typename Scalar>
SphericalCoordinates<Scalar> RadiusLatitudeLongitude(Scalar const& radius,
                                                     Angle const& latitude,
                                                     Angle const& longitude) {
  SphericalCoordinates<Scalar> result;
  result.latitude = latitude;
  result.longitude = longitude;
  result.radius = radius;
  return result;
}

template<typename Scalar>
std::ostream& operator<<(
    std::ostream& out,
    SphericalCoordinates<Scalar> const& spherical_coordinates) {
  out << spherical_coordinates.radius
      << ", " << spherical_coordinates.latitude
      << ", " << spherical_coordinates.longitude;
  return out;
}

template<typename Scalar>
R3Element<Scalar> operator+(R3Element<Scalar> const& right) {
  return R3Element<Scalar>(+right.x, +right.y, +right.z);
}

template<typename Scalar>
R3Element<Scalar> operator-(R3Element<Scalar> const& right) {
  return R3Element<Scalar>(-right.x, -right.y, -right.z);
}

template<typename Scalar>
R3Element<Scalar> operator+(R3Element<Scalar> const& left,
                            R3Element<Scalar> const& right) {
#if PRINCIPIA_USE_SSE3_INTRINSICS()
  return R3Element<Scalar>(_mm_add_pd(left.xy, right.xy),
                           _mm_add_sd(left.zt, right.zt));
#else
  return R3Element<Scalar>(left.x + right.x,
                           left.y + right.y,
                           left.z + right.z);
#endif
}

template<typename Scalar>
R3Element<Scalar> operator-(R3Element<Scalar> const& left,
                            R3Element<Scalar> const& right) {
#if PRINCIPIA_USE_SSE3_INTRINSICS()
  return R3Element<Scalar>(_mm_sub_pd(left.xy, right.xy),
                           _mm_sub_sd(left.zt, right.zt));
#else
  return R3Element<Scalar>(left.x - right.x,
                           left.y - right.y,
                           left.z - right.z);
#endif
}

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<LScalar>
R3Element<Product<LScalar, RScalar>> operator*(
    LScalar const& left,
    R3Element<RScalar> const& right) {
#if PRINCIPIA_USE_SSE3_INTRINSICS()
  __m128d const left_128d = ToM128D(left);
  return R3Element<Product<LScalar, RScalar>>(_mm_mul_pd(right.xy, left_128d),
                                              _mm_mul_sd(right.zt, left_128d));
#else
  return R3Element<Product<LScalar, RScalar>>(left * right.x,
                                              left * right.y,
                                              left * right.z);
#endif
}

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3Element<Product<LScalar, RScalar>> operator*(R3Element<LScalar> const& left,
                                               RScalar const& right) {
#if PRINCIPIA_USE_SSE3_INTRINSICS()
  __m128d const right_128d = ToM128D(right);
  return R3Element<Product<LScalar, RScalar>>(_mm_mul_pd(left.xy, right_128d),
                                              _mm_mul_sd(left.zt, right_128d));
#else
  return R3Element<Product<LScalar, RScalar>>(left.x * right,
                                              left.y * right,
                                              left.z * right);
#endif
}

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3Element<Quotient<LScalar, RScalar>> operator/(R3Element<LScalar> const& left,
                                                RScalar const& right) {
#if PRINCIPIA_USE_SSE3_INTRINSICS()
  __m128d const right_128d = ToM128D(right);
  return R3Element<Quotient<LScalar, RScalar>>(_mm_div_pd(left.xy, right_128d),
                                               _mm_div_sd(left.zt, right_128d));
#else
  return R3Element<Quotient<LScalar, RScalar>>(left.x / right,
                                               left.y / right,
                                               left.z / right);
#endif
}

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3Element<Product<LScalar, RScalar>> FusedMultiplyAdd(
    R3Element<LScalar> const& a,
    RScalar const& b,
    R3Element<Product<LScalar, RScalar>> const& c) {
  if constexpr (CanEmitFMAInstructions) {
    __m128d const b_128d = ToM128D(b);
    return R3Element<Product<LScalar, RScalar>>(
        _mm_fmadd_pd(a.xy, b_128d, c.xy), _mm_fmadd_sd(a.zt, b_128d, c.zt));
  } else {
    LOG(FATAL) << "Clang cannot use FMA without VEX-encoding everything";
  }
}

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3Element<Product<LScalar, RScalar>> FusedMultiplySubtract(
    R3Element<LScalar> const& a,
    RScalar const& b,
    R3Element<Product<LScalar, RScalar>> const& c) {
  if constexpr (CanEmitFMAInstructions) {
    __m128d const b_128d = ToM128D(b);
    return R3Element<Product<LScalar, RScalar>>(
        _mm_fmsub_pd(a.xy, b_128d, c.xy), _mm_fmsub_sd(a.zt, b_128d, c.zt));
  } else {
    LOG(FATAL) << "Clang cannot use FMA without VEX-encoding everything";
  }
}

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3Element<Product<LScalar, RScalar>> FusedNegatedMultiplyAdd(
    R3Element<LScalar> const& a,
    RScalar const& b,
    R3Element<Product<LScalar, RScalar>> const& c) {
  if constexpr (CanEmitFMAInstructions) {
    __m128d const b_128d = ToM128D(b);
    return R3Element<Product<LScalar, RScalar>>(
        _mm_fnmadd_pd(a.xy, b_128d, c.xy), _mm_fnmadd_sd(a.zt, b_128d, c.zt));
  } else {
    LOG(FATAL) << "Clang cannot use FMA without VEX-encoding everything";
  }
}

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3Element<Product<LScalar, RScalar>> FusedNegatedMultiplySubtract(
    R3Element<LScalar> const& a,
    RScalar const& b,
    R3Element<Product<LScalar, RScalar>> const& c) {
  if constexpr (CanEmitFMAInstructions) {
    __m128d const b_128d = ToM128D(b);
    return R3Element<Product<LScalar, RScalar>>(
        _mm_fnmsub_pd(a.xy, b_128d, c.xy), _mm_fnmsub_sd(a.zt, b_128d, c.zt));
  } else {
    LOG(FATAL) << "Clang cannot use FMA without VEX-encoding everything";
  }
}

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<LScalar>
R3Element<Product<LScalar, RScalar>> FusedMultiplyAdd(
    LScalar const& a,
    R3Element<RScalar> const& b,
    R3Element<Product<LScalar, RScalar>> const& c) {
  if constexpr (CanEmitFMAInstructions) {
    __m128d const a_128d = ToM128D(a);
    return R3Element<Product<LScalar, RScalar>>(
        _mm_fmadd_pd(a_128d, b.xy, c.xy), _mm_fmadd_sd(a_128d, b.zt, c.zt));
  } else {
    LOG(FATAL) << "Clang cannot use FMA without VEX-encoding everything";
  }
}

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<LScalar>
R3Element<Product<LScalar, RScalar>> FusedMultiplySubtract(
    LScalar const& a,
    R3Element<RScalar> const& b,
    R3Element<Product<LScalar, RScalar>> const& c) {
  if constexpr (CanEmitFMAInstructions) {
    __m128d const a_128d = ToM128D(a);
    return R3Element<Product<LScalar, RScalar>>(
        _mm_fmsub_pd(a_128d, b.xy, c.xy), _mm_fmsub_sd(a_128d, b.zt, c.zt));
  } else {
    LOG(FATAL) << "Clang cannot use FMA without VEX-encoding everything";
  }
}

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<LScalar>
R3Element<Product<LScalar, RScalar>> FusedNegatedMultiplyAdd(
    LScalar const& a,
    R3Element<RScalar> const& b,
    R3Element<Product<LScalar, RScalar>> const& c) {
  if constexpr (CanEmitFMAInstructions) {
    __m128d const a_128d = ToM128D(a);
    return R3Element<Product<LScalar, RScalar>>(
        _mm_fnmadd_pd(a_128d, b.xy, c.xy), _mm_fnmadd_sd(a_128d, b.zt, c.zt));
  } else {
    LOG(FATAL) << "Clang cannot use FMA without VEX-encoding everything";
  }
}

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<LScalar>
R3Element<Product<LScalar, RScalar>> FusedNegatedMultiplySubtract(
    LScalar const& a,
    R3Element<RScalar> const& b,
    R3Element<Product<LScalar, RScalar>> const& c) {
  if constexpr (CanEmitFMAInstructions) {
    __m128d const a_128d = ToM128D(a);
    return R3Element<Product<LScalar, RScalar>>(
        _mm_fnmsub_pd(a_128d, b.xy, c.xy), _mm_fnmsub_sd(a_128d, b.zt, c.zt));
  } else {
    LOG(FATAL) << "Clang cannot use FMA without VEX-encoding everything";
  }
}

template<typename Scalar>
constexpr bool operator==(R3Element<Scalar> const& left,
                          R3Element<Scalar> const& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

template<typename Scalar>
constexpr bool operator!=(R3Element<Scalar> const& left,
                          R3Element<Scalar> const& right) {
  return left.x != right.x || left.y != right.y || left.z != right.z;
}

template<typename Scalar>
R3Element<double> Normalize(R3Element<Scalar> const& r3_element) {
  Scalar const norm = r3_element.Norm();
#ifdef _DEBUG
  CHECK_NE(Scalar(), norm);
#endif
  return r3_element / norm;
}

template<typename Scalar>
R3Element<double> NormalizeOrZero(R3Element<Scalar> const& r3_element) {
  Scalar const norm = r3_element.Norm();
  if (norm == Scalar()) {
    static R3Element<double> const zeroes = {0, 0, 0};
    return zeroes;
  } else {
    return r3_element / norm;
  }
}

template<typename Scalar>
std::string DebugString(R3Element<Scalar> const& r3_element) {
  using quantities::_quantities::DebugString;
  std::string result = "{";
  result += DebugString(r3_element.x);
  result += ", ";
  result += DebugString(r3_element.y);
  result += ", ";
  result += DebugString(r3_element.z);
  result +="}";
  return result;
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         R3Element<Scalar> const& r3_element) {
  out << DebugString(r3_element);
  return out;
}

template<typename LScalar, typename RScalar>
R3Element<Product<LScalar, RScalar>> Cross(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right) {
  return R3Element<Product<LScalar, RScalar>>(
      left.y * right.z - left.z * right.y,
      left.z * right.x - left.x * right.z,
      left.x * right.y - left.y * right.x);
}

template<typename LScalar, typename RScalar>
Product<LScalar, RScalar> Dot(R3Element<LScalar> const& left,
                              R3Element<RScalar> const& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

// Implementation from [Kah06], §12 "Mangled Angles", p. 47.
template<typename LScalar, typename RScalar>
Angle AngleBetween(R3Element<LScalar> const& left,
                   R3Element<RScalar> const& right) {
  auto const left_norm_right = left * right.Norm();
  auto const right_norm_left = right * left.Norm();
  return 2 * ArcTan((left_norm_right - right_norm_left).Norm(),
                    (left_norm_right + right_norm_left).Norm());
}

inline R3Element<double> BasisVector(int const i) {
  DCHECK_GE(i, 0) << i;
  DCHECK_LT(i, 3) << i;
  return {static_cast<double>(i == 0),
          static_cast<double>(i == 1),
          static_cast<double>(i == 2)};
}

}  // namespace internal
}  // namespace _r3_element
}  // namespace geometry
}  // namespace principia
