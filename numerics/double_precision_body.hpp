#pragma once

#include "numerics/double_precision.hpp"

#include <array>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>

#include "base/not_constructible.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/interval.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/serialization.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/m128d.hpp"
#include "quantities/concepts.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _double_precision {
namespace internal {

using namespace principia::base::_not_constructible;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_interval;
using namespace principia::geometry::_point;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_serialization;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_m128d;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_si;

// This returns the possible exponent interval for representing `x` as a
// bounded floating point number in the sense of [DRT01], section 2.2 and
// hypothesis `pGivesBound`.  Essentially this corresponts to a bounded
// (integral) mantissa not required to be normalized (thus, larger exponents are
// acceptable if the mantissa has trailing zeroes).
inline Interval<int> ExponentInterval(double const x) {
  static constexpr int M = std::numeric_limits<double>::digits;
  int exponent;
  double const double_mantissa = std::abs(std::frexp(x, &exponent));
  DCHECK_LE(0.5, double_mantissa);
  DCHECK_LT(double_mantissa, 1.0);
  std::uint64_t const integer_mantissa =
      static_cast<std::uint64_t>(double_mantissa * (1LL << M));
  DCHECK_LE(1LL << (M - 1), integer_mantissa);
  DCHECK_LT(integer_mantissa, 1LL << M);
  return {.min = exponent - M,
          .max = exponent - M + std::countr_zero(integer_mantissa)};
}

// A helper to check that the preconditions of QuickTwoSum are met.  Annoyingly
// complicated as it needs to peel off all of our abstractions until it reaches
// doubles.
template<typename T, typename U, typename = void>
struct DaumasRideauThéryComparator;

template<typename T, typename U>
struct DaumasRideauThéryComparator<DoublePrecision<T>, DoublePrecision<U>> {
  static bool ProperlyOrdered(DoublePrecision<T> const& biggish,
                              DoublePrecision<U> const& smallish) {
    // This check is incomplete: it doesn't compare the errors component.  To
    // the best of my understanding this code is only used in tests.
    return DaumasRideauThéryComparator<T, U>::ProperlyOrdered(biggish.value,
                                                              smallish.value);
  }
};

template<typename T, typename U>
struct DaumasRideauThéryComparator<Point<T>, U> : not_constructible {
  static bool ProperlyOrdered(Point<T> const& biggish, U const& smallish) {
    // We only care about the coordinates, the geometric structure is
    // irrelevant.
    return DaumasRideauThéryComparator<T, U>::ProperlyOrdered(
        biggish - Point<T>{}, smallish);
  }
};

template<typename T, typename U>
struct DaumasRideauThéryComparator<T, Point<U>> : not_constructible {
  static bool ProperlyOrdered(T const& biggish, Point<U> const& smallish) {
    // We only care about the coordinates, the geometric structure is
    // irrelevant.
    return DaumasRideauThéryComparator<T, U>::ProperlyOrdered(
        biggish, smallish - Point<U>{});
  }
};

template<typename T, typename TFrame, int trank,
         typename U, typename UFrame, int urank>
struct DaumasRideauThéryComparator<Multivector<T, TFrame, trank>,
                                   Multivector<U, UFrame, urank>>
    : not_constructible {
  static bool ProperlyOrdered(Multivector<T, TFrame, trank> const& biggish,
                              Multivector<U, UFrame, urank> const& smallish) {
    // This doesn't handle trivectors.
    return DaumasRideauThéryComparator<R3Element<T>, R3Element<U>>::
        ProperlyOrdered(biggish.coordinates(), smallish.coordinates());
  }
};

template<typename T, typename U>
struct DaumasRideauThéryComparator<R3Element<T>, R3Element<U>>
    : not_constructible {
  static bool ProperlyOrdered(R3Element<T> const& biggish,
                              R3Element<U> const& smallish) {
    bool result = true;
    for (int i = 0; i < 3; ++i) {
      result &= DaumasRideauThéryComparator<T, U>::ProperlyOrdered(biggish[i],
                                                                   smallish[i]);
    }
    return result;
  }
};

template<>
struct DaumasRideauThéryComparator<double, double> {
  static bool ProperlyOrdered(double const biggish, double const smallish) {
    // In the elementary functions, we use NaN to fall back to the CORE-MATH
    // implementation.  We don't want to die because of the weird comparisons of
    // NaNs.
    if (biggish == 0.0 || smallish == 0.0 ||
        biggish != biggish || smallish != smallish) {
      return true;
    } else {
      // We don't use the classical check `|biggish| >= |smallish|` because the
      // Boldo- Daumas-Li reduction needs the more precise check from [DRT01],
      // section 6.1, p. 179.
      auto const biggish_exponent_interval = ExponentInterval(biggish);
      auto const smallish_exponent_interval = ExponentInterval(smallish);
      return smallish_exponent_interval.min <= biggish_exponent_interval.max;
    }
  }
};

template<>
struct DaumasRideauThéryComparator<M128D, M128D> {
  static bool ProperlyOrdered(M128D const& biggish, M128D const& smallish) {
    return DaumasRideauThéryComparator<double, double>::ProperlyOrdered(
        static_cast<double>(biggish), static_cast<double>(smallish));
  }
};

template<typename T, typename U>
  requires convertible_to_quantity<T> && convertible_to_quantity<U>
struct DaumasRideauThéryComparator<T, U> {
  static bool ProperlyOrdered(T const& biggish, U const& smallish) {
    return DaumasRideauThéryComparator<double, double>::ProperlyOrdered(
        biggish / si::Unit<T>, smallish / si::Unit<U>);
  }
};


template<typename T>
constexpr DoublePrecision<T>::DoublePrecision() : value{}, error{} {}

template<typename T>
constexpr DoublePrecision<T>::DoublePrecision(uninitialized_t) {}

template<typename T>
constexpr DoublePrecision<T>::DoublePrecision(T const& value)
    : value(value),
      error() {}

template<typename T>
DoublePrecision<T>& DoublePrecision<T>::operator+=(
    DoublePrecision<Difference<T>> const& right) {
  *this = *this + right;
  return *this;
}

template<typename T>
DoublePrecision<T>& DoublePrecision<T>::operator+=(
    Difference<T> const& right) {
  *this = *this + right;
  return *this;
}

template<typename T>
DoublePrecision<T>& DoublePrecision<T>::operator-=(
    DoublePrecision<Difference<T>> const& right) {
  *this = *this - right;
  return *this;
}

template<typename T>
DoublePrecision<T>& DoublePrecision<T>::operator-=(
    Difference<T> const& right) {
  *this = *this - right;
  return *this;
}

template<typename T>
template<typename U>
DoublePrecision<T>& DoublePrecision<T>::operator*=(U const& right) {
  *this = *this * right;
  return *this;
}

template<typename T>
template<typename U>
DoublePrecision<T>& DoublePrecision<T>::operator/=(U const& right) {
  *this = *this / right;
  return *this;
}

template <typename T>
DoublePrecision<T>& DoublePrecision<T>::Decrement(Difference<T> const& right) {
  // See Higham, Accuracy and Stability of Numerical Algorithms, Algorithm 4.2.
  // This is equivalent to `QuickTwoSum(value, error - right)`.
  T const temp = value;
  Difference<T> const y = error - right;
  value = temp + y;
  error = (temp - value) + y;
  return *this;
}

template <typename T>
DoublePrecision<T>& DoublePrecision<T>::Increment(Difference<T> const& right) {
  // See Higham, Accuracy and Stability of Numerical Algorithms, Algorithm 4.2.
  // This is equivalent to `QuickTwoSum(value, error + right)`.
  T const temp = value;
  Difference<T> const y = error + right;
  value = temp + y;
  error = (temp - value) + y;
  return *this;
}

template<typename T>
void DoublePrecision<T>::WriteToMessage(
    not_null<serialization::DoublePrecision*> const message) const {
  using ValueSerializer = DoubleOrQuantityOrPointOrMultivectorSerializer<
                              T, serialization::DoublePrecision::Value>;
  using ErrorSerializer = DoubleOrQuantityOrMultivectorSerializer<
                              Difference<T>,
                              serialization::DoublePrecision::Error>;
  ValueSerializer::WriteToMessage(value, message->mutable_value());
  ErrorSerializer::WriteToMessage(error, message->mutable_error());
}

template<typename T>
DoublePrecision<T> DoublePrecision<T>::ReadFromMessage(
    serialization::DoublePrecision const& message) {
  using ValueSerializer = DoubleOrQuantityOrPointOrMultivectorSerializer<
                              T, serialization::DoublePrecision::Value>;
  using ErrorSerializer = DoubleOrQuantityOrMultivectorSerializer<
                              Difference<T>,
                              serialization::DoublePrecision::Error>;
  DoublePrecision double_precision;
  double_precision.value = ValueSerializer::ReadFromMessage(message.value());
  double_precision.error = ErrorSerializer::ReadFromMessage(message.error());
  return double_precision;
}

template<typename T, typename U>
DoublePrecision<Product<T, U>> Scale(T const & scale,
                                     DoublePrecision<U> const& right) {
#ifdef _DEBUG
  double const s = scale / si::Unit<T>;
  if (s != 0.0) {
    int exponent;
    double const mantissa = std::frexp(s, &exponent);
    CHECK_EQ(0.5, std::fabs(mantissa)) << scale;
  }
#endif
  DoublePrecision<Product<T, U>> result(uninitialized);
  result.value = right.value * scale;
  result.error = right.error * scale;
  return result;
}

template<std::floating_point T>
DoublePrecision<T> Trunc(DoublePrecision<T> const& a) {
  // The fractional point cannot be in both `value` and `error`.  It is
  // important to round both elements consistently based on the sign of `value`,
  // so we cannot use `trunc`.
  DoublePrecision<T> result(uninitialized);
  result.value = a.value > 0.0 ? std::floor(a.value) : std::ceil(a.value);
  result.error = a.value > 0.0 ? std::floor(a.error) : std::ceil(a.error);
  return result;
}

template<std::floating_point T>
DoublePrecision<T> Frac(DoublePrecision<T> const& a) {
  return a - Trunc(a);
}

template<std::int64_t s, typename T>
constexpr void VeltkampSplitting(T const& x, T& xh, T& xl) {
  // We used to have a bug where this constant was just `1 << s`.
  // [Dek71, p. 234] gives counter examples where this does not work when the
  // number of mantissa bits is even, or when C < 0, but it's unclear in what
  // cases our mistake would have caused incorrect results.
  constexpr double C = (1 << s) + 1;
  T γ;
  if constexpr (real_vector_space<T>) {
    γ = C * x;
  } else if constexpr (ring<T>) {
    T const Ct(C);
    γ = Ct * x;
  } else {
    static_assert(false);
  }
  T const δ = x - γ;
  xh = γ + δ;
  xl = x - xh;
}

template<typename T, typename U>
constexpr DoublePrecision<Product<T, U>> VeltkampDekkerProduct(T const& a,
                                                               U const& b) {
  DoublePrecision<Product<T, U>> result(uninitialized);
  auto const& x = a;
  auto const& y = b;
  auto& z = result.value;
  auto& zz = result.error;
  // Split x and y as in mul12 from [Dek71, p. 241]; see also [Dek71, p. 235].
  constexpr std::int64_t t = std::numeric_limits<double>::digits;
  constexpr std::int64_t s = t - t / 2;
  T hx;
  T tx;
  VeltkampSplitting<s>(x, hx, tx);
  U hy;
  U ty;
  VeltkampSplitting<s>(y, hy, ty);
  // Veltkamp’s 1968 algorithm, as given in [Dek71, p. 234].
  // See also exactmul in [Lin81, p. 278].
  z = x * y;
  zz = (((hx * hy - z) + hx * ty) + tx * hy) + tx * ty;
  // Dekker’s algorithm (5.12) would be
  // z = (hx * hy) + (hx * ty + tx * hy);
  // zz = (hx * hy) - z + (hx * ty + tx * hy) + tx * ty;
  // where the parenthesized expressions are common (overall one more addition
  // and one multiplication fewer), but it then requires an additional
  // QuickTwoSum (5.14) for z to be the correctly-rounded product x * y,
  // which we require, e.g., for equality comparison.
  return result;
}

template<FMAPresence fma_presence, typename T, typename U>
FORCE_INLINE DoublePrecision<Product<T, U>> TwoProduct(T const& a, U const& b) {
  if (fma_presence == FMAPresence::Present ||
      (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA)) {
    // [Mul+10, p. 152], algorithm 2MultFMA.
    using numerics::_elementary_functions::FusedMultiplySubtract;
    DoublePrecision<Product<T, U>> result(uninitialized);
    result.value = a * b;
    result.error = FusedMultiplySubtract(a, b, result.value);
    return result;
  } else {
    return VeltkampDekkerProduct(a, b);
  }
}

template<typename T, typename U>
FORCE_INLINE constexpr
DoublePrecision<Sum<T, U>> QuickTwoSum(T const& a, U const& b) {
#if _DEBUG
  using quantities::_quantities::DebugString;
  using Comparator = DaumasRideauThéryComparator<T, U>;
  CHECK(Comparator::ProperlyOrdered(a, b))
      << "|" << DebugString(a) << "| < |" << DebugString(b) << "|";
#endif
  // [HLB07], Algorithm 3.
  DoublePrecision<Sum<T, U>> result(uninitialized);
  auto& s = result.value;
  auto& e = result.error;
  s = a + b;
  e = b - (s - a);
  return result;
}

template<typename T, typename U>
FORCE_INLINE constexpr
DoublePrecision<Difference<T, U>> QuickTwoDifference(T const& a, U const& b) {
#if _DEBUG
  using quantities::_quantities::DebugString;
  using Comparator = DaumasRideauThéryComparator<T, U>;
  if (!Comparator::ProperlyOrdered(a, b))
     LOG(FATAL) << "|" << DebugString(a) << "| < |" << DebugString(b) << "|";
#endif
  // [HLB07], Algorithm 3.
  DoublePrecision<Sum<T, U>> result(uninitialized);
  auto& s = result.value;
  auto& e = result.error;
  s = a - b;
  e = (a - s) - b;
  return result;
}

template<typename T, typename U>
constexpr DoublePrecision<Sum<T, U>> TwoSum(T const& a, U const& b) {
  // [HLB07], Algorithm 4.
  DoublePrecision<Sum<T, U>> result(uninitialized);
  auto& s = result.value;
  auto& e = result.error;
  s = a + b;
  auto const v = s - a;
  e = (a - (s - v)) + (b - v);
  return result;
}

// Point × Point → Vector.
template<typename T, typename U, typename>
  requires(!std::is_same_v<U, Difference<U>>)
constexpr DoublePrecision<Difference<T, U>> TwoDifference(T const& a,
                                                          U const& b) {
  // [HLB07], Algorithm 4.
  static_assert(std::is_same_v<T, U>,
                "Template metaprogramming went wrong");
  using Point = T;
  using Vector = Difference<T, U>;
  DoublePrecision<Vector> result(uninitialized);
  Vector& s = result.value;
  Vector& e = result.error;
  s = a - b;
  // Corresponds to -v in `TwoSum`.
  Point const w = a - s;
  e = (a - (s + w)) + (w - b);
  return result;
}

// Point × Vector → Point, or Vector × Vector → Vector.
template<typename T, typename U, typename>
constexpr DoublePrecision<Difference<T, U>> TwoDifference(T const& a,
                                                          U const& b) {
  // [HLB07], Algorithm 4.
  DoublePrecision<Sum<T, U>> result(uninitialized);
  auto& s = result.value;
  auto& e = result.error;
  s = a - b;
  auto const v = s - a;
  e = (a - (s - v)) - (b + v);
  return result;
}

template<typename T>
DoublePrecision<Difference<T>> operator+(DoublePrecision<T> const& left) {
  static_assert(std::is_same_v<Difference<T>, T>,
                "Unary + must be used on a vector");
  return left;
}

template<typename T>
DoublePrecision<Difference<T>> operator-(DoublePrecision<T> const& left) {
  static_assert(std::is_same_v<Difference<T>, T>,
                "Unary - must be used on a vector");
  DoublePrecision<Difference<T>> result(uninitialized);
  result.value = -left.value;
  result.error = -left.error;
  return result;
}

template<typename T, typename U>
DoublePrecision<Sum<T, U>> operator+(T const& left,
                                     DoublePrecision<U> const& right) {
  // [MR22], algorithm 5, relative error 2 u^2.
  auto const s = TwoSum(left, right.value);
  auto const v = s.error + right.error;
  return QuickTwoSum(s.value, v);
}

template<typename T, typename U>
DoublePrecision<Sum<T, U>> operator+(DoublePrecision<T> const& left,
                                     U const& right) {
  // [MR22], algorithm 5, relative error 2 u^2.
  auto const s = TwoSum(left.value, right);
  auto const v = s.error + left.error;
  return QuickTwoSum(s.value, v);
}

template<typename T, typename U>
DoublePrecision<Sum<T, U>> operator+(DoublePrecision<T> const& left,
                                     DoublePrecision<U> const& right) {
  // [MR22], algorithm 4, relative error 3 u^2 + 13 u^3.
  auto const s = TwoSum(left.value, right.value);
  auto const t = TwoSum(left.error, right.error);
  auto const c = s.error + t.value;
  auto const v = QuickTwoSum(s.value, c);
  auto const w = t.error + v.error;
  return QuickTwoSum(v.value, w);
}

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(T const& left,
                                            DoublePrecision<U> const& right) {
  // [MR22], algorithm 5, relative error 2 u^2.
  auto const s = TwoDifference(left, right.value);
  auto const v = s.error - right.error;
  return QuickTwoSum(s.value, v);
}

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(DoublePrecision<T> const& left,
                                            U const& right) {
  // [MR22], algorithm 5, relative error 2 u^2.
  auto const s = TwoDifference(left.value, right);
  auto const v = s.error + left.error;
  return QuickTwoSum(s.value, v);
}

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(DoublePrecision<T> const& left,
                                            DoublePrecision<U> const& right) {
  // [MR22], algorithm 4, relative error 3 u^2 + 13 u^3.
  auto const s = TwoDifference(left.value, right.value);
  auto const t = TwoDifference(left.error, right.error);
  auto const c = s.error + t.value;
  auto const v = QuickTwoSum(s.value, c);
  auto const w = t.error + v.error;
  return QuickTwoSum(v.value, w);
}

template<typename T, typename U>
FORCE_INLINE
DoublePrecision<Product<T, U>> operator*(T const& left,
                                         DoublePrecision<U> const& right) {
  // [Lin81], algorithm longmul.
  auto product = TwoProduct<FMAPresence::Unknown>(left, right.value);
  product.error += left * right.error;
  return QuickTwoSum(product.value, product.error);
}

template<typename T, typename U>
FORCE_INLINE
DoublePrecision<Product<T, U>> operator*(DoublePrecision<T> const& left,
                                         U const& right) {
  // [Lin81], algorithm longmul.
  auto product = TwoProduct<FMAPresence::Unknown>(left.value, right);
  product.error += +left.error * right;
  return QuickTwoSum(product.value, product.error);
}

template<typename T, typename U>
FORCE_INLINE
DoublePrecision<Product<T, U>> operator*(DoublePrecision<T> const& left,
                                         DoublePrecision<U> const& right) {
  // [Lin81], algorithm longmul.
  auto product = TwoProduct<FMAPresence::Unknown>(left.value, right.value);
  product.error +=
      (left.value + left.error) * right.error + left.error * right.value;
  return QuickTwoSum(product.value, product.error);
}

template<typename T, typename U>
DoublePrecision<Quotient<T, U>> operator/(T const& left,
                                          DoublePrecision<U> const& right) {
  // [Lin81], algorithm longdiv.
  auto const z = left / right.value;
  auto const product = TwoProduct<FMAPresence::Unknown>(right.value, z);
  auto const zz = (((left - product.value) - product.error) - z * right.error) /
                  (right.value + right.error);
  return QuickTwoSum(z, zz);
}

template<typename T, typename U>
DoublePrecision<Quotient<T, U>> operator/(DoublePrecision<T> const& left,
                                          U const& right) {
  // [Lin81], algorithm longdiv.
  auto const z = left.value / right;
  auto const product = TwoProduct<FMAPresence::Unknown>(right, z);
  auto const zz =
      (((left.value - product.value) - product.error) + left.error) / right;
  return QuickTwoSum(z, zz);
}

template<typename T, typename U>
DoublePrecision<Quotient<T, U>> operator/(DoublePrecision<T> const& left,
                                          DoublePrecision<U> const& right) {
  // [Lin81], algorithm longdiv.
  auto const z = left.value / right.value;
  auto const product = TwoProduct<FMAPresence::Unknown>(right.value, z);
  auto const zz =
      ((((left.value - product.value) - product.error) + left.error) -
       z * right.error) /
      (right.value + right.error);
  return QuickTwoSum(z, zz);
}

template<typename T>
bool IsFinite(DoublePrecision<T> const& double_precision) {
  return IsFinite(double_precision.value) && IsFinite(double_precision.error);
}

template<typename T>
std::string DebugString(DoublePrecision<T> const& double_precision) {
  // We use `DebugString` to get all digits when `T` is `double`.  In that case
  // ADL will not find it, so we need the `using`.  For some values of `T`,
  // `DebugString` will come from elsewhere, so we cannot directly call
  // `quantities::Multivector`.
  using quantities::_quantities::DebugString;
  return DebugString(double_precision.value) + "|" +
         DebugString(double_precision.error);
}

template<typename T>
std::ostream& operator<<(std::ostream& os,
                         const DoublePrecision<T>& double_precision) {
  os << DebugString(double_precision);
  return os;
}

}  // namespace internal
}  // namespace _double_precision
}  // namespace numerics
}  // namespace principia
