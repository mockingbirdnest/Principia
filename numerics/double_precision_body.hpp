#pragma once

#include "numerics/double_precision.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>

#include "base/not_constructible.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/serialization.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/concepts.hpp"
#include "quantities/m128d.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _double_precision {
namespace internal {

using namespace principia::base::_not_constructible;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_point;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_serialization;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_m128d;
using namespace principia::quantities::_si;

// A helper to check that the preconditions of QuickTwoSum are met.  Annoyingly
// complicated as it needs to peel off all of our abstractions until it reaches
// doubles.
template<typename T, typename U, typename = void>
struct ComponentwiseComparator;

template<typename T, typename U>
struct ComponentwiseComparator<DoublePrecision<T>, DoublePrecision<U>> {
  static bool GreaterThanOrEqualOrZero(DoublePrecision<T> const& left,
                                       DoublePrecision<U> const& right) {
    // This check is incomplete: it doesn't compare the errors component.  To
    // the best of my understanding this code is only used in tests.
    return ComponentwiseComparator<T, U>::GreaterThanOrEqualOrZero(left.value,
                                                                   right.value);
  }
};

template<typename T, typename U>
struct ComponentwiseComparator<Point<T>, U> : not_constructible {
  static bool GreaterThanOrEqualOrZero(Point<T> const& left,
                                       U const& right) {
    // We only care about the coordinates, the geometric structure is
    // irrelevant.
    return ComponentwiseComparator<T, U>::GreaterThanOrEqualOrZero(
        left - Point<T>{}, right);
  }
};

template<typename T, typename U>
struct ComponentwiseComparator<T, Point<U>> : not_constructible {
  static bool GreaterThanOrEqualOrZero(T const& left,
                                       Point<U> const& right) {
    // We only care about the coordinates, the geometric structure is
    // irrelevant.
    return ComponentwiseComparator<T, U>::GreaterThanOrEqualOrZero(
        left, right - Point<U>{});
  }
};

template<typename T, typename TFrame, int trank,
         typename U, typename UFrame, int urank>
struct ComponentwiseComparator<Multivector<T, TFrame, trank>,
                               Multivector<U, UFrame, urank>>
    : not_constructible {
  static bool GreaterThanOrEqualOrZero(
      Multivector<T, TFrame, trank> const& left,
      Multivector<U, UFrame, urank> const& right) {
    // This doesn't handle trivectors.
    return ComponentwiseComparator<R3Element<T>, R3Element<U>>::
        GreaterThanOrEqualOrZero(left.coordinates(), right.coordinates());
  }
};

template<typename T, typename U>
struct ComponentwiseComparator<R3Element<T>, R3Element<U>> : not_constructible {
  static bool GreaterThanOrEqualOrZero(R3Element<T> const& left,
                                       R3Element<U> const& right) {
    bool result = true;
    for (int i = 0; i < 3; ++i) {
      result &= ComponentwiseComparator<T, U>::GreaterThanOrEqualOrZero(
          left[i], right[i]);
    }
    return result;
  }
};

template<>
struct ComponentwiseComparator<M128D, M128D> {
  static bool GreaterThanOrEqualOrZero(M128D const& left, M128D const& right) {
    static M128D const zero(0.0);
    return Abs(left) >= Abs(right) || left == zero ||
           left != left || right != right;
  }
};

template<typename T, typename U>
  requires convertible_to_quantity<T> && convertible_to_quantity<U>
struct ComponentwiseComparator<T, U> {
  static bool GreaterThanOrEqualOrZero(T const& left, U const& right) {
    // In the elementary functions, we use NaN to fall back to the CORE-MATH
    // implementation.  We don't want to die because of the weird comparisons of
    // NaNs.
    return Abs(left) >= Abs(right) || left == T{} ||
           !IsFinite(left) || !IsFinite(right);
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
  *this = *this + DoublePrecision(right);
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
  *this = *this - DoublePrecision(right);
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
  constexpr double c = 1 << (t - t / 2);
  T const px = x * c;
  T const hx = x - px + px;
  T const tx = x - hx;
  U const py = y * c;
  U const hy = y - py + py;
  U const ty = y - hy;
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

template<FMAPolicy fma_policy, typename T, typename U>
FORCE_INLINE(inline)
DoublePrecision<Product<T, U>> TwoProduct(T const& a, U const& b) {
  if ((fma_policy == FMAPolicy::Force && CanEmitFMAInstructions) ||
      (fma_policy == FMAPolicy::Auto && UseHardwareFMA)) {
    using numerics::_elementary_functions::FusedMultiplySubtract;
    DoublePrecision<Product<T, U>> result(uninitialized);
    result.value = a * b;
    result.error = FusedMultiplySubtract(a, b, result.value);
    return result;
  } else {
    return VeltkampDekkerProduct(a, b);
  }
}

template<FMAPolicy fma_policy, typename T, typename U>
FORCE_INLINE(inline)
DoublePrecision<Product<T, U>> TwoProductAdd(T const& a,
                                             U const& b,
                                             Product<T, U> const& c) {
  if ((fma_policy == FMAPolicy::Force && CanEmitFMAInstructions) ||
      (fma_policy == FMAPolicy::Auto && UseHardwareFMA)) {
    using numerics::_elementary_functions::FusedMultiplyAdd;
    using numerics::_elementary_functions::FusedMultiplySubtract;
    DoublePrecision<Product<T, U>> result(uninitialized);
    result.value = FusedMultiplyAdd(a, b, c);
    result.error = FusedMultiplySubtract(a, b, result.value - c);
    return result;
  } else {
    auto result = VeltkampDekkerProduct(a, b);
    result += c;
    return result;
  }
}

template<FMAPolicy fma_policy, typename T, typename U>
FORCE_INLINE(inline)
DoublePrecision<Product<T, U>> TwoProductSubtract(T const& a,
                                                  U const& b,
                                                  Product<T, U> const& c) {
  if ((fma_policy == FMAPolicy::Force && CanEmitFMAInstructions) ||
      (fma_policy == FMAPolicy::Auto && UseHardwareFMA)) {
    using numerics::_elementary_functions::FusedMultiplySubtract;
    DoublePrecision<Product<T, U>> result(uninitialized);
    result.value = FusedMultiplySubtract(a, b, c);
    result.error = FusedMultiplySubtract(a, b, result.value + c);
    return result;
  } else {
    auto result = VeltkampDekkerProduct(a, b);
    result -= c;
    return result;
  }
}

template<FMAPolicy fma_policy, typename T, typename U>
FORCE_INLINE(inline)
DoublePrecision<Product<T, U>> TwoProductNegatedAdd(T const& a,
                                                    U const& b,
                                                    Product<T, U> const& c) {
  if ((fma_policy == FMAPolicy::Force && CanEmitFMAInstructions) ||
      (fma_policy == FMAPolicy::Auto && UseHardwareFMA)) {
    using numerics::_elementary_functions::FusedNegatedMultiplyAdd;
    using numerics::_elementary_functions::FusedNegatedMultiplySubtract;
    DoublePrecision<Product<T, U>> result(uninitialized);
    result.value = FusedNegatedMultiplyAdd(a, b, c);
    result.error = FusedNegatedMultiplySubtract(a, b, result.value - c);
    return result;
  } else {
    auto result = VeltkampDekkerProduct(-a, b);
    result += c;
    return result;
  }
}

template<FMAPolicy fma_policy, typename T, typename U>
FORCE_INLINE(inline)
DoublePrecision<Product<T, U>>
TwoProductNegatedSubtract(T const& a,
                          U const& b,
                          Product<T, U> const& c) {
  if ((fma_policy == FMAPolicy::Force && CanEmitFMAInstructions) ||
      (fma_policy == FMAPolicy::Auto && UseHardwareFMA)) {
    using numerics::_elementary_functions::FusedNegatedMultiplySubtract;
    DoublePrecision<Product<T, U>> result(uninitialized);
    result.value = FusedNegatedMultiplySubtract(a, b, c);
    result.error = FusedNegatedMultiplySubtract(a, b, result.value + c);
    return result;
  } else {
    auto result = VeltkampDekkerProduct(-a, b);
    result -= c;
    return result;
  }
}

template<typename T, typename U>
FORCE_INLINE(constexpr)
DoublePrecision<Sum<T, U>> QuickTwoSum(T const& a, U const& b) {
#if _DEBUG
  using quantities::_quantities::DebugString;
  using Comparator = ComponentwiseComparator<T, U>;
  CHECK(Comparator::GreaterThanOrEqualOrZero(a, b))
      << "|" << DebugString(a) << "| < |" << DebugString(b) << "|";
#endif
  // [HLB07].
  DoublePrecision<Sum<T, U>> result{uninitialized};
  auto& s = result.value;
  auto& e = result.error;
  s = a + b;
  e = b - (s - a);
  return result;
}

template<typename T, typename U>
FORCE_INLINE(constexpr)
DoublePrecision<Difference<T, U>> QuickTwoDifference(T const& a, U const& b) {
#if _DEBUG
  using quantities::_quantities::DebugString;
  using Comparator = ComponentwiseComparator<T, U>;
  CHECK(Comparator::GreaterThanOrEqualOrZero(a, b))
      << "|" << DebugString(a) << "| < |" << DebugString(b) << "|";
#endif
  // [HLB07].
  DoublePrecision<Sum<T, U>> result{uninitialized};
  auto& s = result.value;
  auto& e = result.error;
  s = a - b;
  e = (a - s) - b;
  return result;
}

template<typename T, typename U>
constexpr DoublePrecision<Sum<T, U>> TwoSum(T const& a, U const& b) {
  // [HLB07].
  DoublePrecision<Sum<T, U>> result{uninitialized};
  auto& s = result.value;
  auto& e = result.error;
  s = a + b;
  auto const v = s - a;
  e = (a - (s - v)) + (b - v);
  return result;
}

// Point × Point → Vector.
template<typename T, typename U, typename, typename>
constexpr DoublePrecision<Difference<T, U>> TwoDifference(T const& a,
                                                          U const& b) {
  static_assert(std::is_same<T, U>::value,
                "Template metaprogramming went wrong");
  using Point = T;
  using Vector = Difference<T, U>;
  DoublePrecision<Vector> result{uninitialized};
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
  DoublePrecision<Sum<T, U>> result{uninitialized};
  auto& s = result.value;
  auto& e = result.error;
  s = a - b;
  auto const v = s - a;
  e = (a - (s - v)) - (b + v);
  return result;
}

template<typename T>
DoublePrecision<Difference<T>> operator+(DoublePrecision<T> const& left) {
  static_assert(std::is_same<Difference<T>, T>::value,
                "Unary + must be used on a vector");
  return left;
}

template<typename T>
DoublePrecision<Difference<T>> operator-(DoublePrecision<T> const& left) {
  static_assert(std::is_same<Difference<T>, T>::value,
                "Unary - must be used on a vector");
  DoublePrecision<Difference<T>> result;
  result.value = -left.value;
  result.error = -left.error;
  return result;
}

template<typename T, typename U>
DoublePrecision<Sum<T, U>> operator+(T const& left,
                                     DoublePrecision<U> const& right) {
  // [Lin81], algorithm longadd.
  auto const sum = TwoSum(left, right.value);
  return QuickTwoSum(sum.value, sum.error + right.error);
}

template<typename T, typename U>
DoublePrecision<Sum<T, U>> operator+(DoublePrecision<T> const& left,
                                     U const& right) {
  // [Lin81], algorithm longadd.
  auto const sum = TwoSum(left.value, right);
  return QuickTwoSum(sum.value, sum.error + left.error);
}

template<typename T, typename U>
DoublePrecision<Sum<T, U>> operator+(DoublePrecision<T> const& left,
                                     DoublePrecision<U> const& right) {
  // [Lin81], algorithm longadd.
  auto const sum = TwoSum(left.value, right.value);
  return QuickTwoSum(sum.value, (sum.error + left.error) + right.error);
}

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(T const& left,
                                            DoublePrecision<U> const& right) {
  // [Lin81], algorithm longadd.
  auto const sum = TwoDifference(left, right.value);
  return QuickTwoSum(sum.value, sum.error - right.error);
}

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(DoublePrecision<T> const& left,
                                            U const& right) {
  // [Lin81], algorithm longadd.
  auto const sum = TwoDifference(left.value, right);
  return QuickTwoSum(sum.value, sum.error + left.error);
}

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(DoublePrecision<T> const& left,
                                            DoublePrecision<U> const& right) {
  // [Lin81], algorithm longadd.
  auto const sum = TwoDifference(left.value, right.value);
  return QuickTwoSum(sum.value, (sum.error + left.error) - right.error);
}

template<typename T, typename U>
FORCE_INLINE(inline)
DoublePrecision<Product<T, U>> operator*(T const& left,
                                         DoublePrecision<U> const& right) {
  // [Lin81], algorithm longmul.
  auto product = TwoProduct(left, right.value);
  product.error += left * right.error;
  return QuickTwoSum(product.value, product.error);
}

template<typename T, typename U>
FORCE_INLINE(inline)
DoublePrecision<Product<T, U>> operator*(DoublePrecision<T> const& left,
                                         U const& right) {
  // [Lin81], algorithm longmul.
  auto product = TwoProduct(left.value, right);
  product.error += +left.error * right;
  return QuickTwoSum(product.value, product.error);
}

template<typename T, typename U>
FORCE_INLINE(inline)
DoublePrecision<Product<T, U>> operator*(DoublePrecision<T> const& left,
                                         DoublePrecision<U> const& right) {
  // [Lin81], algorithm longmul.
  auto product = TwoProduct(left.value, right.value);
  product.error +=
      (left.value + left.error) * right.error + left.error * right.value;
  return QuickTwoSum(product.value, product.error);
}

template<typename T, typename U>
DoublePrecision<Quotient<T, U>> operator/(T const& left,
                                          DoublePrecision<U> const& right) {
  // [Lin81], algorithm longdiv.
  auto const z = left / right.value;
  auto const product = TwoProduct(right.value, z);
  auto const zz = (((left - product.value) - product.error) - z * right.error) /
                  (right.value + right.error);
  return QuickTwoSum(z, zz);
}

template<typename T, typename U>
DoublePrecision<Quotient<T, U>> operator/(DoublePrecision<T> const& left,
                                          U const& right) {
  // [Lin81], algorithm longdiv.
  auto const z = left.value / right;
  auto const product = TwoProduct(right, z);
  auto const zz =
      (((left.value - product.value) - product.error) + left.error) / right;
  return QuickTwoSum(z, zz);
}

template<typename T, typename U>
DoublePrecision<Quotient<T, U>> operator/(DoublePrecision<T> const& left,
                                          DoublePrecision<U> const& right) {
  // [Lin81], algorithm longdiv.
  auto const z = left.value / right.value;
  auto const product = TwoProduct(right.value, z);
  auto const zz =
      ((((left.value - product.value) - product.error) + left.error) -
       z * right.error) /
      (right.value + right.error);
  return QuickTwoSum(z, zz);
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
