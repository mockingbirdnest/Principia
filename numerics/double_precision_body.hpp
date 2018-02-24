
#pragma once

#include "numerics/double_precision.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <string>

#include "geometry/serialization.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace internal_double_precision {

using geometry::DoubleOrQuantityOrPointOrMultivectorSerializer;
using geometry::DoubleOrQuantityOrMultivectorSerializer;
using quantities::Abs;
using quantities::FusedMultiplyAdd;
using quantities::Quantity;
using quantities::SIUnit;

// Assumes that |T| and |U| have a memory representation that is a sequence of
// |double|s, and returns the conjunction of componentwise
// |left[i]| >= |right[i]| or left[i] == 0.
template<typename T, typename U>
bool ComponentwiseGreaterThanOrEqualOrZero(T const& left, U const& right) {
  static_assert(sizeof(left) == sizeof(right),
                "Comparing types of different sizes");
  static_assert(sizeof(left) % sizeof(double) == 0,
                "Types are not sequences of doubles");
  constexpr int size = sizeof(left) / sizeof(double);
  std::array<double, size> left_doubles;
  std::memcpy(left_doubles.data(), &left, sizeof(left));
  std::array<double, size> right_doubles;
  std::memcpy(right_doubles.data(), &right, sizeof(right));
  bool result = true;
  for (int i = 0; i < size; ++i) {
    result &=
        (Abs(left_doubles[i]) >= Abs(right_doubles[i]) || left_doubles[i] == 0);
  }
  return result;
}

template<typename T>
constexpr DoublePrecision<T>::DoublePrecision(T const& value)
    : value(value),
      error() {}

template <typename T>
DoublePrecision<T>& DoublePrecision<T>::Decrement(Difference<T> const& right) {
  // See Higham, Accuracy and Stability of Numerical Algorithms, Algorithm 4.2.
  // This is equivalent to |QuickTwoSum(value, error - right)|.
  T const temp = value;
  Difference<T> const y = error - right;
  value = temp + y;
  error = (temp - value) + y;
  return *this;
}

template <typename T>
DoublePrecision<T>& DoublePrecision<T>::Increment(Difference<T> const& right) {
  // See Higham, Accuracy and Stability of Numerical Algorithms, Algorithm 4.2.
  // This is equivalent to |QuickTwoSum(value, error + right)|.
  T const temp = value;
  Difference<T> const y = error + right;
  value = temp + y;
  error = (temp - value) + y;
  return *this;
}

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
  double const s = scale / quantities::SIUnit<T>();
  if (s != 0.0) {
    int exponent;
    double const mantissa = std::frexp(s, &exponent);
    CHECK_EQ(0.5, std::fabs(mantissa)) << scale;
  }
#endif
  DoublePrecision<Product<T, U>> result;
  result.value = right.value * scale;
  result.error = right.error * scale;
  return result;
}

template<typename T, typename U>
DoublePrecision<Product<T, U>> TwoProduct(T const& a, U const& b) {
  DoublePrecision<Product<T, U>> result(a * b);
  result.error = FusedMultiplyAdd(a, b, -result.value);
  return result;
}

template<typename T, typename U>
FORCE_INLINE(inline)
DoublePrecision<Sum<T, U>> QuickTwoSum(T const& a, U const& b) {
#if _DEBUG
  using quantities::DebugString;
  CHECK(ComponentwiseGreaterThanOrEqualOrZero(a, b))
      << "|" << DebugString(a) << "| < |" << DebugString(b) << "|";
#endif
  // Hida, Li and Bailey (2007), Library for Double-Double and Quad-Double
  // Arithmetic.
  DoublePrecision<Sum<T, U>> result;
  auto& s = result.value;
  auto& e = result.error;
  s = a + b;
  e = b - (s - a);
  return result;
}

template<typename T, typename U>
DoublePrecision<Sum<T, U>> TwoSum(T const& a, U const& b) {
  // Hida, Li and Bailey (2007), Library for Double-Double and Quad-Double
  // Arithmetic.
  DoublePrecision<Sum<T, U>> result;
  auto& s = result.value;
  auto& e = result.error;
  s = a + b;
  auto const v = s - a;
  e = (a - (s - v)) + (b - v);
  return result;
}

// Point × Point → Vector.
template<typename T, typename U, typename, typename>
DoublePrecision<Difference<T, U>> TwoDifference(T const& a, U const& b) {
  static_assert(std::is_same<T, U>::value,
                "Template metaprogramming went wrong");
  using Point = T;
  using Vector = Difference<T, U>;
  DoublePrecision<Vector> result;
  Vector& s = result.value;
  Vector& e = result.error;
  s = a - b;
  // Corresponds to -v in |TwoSum|.
  Point const w = a - s;
  e = (a - (s + w)) + (w - b);
  return result;
}

// Point × Vector → Point, or Vector × Vector → Vector.
template<typename T, typename U, typename>
DoublePrecision<Difference<T, U>> TwoDifference(T const& a, U const& b) {
  return TwoSum(a, -b);
}

template<typename T>
bool operator==(DoublePrecision<T> const& left,
                DoublePrecision<T> const& right) {
  // This is correct assuming that left and right have non-overlapping
  // mantissas.
  return left.value == right.value && left.error == right.error;
}

template<typename T>
bool operator!=(DoublePrecision<T> const& left,
                DoublePrecision<T> const& right) {
  // This is correct assuming that left and right have non-overlapping
  // mantissas.
  return left.value != right.value || left.error != right.error;
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
DoublePrecision<Sum<T, U>> operator+(DoublePrecision<T> const& left,
                                     DoublePrecision<U> const& right) {
  // Linnainmaa (1981), Software for Doubled-Precision Floating-Point
  // Computations, algorithm longadd.
  auto const sum = TwoSum(left.value, right.value);
  return QuickTwoSum(sum.value, (sum.error + left.error) + right.error);
}

template<typename T, typename U>
DoublePrecision<Sum<T, U>> operator+(DoublePrecision<T> const& left,
                                     U const& right) {
  auto const sum = TwoSum(left.value, right);
  return QuickTwoSum(sum.value, sum.error + left.error);
}

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(DoublePrecision<T> const& left,
                                            DoublePrecision<U> const& right) {
  // Linnainmaa (1981), Software for Doubled-Precision Floating-Point
  // Computations, algorithm longadd.
  auto const sum = TwoDifference(left.value, right.value);
  return QuickTwoSum(sum.value, (sum.error + left.error) - right.error);
}

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(DoublePrecision<T> const& left,
                                            U const& right) {
  auto const sum = TwoDifference(left.value, right);
  return QuickTwoSum(sum.value, sum.error + left.error);
}

template<typename T>
std::string DebugString(DoublePrecision<T> const& double_precision) {
  // We use |DebugString| to get all digits when |T| is |double|.  In that case
  // ADL will not find it, so we need the |using|.  For some values of |T|,
  // |DebugString| will come from elsewhere, so we cannot directly call
  // |quantities::Multivector|.
  using quantities::DebugString;
  return DebugString(double_precision.value) + "|" +
         DebugString(double_precision.error);
}

template<typename T>
std::ostream& operator<<(std::ostream& os,
                         const DoublePrecision<T>& double_precision) {
  os << DebugString(double_precision);
  return os;
}

}  // namespace internal_double_precision
}  // namespace numerics
}  // namespace principia
