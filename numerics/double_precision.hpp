
#pragma once

#include "quantities/quantities.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
namespace internal_double_precision {

using base::not_null;
using quantities::Difference;
using quantities::Product;
using quantities::Sum;

// A simple container for accumulating a value using double precision.  The
// type of the value must be an affine space.  The notations follow
// Library for Double-Double and Quad-Double Arithmetic, Hida, Li and Bailey,
// 2007.
template<typename T>
struct DoublePrecision final {
  constexpr DoublePrecision() = default;

  // This constructor is not explicit to make it easy to construct an object
  // with no error.
  constexpr DoublePrecision(T const& value);  // NOLINT(runtime/explicit)

  DoublePrecision<T>& operator+=(Difference<T> const& right);
  DoublePrecision<T>& operator+=(DoublePrecision<Difference<T>> const& right);
  DoublePrecision<T>& operator-=(Difference<T> const& right);
  DoublePrecision<T>& operator-=(DoublePrecision<Difference<T>> const& right);

  void WriteToMessage(not_null<serialization::DoublePrecision*> message) const;
  static DoublePrecision ReadFromMessage(
      serialization::DoublePrecision const& message);

  T value;
  Difference<T> error;
};

// |scale| must be a signed power of two or zero.
template<typename T, typename U>
DoublePrecision<Product<T, U>> Scale(T const& scale,
                                     DoublePrecision<U> const& right);

// The arguments must be such that |a| >= |b|.
template<typename T, typename U>
DoublePrecision<Sum<T, U>> QuickTwoSum(T const& a, U const& b);

template<typename T, typename U>
DoublePrecision<Sum<T, U>> TwoSum(T const& a, U const& b);

template<typename T>
DoublePrecision<Difference<T>> operator+(
    DoublePrecision<Difference<T>> const& left);

template<typename T>
DoublePrecision<Difference<T>> operator-(
    DoublePrecision<Difference<T>> const& left);

template<typename T, typename U>
DoublePrecision<Sum<T, U>> operator+(DoublePrecision<T> const& left,
                                     DoublePrecision<U> const& right);

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(DoublePrecision<T> const& left,
                                            DoublePrecision<U> const& right);

template<typename T>
std::ostream& operator<<(std::ostream& os,
                         const DoublePrecision<T>& double_precision);

}  // namespace internal_double_precision

using internal_double_precision::DoublePrecision;
using internal_double_precision::TwoSum;

}  // namespace numerics
}  // namespace principia

#include "numerics/double_precision_body.hpp"
