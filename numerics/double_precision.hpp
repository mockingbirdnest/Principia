
#pragma once

#include <string>

#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
namespace internal_double_precision {

using base::not_null;
using quantities::Angle;
using quantities::Difference;
using quantities::Product;
using quantities::Quotient;
using quantities::Sum;

// A simple container for accumulating a value using double precision.  The
// type of the value must be an affine space.  The notations follow [HLB08].
template<typename T>
struct DoublePrecision final {
  constexpr DoublePrecision() = default;

  explicit constexpr DoublePrecision(T const& value);

  // Compensated summation.  This is less precise, but more efficient, than
  // |operator-=| or |operator+=|.  Unlike |QuickTwoSum|, these functions don't
  // DCHECK their argument, so the caller must ensure that |right| is small
  // enough.
  DoublePrecision<T>& Decrement(Difference<T> const& right);
  DoublePrecision<T>& Increment(Difference<T> const& right);

  DoublePrecision<T>& operator+=(DoublePrecision<Difference<T>> const& right);
  DoublePrecision<T>& operator+=(Difference<T> const& right);
  DoublePrecision<T>& operator-=(DoublePrecision<Difference<T>> const& right);
  DoublePrecision<T>& operator-=(Difference<T> const& right);

  void WriteToMessage(not_null<serialization::DoublePrecision*> message) const;
  static DoublePrecision ReadFromMessage(
      serialization::DoublePrecision const& message);

  T value{};
  Difference<T> error{};
};

// |scale| must be a signed power of two or zero.
template<typename T, typename U>
DoublePrecision<Product<T, U>> Scale(T const& scale,
                                     DoublePrecision<U> const& right);

// Returns the exact product of its arguments.  Note that this function checks whether
// |UseHardwareFMA| is true.  If the value of that flag is already known from
// context, it may be preferable to either:
// — use VeltkampDekkerProduct(a, b) below;
// — directly compute value = a * b, error = FusedMultiplySubtract(a, b, value).
template<typename T, typename U>
DoublePrecision<Product<T, U>> TwoProduct(T const& a, U const& b);

// Same as |TwoProduct|, but never uses FMA.
template<typename T, typename U>
DoublePrecision<Product<T, U>> VeltkampDekkerProduct(T const& a, U const& b);

// Computes the exact sum of a and b.  The arguments must be such that
// |a| >= |b| or a == 0.
template<typename T, typename U>
DoublePrecision<Sum<T, U>> QuickTwoSum(T const& a, U const& b);

// Computes the exact sum of a and b.
template<typename T, typename U>
DoublePrecision<Sum<T, U>> TwoSum(T const& a, U const& b);

// |TwoDifference| may have any of the following signatures:
//   1. Point × Point → Vector;
//   2. Point × Vector → Point;
//   3. Vector × Vector → Vector;
// The first overload handles the first case, and the second handles the last
// two.
template<typename T,
         typename U,
         typename = Difference<T, Difference<T, U>>,
         typename = std::enable_if_t<!std::is_same<U, Difference<U>>::value>>
DoublePrecision<Difference<T, U>> TwoDifference(T const& a, U const& b);

template<typename T, typename U, typename = Difference<Difference<T, U>, T>>
DoublePrecision<Difference<T, U>> TwoDifference(T const& a, U const& b);

DoublePrecision<Angle> Mod2π(DoublePrecision<Angle> const& θ);

template<typename T>
bool operator==(DoublePrecision<T> const& left,
                DoublePrecision<T> const& right);

template<typename T>
bool operator!=(DoublePrecision<T> const& left,
                DoublePrecision<T> const& right);

// |T| must be a vector.
template<typename T>
DoublePrecision<Difference<T>> operator+(DoublePrecision<T> const& left);

// |T| must be a vector.
template<typename T>
DoublePrecision<Difference<T>> operator-(DoublePrecision<T> const& left);

template<typename T, typename U>
DoublePrecision<Sum<T, U>> operator+(DoublePrecision<T> const& left,
                                     DoublePrecision<U> const& right);

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(DoublePrecision<T> const& left,
                                            DoublePrecision<U> const& right);

template<typename T, typename U>
DoublePrecision<Product<T, U>> operator*(DoublePrecision<T> const& left,
                                         DoublePrecision<U> const& right);

template<typename T, typename U>
DoublePrecision<Quotient<T, U>> operator/(DoublePrecision<T> const& left,
                                          DoublePrecision<U> const& right);

template<typename T>
std::string DebugString(DoublePrecision<T> const& double_precision);

template<typename T>
std::ostream& operator<<(std::ostream& os,
                         DoublePrecision<T> const& double_precision);

}  // namespace internal_double_precision

using internal_double_precision::DoublePrecision;
using internal_double_precision::Mod2π;
using internal_double_precision::TwoDifference;
using internal_double_precision::TwoProduct;
using internal_double_precision::TwoSum;
using internal_double_precision::VeltkampDekkerProduct;

}  // namespace numerics
}  // namespace principia

#include "numerics/double_precision_body.hpp"
