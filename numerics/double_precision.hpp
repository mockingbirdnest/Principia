#pragma once

#include <string>

#include "base/not_null.hpp"
#include "numerics/fma.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
namespace _double_precision {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::numerics::_fma;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// A simple container for accumulating a value using double precision.  The
// type of the value must be an affine space.  The notations follow [HLB08].
template<typename T>
struct DoublePrecision final {
  constexpr DoublePrecision() = default;

  explicit constexpr DoublePrecision(T const& value);

  // This is correct assuming that left and right have non-overlapping
  // mantissas.
  friend auto operator<=>(DoublePrecision const& left,
                          DoublePrecision const& right) = default;

  DoublePrecision<T>& operator+=(DoublePrecision<Difference<T>> const& right);
  DoublePrecision<T>& operator+=(Difference<T> const& right);
  DoublePrecision<T>& operator-=(DoublePrecision<Difference<T>> const& right);
  DoublePrecision<T>& operator-=(Difference<T> const& right);

  // Compensated summation.  This is less precise, but more efficient, than
  // `operator-=` or `operator+=`.  Unlike `QuickTwoSum`, these functions don't
  // DCHECK their argument, so the caller must ensure that `right` is small
  // enough.
  DoublePrecision<T>& Decrement(Difference<T> const& right);
  DoublePrecision<T>& Increment(Difference<T> const& right);

  void WriteToMessage(not_null<serialization::DoublePrecision*> message) const;
  static DoublePrecision ReadFromMessage(
      serialization::DoublePrecision const& message);

  T value{};
  Difference<T> error{};
};

// `scale` must be a signed power of two or zero.
template<typename T, typename U>
DoublePrecision<Product<T, U>> Scale(T const& scale,
                                     DoublePrecision<U> const& right);

// Returns the exact product of its arguments.  Note that this function checks
// whether `UseHardwareFMA` is true.  If the value of that flag is already known
// from context, it may be preferable to either:
// — use VeltkampDekkerProduct(a, b) below;
// — directly compute value = a * b, error = FusedMultiplySubtract(a, b, value).
template<FMAPolicy fma_policy = FMAPolicy::Auto, typename T, typename U>
DoublePrecision<Product<T, U>> TwoProduct(T const& a, U const& b);

// Returns the exact value of `a * b + c` if `|a * b|` is small compared to
// `|c|`. See [SZ05], section 2.1.
template<FMAPolicy fma_policy = FMAPolicy::Auto, typename T, typename U>
DoublePrecision<Product<T, U>> TwoProductAdd(T const& a,
                                             U const& b,
                                             Product<T, U> const& c);

// Returns the exact value of `a * b - c` if `|a * b|` is small compared to
// `|c|`. See [SZ05], section 2.1.
template<FMAPolicy fma_policy = FMAPolicy::Auto, typename T, typename U>
DoublePrecision<Product<T, U>> TwoProductSubtract(T const& a,
                                                  U const& b,
                                                  Product<T, U> const& c);

// Returns the exact value of `-a * b + c` if `|a * b|` is small compared to
// `|c|`. See [SZ05], section 2.1.
template<FMAPolicy fma_policy = FMAPolicy::Auto, typename T, typename U>
DoublePrecision<Product<T, U>> TwoProductNegatedAdd(T const& a,
                                                    U const& b,
                                                    Product<T, U> const& c);

// Returns the exact value of `-a * b - c` if `|a * b|` is small compared to
// `|c|`. See [SZ05], section 2.1.
template<FMAPolicy fma_policy = FMAPolicy::Auto, typename T, typename U>
DoublePrecision<Product<T, U>>
TwoProductNegatedSubtract(T const& a,
                          U const& b,
                          Product<T, U> const& c);

// Same as `TwoProduct`, but never uses FMA.
template<typename T, typename U>
constexpr DoublePrecision<Product<T, U>> VeltkampDekkerProduct(T const& a,
                                                               U const& b);

// Computes the exact sum of a and b.  The arguments must be such that
// |a| >= |b| or a == 0.
template<typename T, typename U>
constexpr DoublePrecision<Sum<T, U>> QuickTwoSum(T const& a, U const& b);

// Computes the exact sum of a and b.
template<typename T, typename U>
constexpr DoublePrecision<Sum<T, U>> TwoSum(T const& a, U const& b);

// `TwoDifference` may have any of the following signatures:
//   1. Point × Point → Vector;
//   2. Point × Vector → Point;
//   3. Vector × Vector → Vector;
// The first overload handles the first case, and the second handles the last
// two.
template<typename T,
         typename U,
         typename = Difference<T, Difference<T, U>>,
         typename = std::enable_if_t<!std::is_same<U, Difference<U>>::value>>
constexpr DoublePrecision<Difference<T, U>> TwoDifference(T const& a,
                                                          U const& b);

template<typename T, typename U, typename = Difference<Difference<T, U>, T>>
constexpr DoublePrecision<Difference<T, U>> TwoDifference(T const& a,
                                                          U const& b);

// `T` must be a vector.
template<typename T>
DoublePrecision<Difference<T>> operator+(DoublePrecision<T> const& left);

// `T` must be a vector.
template<typename T>
DoublePrecision<Difference<T>> operator-(DoublePrecision<T> const& left);

template<typename T, typename U>
DoublePrecision<Sum<T, U>> operator+(T const& left,
                                     DoublePrecision<U> const& right);

template<typename T, typename U>
DoublePrecision<Sum<T, U>> operator+(DoublePrecision<T> const& left,
                                     U const& right);

template<typename T, typename U>
DoublePrecision<Sum<T, U>> operator+(DoublePrecision<T> const& left,
                                     DoublePrecision<U> const& right);

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(T const& left,
                                            DoublePrecision<U> const& right);

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(DoublePrecision<T> const& left,
                                            U const& right);

template<typename T, typename U>
DoublePrecision<Difference<T, U>> operator-(DoublePrecision<T> const& left,
                                            DoublePrecision<U> const& right);

template<typename T, typename U>
DoublePrecision<Product<T, U>> operator*(T const& left,
                                         DoublePrecision<U> const& right);

template<typename T, typename U>
DoublePrecision<Product<T, U>> operator*(DoublePrecision<T> const& left,
                                         U const& right);

template<typename T, typename U>
DoublePrecision<Product<T, U>> operator*(DoublePrecision<T> const& left,
                                         DoublePrecision<U> const& right);

template<typename T, typename U>
DoublePrecision<Quotient<T, U>> operator/(T const& left,
                                          DoublePrecision<U> const& right);

template<typename T, typename U>
DoublePrecision<Quotient<T, U>> operator/(DoublePrecision<T> const& left,
                                          U const& right);

template<typename T, typename U>
DoublePrecision<Quotient<T, U>> operator/(DoublePrecision<T> const& left,
                                          DoublePrecision<U> const& right);

template<typename T>
std::string DebugString(DoublePrecision<T> const& double_precision);

template<typename T>
std::ostream& operator<<(std::ostream& os,
                         DoublePrecision<T> const& double_precision);

}  // namespace internal

using internal::DoublePrecision;
using internal::QuickTwoSum;
using internal::TwoDifference;
using internal::TwoProduct;
using internal::TwoProductAdd;
using internal::TwoProductNegatedAdd;
using internal::TwoProductNegatedSubtract;
using internal::TwoProductSubtract;
using internal::TwoSum;
using internal::VeltkampDekkerProduct;

}  // namespace _double_precision
}  // namespace numerics
}  // namespace principia

#include "numerics/double_precision_body.hpp"
