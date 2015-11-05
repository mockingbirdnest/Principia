#pragma once

// We use ostream for logging purposes.
#include <iostream>  // NOLINT(readability/streams)
#include <limits>
#include <string>
#include <type_traits>

#include "base/not_null.hpp"
#include "serialization/quantities.pb.h"

namespace principia {

using base::not_null;

namespace quantities {

template<int64_t LengthExponent, int64_t MassExponent, int64_t TimeExponent,
         int64_t CurrentExponent, int64_t TemperatureExponent,
         int64_t AmountExponent, int64_t LuminousIntensityExponent,
         int64_t WindingExponent, int64_t AngleExponent,
         int64_t SolidAngleExponent>
struct Dimensions;
template<typename D> class Quantity;

using NoDimensions = Dimensions<0, 0, 0, 0, 0, 0, 0, 0, 0, 0>;

// Base quantities
using Length            = Quantity<Dimensions<1, 0, 0, 0, 0, 0, 0, 0, 0, 0>>;
using Mass              = Quantity<Dimensions<0, 1, 0, 0, 0, 0, 0, 0, 0, 0>>;
using Time              = Quantity<Dimensions<0, 0, 1, 0, 0, 0, 0, 0, 0, 0>>;
using Current           = Quantity<Dimensions<0, 0, 0, 1, 0, 0, 0, 0, 0, 0>>;
using Temperature       = Quantity<Dimensions<0, 0, 0, 0, 1, 0, 0, 0, 0, 0>>;
using Amount            = Quantity<Dimensions<0, 0, 0, 0, 0, 1, 0, 0, 0, 0>>;
using LuminousIntensity = Quantity<Dimensions<0, 0, 0, 0, 0, 0, 1, 0, 0, 0>>;
// Nonstandard; winding is a dimensionless quantity counting cycles, in order to
// strongly type the distinction between Frequency = Winding/Time and
// AngularFrequency = Angle/Time. We also strongly type angles.
using Winding           = Quantity<Dimensions<0, 0, 0, 0, 0, 0, 0, 1, 0, 0>>;
using Angle             = Quantity<Dimensions<0, 0, 0, 0, 0, 0, 0, 0, 1, 0>>;
using SolidAngle        = Quantity<Dimensions<0, 0, 0, 0, 0, 0, 0, 0, 0, 1>>;

namespace internal {

template<typename Left, typename Right> struct ProductGenerator;
template<typename Left, typename Right> struct QuotientGenerator;
template<typename Q, typename = void> struct SquareRootGenerator;
template<typename T, int exponent, typename = void>
struct ExponentiationGenerator;

template<typename Left, typename Right>
using Product = typename ProductGenerator<Left, Right>::Type;
template<typename Left, typename Right>
using Quotient = typename QuotientGenerator<Left, Right>::Type;

}  // namespace internal

// The result type of +, -, * and / on arguments of types |Left| and |Right|.
template<typename Left, typename Right>
using Sum = decltype(std::declval<Left>() = std::declval<Right>());
template<typename Left, typename Right = Left>
using Difference = decltype(std::declval<Left>() - std::declval<Right>());
template<typename Left, typename Right>
using Product = decltype(std::declval<Left>() * std::declval<Right>());
template<typename Left, typename Right>
using Quotient = decltype(std::declval<Left>() / std::declval<Right>());

// |Exponentiation<T, n>| is an alias for the following, where t is a value of
// type |T|:
//   The type of ( ... (t * t) * ... * t), with n factors, if n >= 1;
//   The type of t / ( ... (t * t) * ... * t), with n + 1 factors in the
//   denominator, if n < 1.
template<typename T, int exponent>
using Exponentiation =
    typename internal::ExponentiationGenerator<T, exponent>::Type;

// |SquareRoot<T>| is only defined if |T| is an instance of |Quantity| with only
// even dimensions.  In that case, it is the unique instance |S| of |Quantity|
// such that |Product<S, S>| is |T|.
template<typename Q>
using SquareRoot = typename internal::SquareRootGenerator<Q>::Type;

// Returns the base or derived SI Unit of |Q|.
// For instance, |SIUnit<Action>() == Joule * Second|.
template<typename Q>
constexpr Q SIUnit();
// Returns 1.
template<>
constexpr double SIUnit<double>();

template<typename LDimensions, typename RDimensions>
constexpr internal::Product<Quantity<LDimensions>, Quantity<RDimensions>>
operator*(Quantity<LDimensions> const&, Quantity<RDimensions> const&);
template<typename LDimensions, typename RDimensions>
constexpr internal::Quotient<Quantity<LDimensions>, Quantity<RDimensions>>
operator/(Quantity<LDimensions> const&, Quantity<RDimensions> const&);
template<typename RDimensions>
constexpr Quantity<RDimensions>
operator*(double const, Quantity<RDimensions> const&);
template<typename RDimensions>
constexpr typename Quantity<RDimensions>::Inverse
operator/(double const, Quantity<RDimensions> const&);

// Equivalent to |std::pow(x, exponent)| unless -3 ≤ x ≤ 3, in which case
// explicit specialization yields multiplications statically.
template<int exponent>
constexpr double Pow(double x);
template<int exponent, typename D>
constexpr Exponentiation<Quantity<D>, exponent> Pow(Quantity<D> const& x);

template<typename D>
std::ostream& operator<<(std::ostream& out, Quantity<D> const& quantity);

// Equivalent to |std::abs(x)|.
double Abs(double const x);
template<typename D>
Quantity<D> Abs(Quantity<D> const& x);

template<typename D>
SquareRoot<Quantity<D>> Sqrt(Quantity<D> const& x);

template<typename D>
Angle ArcTan(Quantity<D> const& y, Quantity<D> const& x);

std::string DebugString(
    double const number,
    int const precision = std::numeric_limits<double>::max_digits10);
template<typename D>
std::string DebugString(
    Quantity<D> const& quantity,
    int const precision = std::numeric_limits<double>::max_digits10);

template<typename D>
class Quantity {
 public:
  using Dimensions = D;
  using Inverse = internal::Quotient<double, Quantity>;

  constexpr Quantity();
  ~Quantity() = default;

  constexpr Quantity operator+() const;
  constexpr Quantity operator-() const;
  constexpr Quantity operator+(Quantity const& right) const;
  constexpr Quantity operator-(Quantity const& right) const;

  constexpr Quantity operator*(double const right) const;
  constexpr Quantity operator/(double const right) const;

  Quantity& operator+=(Quantity const&);
  Quantity& operator-=(Quantity const&);
  Quantity& operator*=(double const);
  Quantity& operator/=(double const);

  constexpr bool operator>(Quantity const& right) const;
  constexpr bool operator<(Quantity const& right) const;
  constexpr bool operator>=(Quantity const& right) const;
  constexpr bool operator<=(Quantity const& right) const;
  constexpr bool operator==(Quantity const& right) const;
  constexpr bool operator!=(Quantity const& right) const;

  void WriteToMessage(not_null<serialization::Quantity*> const message) const;
  static Quantity ReadFromMessage(serialization::Quantity const& message);

 private:
  explicit constexpr Quantity(double const magnitude);
  double magnitude_;

  template<typename LDimensions, typename RDimensions>
  friend constexpr internal::Product<Quantity<LDimensions>,
                                     Quantity<RDimensions>> operator*(
      Quantity<LDimensions> const& left,
      Quantity<RDimensions> const& right);
  template<typename LDimensions, typename RDimensions>
  friend constexpr internal::Quotient<Quantity<LDimensions>,
                                      Quantity<RDimensions>> operator/(
      Quantity<LDimensions> const& left, 
      Quantity<RDimensions> const& right);
  template<typename RDimensions>
  friend constexpr Quantity<RDimensions> operator*(
      double const left,
      Quantity<RDimensions> const& right);
  template<typename RDimensions>
  friend constexpr typename Quantity<RDimensions>::Inverse operator/(
      double const left,
      Quantity<RDimensions> const& right);

  template<typename Q>
  friend constexpr Q SIUnit();

  template<int exponent, typename BaseDimensions>
  friend constexpr Exponentiation<Quantity<BaseDimensions>, exponent> Pow(
      Quantity<BaseDimensions> const& x);

  friend Quantity<D> Abs<>(Quantity<D> const&);

  template<typename ArgumentDimensions>
  friend SquareRoot<Quantity<ArgumentDimensions>> Sqrt(
      Quantity<ArgumentDimensions> const& x);
  friend Angle ArcTan<>(Quantity<D> const& y, Quantity<D> const& x);

  friend std::string DebugString<>(Quantity<D> const&, int const);
};

// A type trait for testing if a type is a quantity.
template<typename T>
struct is_quantity : std::is_floating_point<T> {};
template<typename D>
struct is_quantity<Quantity<D>> : std::true_type {};

}  // namespace quantities
}  // namespace principia

#include "quantities/quantities_body.hpp"
