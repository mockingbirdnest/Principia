
#pragma once

#include <nmmintrin.h>

#include <iostream>
#include <limits>
#include <string>
#include <type_traits>

#include "base/not_constructible.hpp"
#include "base/not_null.hpp"
#include "quantities/dimensions.hpp"
#include "quantities/generators.hpp"
#include "serialization/quantities.pb.h"

namespace principia {
namespace quantities {
namespace internal_quantities {

using base::not_constructible;
using base::not_null;
using internal_dimensions::Dimensions;
using internal_generators::ExponentiationGenerator;
using internal_generators::NthRootGenerator;
using internal_generators::ProductGenerator;
using internal_generators::QuotientGenerator;

template<typename D>
class Quantity;

// Base quantities
using Length            = Quantity<Dimensions<1, 0, 0, 0, 0, 0, 0, 0>>;
using Mass              = Quantity<Dimensions<0, 1, 0, 0, 0, 0, 0, 0>>;
using Time              = Quantity<Dimensions<0, 0, 1, 0, 0, 0, 0, 0>>;
using Current           = Quantity<Dimensions<0, 0, 0, 1, 0, 0, 0, 0>>;
using Temperature       = Quantity<Dimensions<0, 0, 0, 0, 1, 0, 0, 0>>;
using Amount            = Quantity<Dimensions<0, 0, 0, 0, 0, 1, 0, 0>>;
using LuminousIntensity = Quantity<Dimensions<0, 0, 0, 0, 0, 0, 1, 0>>;
// We strongly type angles.
using Angle             = Quantity<Dimensions<0, 0, 0, 0, 0, 0, 0, 1>>;

// |Product| and |Quotient| are not exported from this namespace.  Instead they
// are defined as the result types of |operator*| and |operator/|.
template<typename Left, typename Right>
using Product = typename ProductGenerator<Left, Right>::Type;
template<typename Left, typename Right>
using Quotient = typename QuotientGenerator<Left, Right>::Type;

template<typename T, int exponent>
using Exponentiation = typename ExponentiationGenerator<T, exponent>::Type;
template<typename Q>
using Square = Exponentiation<Q, 2>;
template<typename Q>
using Cube = Exponentiation<Q, 3>;

template<typename Q, int n>
using NthRoot = typename NthRootGenerator<Q, n>::Type;
template<typename Q>
using SquareRoot = NthRoot<Q, 2>;
template<typename Q>
using CubeRoot = NthRoot<Q, 3>;

template<typename D>
class Quantity final {
 public:
  using Dimensions = D;
  using Inverse = Quotient<double, Quantity>;

  constexpr Quantity();

  constexpr Quantity operator+() const;
  constexpr Quantity operator-() const;
  constexpr Quantity operator+(Quantity const& right) const;
  constexpr Quantity operator-(Quantity const& right) const;

  constexpr Quantity operator*(double right) const;
  constexpr Quantity operator/(double right) const;

  Quantity& operator+=(Quantity const& right);
  Quantity& operator-=(Quantity const& right);
  Quantity& operator*=(double right);
  Quantity& operator/=(double right);

  constexpr bool operator>(Quantity const& right) const;
  constexpr bool operator<(Quantity const& right) const;
  constexpr bool operator>=(Quantity const& right) const;
  constexpr bool operator<=(Quantity const& right) const;
  constexpr bool operator==(Quantity const& right) const;
  constexpr bool operator!=(Quantity const& right) const;

  void WriteToMessage(not_null<serialization::Quantity*> message) const;
  static Quantity ReadFromMessage(serialization::Quantity const& message);

 private:
  explicit constexpr Quantity(double magnitude);
  double magnitude_;

  template<typename LDimensions, typename RDimensions>
  friend constexpr Product<Quantity<LDimensions>,
                           Quantity<RDimensions>> operator*(
      Quantity<LDimensions> const& left,
      Quantity<RDimensions> const& right);
  template<typename LDimensions, typename RDimensions>
  friend constexpr Quotient<Quantity<LDimensions>,
                            Quantity<RDimensions>> operator/(
      Quantity<LDimensions> const& left,
      Quantity<RDimensions> const& right);
  template<typename RDimensions>
  friend constexpr Quantity<RDimensions> operator*(
      double left,
      Quantity<RDimensions> const& right);
  template<typename RDimensions>
  friend constexpr typename Quantity<RDimensions>::Inverse operator/(
      double left,
      Quantity<RDimensions> const& right);

  template<typename Q>
  friend constexpr Q Infinity();
  template<typename Q>
  friend constexpr Q NaN();
  template<typename Q>
  friend constexpr Q SIUnit();
  template<typename Q>
  friend Q FromM128D(__m128d x);
  template<typename Q>
  friend __m128d ToM128D(Q x);
};

template<typename LDimensions, typename RDimensions>
constexpr Product<Quantity<LDimensions>, Quantity<RDimensions>>
operator*(Quantity<LDimensions> const&, Quantity<RDimensions> const&);
template<typename LDimensions, typename RDimensions>
constexpr Quotient<Quantity<LDimensions>, Quantity<RDimensions>>
operator/(Quantity<LDimensions> const&, Quantity<RDimensions> const&);
template<typename RDimensions>
constexpr Quantity<RDimensions>
operator*(double, Quantity<RDimensions> const&);
template<typename RDimensions>
constexpr typename Quantity<RDimensions>::Inverse
operator/(double, Quantity<RDimensions> const&);

// Returns the base or derived SI Unit of |Q|.
// For instance, |SIUnit<Action>() == Joule * Second|.
template<typename Q>
constexpr Q SIUnit();
// Returns 1.
template<>
constexpr double SIUnit<double>();

// Conversion to and from intrinsic types.  ToM128D fills both halves of the
// result.
template<typename Q>
Q FromM128D(__m128d x);
template<typename Q>
__m128d ToM128D(Q);
template<>
double FromM128D(__m128d x);
template<>
__m128d ToM128D(double x);

// A type trait for testing if a type is a quantity.
template<typename T>
struct is_quantity : std::is_arithmetic<T>, not_constructible {};
template<typename D>
struct is_quantity<Quantity<D>> : std::true_type, not_constructible {};

// Returns a positive infinity of |Q|.
template<typename Q, typename = std::enable_if<is_quantity<Q>::value>>
constexpr Q Infinity();
template<typename Q, typename = std::enable_if<is_quantity<Q>::value>>
constexpr bool IsFinite(Q const& x);

// Returns a quiet NaN of |Q|.
template<typename Q, typename = std::enable_if<is_quantity<Q>::value>>
constexpr Q NaN();

std::string DebugString(
    double number,
    int precision = std::numeric_limits<double>::max_digits10);
template<typename D>
std::string DebugString(
    Quantity<D> const& quantity,
    int precision = std::numeric_limits<double>::max_digits10);

template<typename D>
std::ostream& operator<<(std::ostream& out, Quantity<D> const& quantity);

}  // namespace internal_quantities

using internal_quantities::Amount;
using internal_quantities::Angle;
using internal_quantities::Cube;
using internal_quantities::CubeRoot;
using internal_quantities::Current;
using internal_quantities::DebugString;
using internal_quantities::Exponentiation;
using internal_quantities::FromM128D;
using internal_quantities::Infinity;
using internal_quantities::IsFinite;
using internal_quantities::is_quantity;
using internal_quantities::Length;
using internal_quantities::LuminousIntensity;
using internal_quantities::Mass;
using internal_quantities::NaN;
using internal_quantities::Quantity;
using internal_quantities::SIUnit;
using internal_quantities::Square;
using internal_quantities::SquareRoot;
using internal_quantities::Temperature;
using internal_quantities::Time;
using internal_quantities::ToM128D;

}  // namespace quantities
}  // namespace principia

#include "quantities/quantities_body.hpp"
