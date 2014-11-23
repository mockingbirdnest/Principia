#pragma once

// We use ostream for logging purposes.
#include <iostream>  // NOLINT(readability/streams)
#include <limits>
#include <string>

namespace principia {
namespace quantities {
template<int LengthExponent, int MassExponent, int TimeExponent,
         int CurrentExponent, int TemperatureExponent, int AmountExponent,
         int LuminousIntensityExponent, int WindingExponent,
         int AngleExponent, int SolidAngleExponent>
struct Dimensions;
template<typename D> class Quantity;

using NoDimensions = Dimensions<0, 0, 0, 0, 0, 0, 0, 0, 0, 0>;

#pragma region Base quantities
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
#pragma endregion

namespace type_generators {
template<typename Left, typename Right> struct ProductGenerator;
template<typename Left, typename Right> struct QuotientGenerator;
template<bool> struct Range;
template<typename Q, int Exponent, typename = Range<true>>
struct PowerGenerator;
template<bool> struct Condition;
template<typename Q, typename = Condition<true>> struct SquareRootGenerator;
}  // namespace type_generators

template<typename Left, typename Right>
using Quotient =
    typename type_generators::QuotientGenerator<Left, Right>::ResultType;
template<typename Left, typename Right>
using Product =
    typename type_generators::ProductGenerator<Left, Right>::ResultType;
template<typename Left, int Exponent>
using Exponentiation =
    typename type_generators::PowerGenerator<Left, Exponent>::ResultType;
template<typename Q>
using SquareRoot = typename type_generators::SquareRootGenerator<Q>::ResultType;

// Returns the base or derived SI Unit of |Q|.
// For instance, |SIUnit<Action>() == Joule * Second|.
template<typename Q>
Q SIUnit();
// Returns 1.
template<>
double SIUnit<double>();

template<typename LDimensions, typename RDimensions>
Product<Quantity<LDimensions>, Quantity<RDimensions>> operator*(
    Quantity<LDimensions> const&,
    Quantity<RDimensions> const&);
template<typename LDimensions, typename RDimensions>
Quotient<Quantity<LDimensions>, Quantity<RDimensions>> operator/(
    Quantity<LDimensions> const&,
    Quantity<RDimensions> const&);
template<typename RDimensions>
Quantity<RDimensions> operator*(double const, Quantity<RDimensions> const&);
template<typename RDimensions>
typename Quantity<RDimensions>::Inverse operator/(double const,
                                                  Quantity<RDimensions> const&);

// Equivalent to |std::pow(x, exponent)| unless -3 ≤ x ≤ 3, in which case
// explicit specialization yields multiplications statically.
template<int exponent>
double Pow(double x);
template<int exponent, typename D>
Exponentiation<Quantity<D>, exponent> Pow(Quantity<D> const& x);

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
  using Inverse = Quotient<double, Quantity>;

  Quantity();
  ~Quantity() = default;

  Quantity operator+() const;
  Quantity operator-() const;
  Quantity operator+(Quantity const& right) const;
  Quantity operator-(Quantity const& right) const;

  Quantity operator*(double const right) const;
  Quantity operator/(double const right) const;

  Quantity& operator+=(Quantity const&);
  Quantity& operator-=(Quantity const&);
  Quantity& operator*=(double const);
  Quantity& operator/=(double const);

  bool operator>(Quantity const& right) const;
  bool operator<(Quantity const& right) const;
  bool operator>=(Quantity const& right) const;
  bool operator<=(Quantity const& right) const;
  bool operator==(Quantity const& right) const;
  bool operator!=(Quantity const& right) const;

 private:
  explicit Quantity(double const magnitude);
  double magnitude_;

  template<typename LDimensions, typename RDimensions>
  friend Product<Quantity<LDimensions>, Quantity<RDimensions>> operator*(
      Quantity<LDimensions> const& left,
      Quantity<RDimensions> const& right);
  template<typename LDimensions, typename RDimensions>
  friend Quotient<Quantity<LDimensions>, Quantity<RDimensions>> operator/(
      Quantity<LDimensions> const& left,
      Quantity<RDimensions> const& right);
  template<typename RDimensions>
  friend Quantity<RDimensions> operator*(double const left,
                                         Quantity<RDimensions> const& right);
  template<typename RDimensions>
  friend typename Quantity<RDimensions>::Inverse operator/(
      double const left,
      Quantity<RDimensions> const& right);

  template<typename Q>
  friend Q SIUnit();

  template<int exponent, typename BaseDimensions>
  friend Exponentiation<Quantity<BaseDimensions>, exponent> Pow(
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
