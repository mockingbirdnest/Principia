#pragma once

#include <cfloat>
// We use ostream for logging purposes.
#include <iostream>  // NOLINT(readability/streams)
#include <string>

namespace principia {
namespace quantities {
template<int LengthExponent, int MassExponent, int TimeExponent,
         int CurrentExponent, int TemperatureExponent, int AmountExponent,
         int LuminousIntensityExponent, int WindingExponent,
         int AngleExponent, int SolidAngleExponent>
struct Dimensions;
template<typename D> class Quantity;

typedef Dimensions<0, 0, 0, 0, 0, 0, 0, 0, 0, 0> NoDimensions;

#pragma region Base quantities
typedef Quantity<Dimensions<1, 0, 0, 0, 0, 0, 0, 0, 0, 0>> Length;
typedef Quantity<Dimensions<0, 1, 0, 0, 0, 0, 0, 0, 0, 0>> Mass;
typedef Quantity<Dimensions<0, 0, 1, 0, 0, 0, 0, 0, 0, 0>> Time;
typedef Quantity<Dimensions<0, 0, 0, 1, 0, 0, 0, 0, 0, 0>> Current;
typedef Quantity<Dimensions<0, 0, 0, 0, 1, 0, 0, 0, 0, 0>> Temperature;
typedef Quantity<Dimensions<0, 0, 0, 0, 0, 1, 0, 0, 0, 0>> Amount;
typedef Quantity<Dimensions<0, 0, 0, 0, 0, 0, 1, 0, 0, 0>> LuminousIntensity;
// Nonstandard; winding is a dimensionless quantity counting cycles, in order to
// strongly type the distinction between Frequency = Winding/Time and
// AngularFrequency = Angle/Time. We also strongly type angles.
typedef Quantity<Dimensions<0, 0, 0, 0, 0, 0, 0, 1, 0, 0>> Winding;
typedef Quantity<Dimensions<0, 0, 0, 0, 0, 0, 0, 0, 1, 0>> Angle;
typedef Quantity<Dimensions<0, 0, 0, 0, 0, 0, 0, 0, 0, 1>> SolidAngle;
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

template<typename D>
class Quantity {
 public:
  typedef D Dimensions;
  typedef Quotient<double, Quantity> Inverse;

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

  template<typename ArgumentDimensions>
  friend Quantity<ArgumentDimensions> Abs(Quantity<ArgumentDimensions> const&);
  template<typename ArgumentDimensions>
  friend SquareRoot<Quantity<ArgumentDimensions>> Sqrt(
      Quantity<ArgumentDimensions> const& x);
  template<typename ArgumentDimensions>
  friend Angle ArcTan(Quantity<ArgumentDimensions> const& y,
                      Quantity<ArgumentDimensions> const& x);

  template<typename ArgumentDimensions>
  friend std::string DebugString(Quantity<ArgumentDimensions> const&,
                                 unsigned char const);
};

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

// Returns the base or derived SI Unit of |Q|.
// For instance, |SIUnit<Action>() == Joule * Second|.
template<typename Q>
Q SIUnit();
// Returns 1.
template<>
double SIUnit<double>();

// Equivalent to |std::pow(x, exponent)| unless -3 ≤ x ≤ 3, in which case
// explicit specialisation yields multiplications statically.
template<int exponent>
double Pow(double x);
template<int exponent, typename BaseDimensions>
Exponentiation<Quantity<BaseDimensions>, exponent> Pow(
    Quantity<BaseDimensions> const& x);

// Equivalent to |std::abs(x)|.
double Abs(double const x);
template<typename ArgumentDimensions>
Quantity<ArgumentDimensions> Abs(Quantity<ArgumentDimensions> const& x);

template<typename ArgumentDimensions>
SquareRoot<Quantity<ArgumentDimensions>> Sqrt(
    Quantity<ArgumentDimensions> const& x);

template<typename ArgumentDimensions>
Angle ArcTan(Quantity<ArgumentDimensions> const& y,
             Quantity<ArgumentDimensions> const& x);

std::string DebugString(double const number,
                        unsigned char const precision = DBL_DIG + 1);
template<typename ArgumentDimensions>
std::string DebugString(Quantity<ArgumentDimensions> const& quantity,
                        unsigned char const precision = DBL_DIG + 1);

template<typename ArgumentDimensions>
std::ostream& operator<<(std::ostream& out,
                         Quantity<ArgumentDimensions> const& quantity);

}  // namespace quantities
}  // namespace principia

#include "quantities/quantities_body.hpp"
