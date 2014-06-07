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
std::string ToString(Quantity<D> const& quantity,
                     unsigned char const precision = 16);

template<typename D>
std::ostream& operator<<(std::ostream& out, Quantity<D> const& quantity);

template<typename D>
class Quantity {
 public:
  typedef typename D Dimensions;
  typedef Quotient<double, Quantity> Inverse;

  Quantity();
  ~Quantity() = default;

  Quantity& operator+=(Quantity const&);
  Quantity& operator-=(Quantity const&);
  Quantity& operator*=(double const);
  Quantity& operator/=(double const);

 private:
  explicit Quantity(double const magnitude);
  double magnitude_;

  template<typename D>
  friend class Quantity;

  template<typename D>
  friend Quantity<D> operator+(Quantity<D> const&);
  template<typename D>
  friend Quantity<D> operator-(Quantity<D> const&);
  template<typename D>
  friend Quantity<D> operator+(Quantity<D> const&, Quantity<D> const&);
  template<typename D>
  friend Quantity<D> operator-(Quantity<D> const&, Quantity<D> const&);
  template<typename DLeft, typename DRight>
  friend Product<typename Quantity<DLeft>,
                 typename Quantity<DRight>> operator*(Quantity<DLeft> const&,
                                                      Quantity<DRight> const&);
  template<typename DLeft, typename DRight>
  friend Quotient<typename Quantity<DLeft>,
                  typename Quantity<DRight>> operator/(Quantity<DLeft> const&,
                                                       Quantity<DRight> const&);
  template<typename D>
  friend Quantity<D> operator*(Quantity<D> const&, double const);
  template<typename D>
  friend Quantity<D> operator*(double const, Quantity<D> const&);
  template<typename D>
  friend Quantity<D> operator/(Quantity<D> const&, double const);
  template<typename D>
  friend typename Quantity<D>::Inverse operator/(double const,
                                                 Quantity<D> const&);

  template<typename D>
  friend bool operator>(Quantity<D> const&, Quantity<D> const&);
  template<typename D>
  friend bool operator<(Quantity<D> const&, Quantity<D> const&);
  template<typename D>
  friend bool operator>=(Quantity<D> const&, Quantity<D> const&);
  template<typename D>
  friend bool operator<=(Quantity<D> const&, Quantity<D> const&);
  template<typename D>
  friend bool operator==(Quantity<D> const&, Quantity<D> const&);
  template<typename D>
  friend bool operator!=(Quantity<D> const&, Quantity<D> const&);

  template<typename Q>
  friend Q SIUnit();

  template<int exponent, typename D>
  friend Exponentiation<Quantity<D>, exponent> Pow(Quantity<D> const& x);

  template<typename D>
  friend Quantity<D> Abs(Quantity<D> const&);
  template<typename D>
  friend SquareRoot<Quantity<D>> Sqrt(Quantity<D> const& x);
  template<typename D>
  friend Angle ArcTan(Quantity<D> const& y, Quantity<D> const& x);

  template<typename D>
  friend std::string ToString(Quantity<D> const&, unsigned char const);
};

template<typename D>
Quantity<D> operator+(Quantity<D> const&);
template<typename D>
Quantity<D> operator-(Quantity<D> const&);
template<typename D>
Quantity<D> operator+(Quantity<D> const&, Quantity<D> const&);
template<typename D>
Quantity<D> operator-(Quantity<D> const&, Quantity<D> const&);
template<typename DLeft, typename DRight>
Product<typename Quantity<DLeft>,
                typename Quantity<DRight>> operator*(Quantity<DLeft> const&,
                                                    Quantity<DRight> const&);
template<typename DLeft, typename DRight>
Quotient<typename Quantity<DLeft>,
                typename Quantity<DRight>> operator/(Quantity<DLeft> const&,
                                                      Quantity<DRight> const&);
template<typename D>
Quantity<D> operator*(Quantity<D> const&, double const);
template<typename D>
Quantity<D> operator*(double const, Quantity<D> const&);
template<typename D>
Quantity<D> operator/(Quantity<D> const&, double const);
template<typename D>
typename Quantity<D>::Inverse operator/(double const,
                                                Quantity<D> const&);

template<typename D>
bool operator>(Quantity<D> const&, Quantity<D> const&);
template<typename D>
bool operator<(Quantity<D> const&, Quantity<D> const&);
template<typename D>
bool operator>=(Quantity<D> const&, Quantity<D> const&);
template<typename D>
bool operator<=(Quantity<D> const&, Quantity<D> const&);
template<typename D>
bool operator==(Quantity<D> const&, Quantity<D> const&);
template<typename D>
bool operator!=(Quantity<D> const&, Quantity<D> const&);

// Returns the base or derived SI Unit of |Q|.
// For instance, |SIUnit<Action>() == Joule * Second|.
template<typename Q>
Q SIUnit();
template<>
double SIUnit<double>();

// Equivalent to |std::pow(x, exponent)| unles -3 ≤ x ≤ 3, in which case
// explicit specialisation yields multiplications statically.
template<int exponent>
double Pow(double x);
template<int exponent, typename D>
Exponentiation<Quantity<D>, exponent> Pow(Quantity<D> const& x);

// Equivalent to |std::abs(x)|.
double Abs(double const x);
template<typename D>
Quantity<D> Abs(Quantity<D> const& x);

template<typename D>
SquareRoot<Quantity<D>> Sqrt(Quantity<D> const& x);
template<typename D>
Angle ArcTan(Quantity<D> const& y, Quantity<D> const& x);

std::string ToString(double const number,
                     unsigned char const precision = DBL_DIG + 1);
template<typename D>
std::string ToString(Quantity<D> const& quantity,
                     unsigned char const precision);

template<typename D>
std::ostream& operator<<(::std::ostream& out, Quantity<D> const& quantity);

}  // namespace quantities
}  // namespace principia

#include "quantities/quantities_body.hpp"
