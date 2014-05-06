#pragma once

#include <string>

#include "Quantities/Dimensionless.hpp"

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
using Quotient = typename type_generators::QuotientGenerator<Left,
                                                             Right>::ResultType;
template<typename Left, typename Right>
using Product = typename type_generators::ProductGenerator<Left,
                                                           Right>::ResultType;
template<typename Left, int Exponent>
using Exponentiation =
    typename type_generators::PowerGenerator<Left, Exponent>::ResultType;
template<typename Q>
using SquareRoot = typename type_generators::SquareRootGenerator<Q>::ResultType;
template<typename Right>
using Inverse = Quotient<Dimensionless, Right>;

namespace factories {
Length            Metres(Dimensionless const&);
Mass              Kilograms(Dimensionless const&);
Time              Seconds(Dimensionless const&);
Current           Amperes(Dimensionless const&);
Temperature       Kelvins(Dimensionless const&);
Amount            Moles(Dimensionless const&);
LuminousIntensity Candelas(Dimensionless const&);
Winding           Cycles(Dimensionless const&);
Angle             Radians(Dimensionless const&);
SolidAngle        Steradians(Dimensionless const&);
}  // namespace factories

template<typename D>
std::wstring ToString(Quantity<D> const& quantity,
                      unsigned char const precision = 16);

template<typename D>
class Quantity {
 public:
  typedef typename D Dimensions;
  Quantity();
  template<int Exponent>
  Exponentiation<Quantity<D>, Exponent> Pow() const;
 private:
  explicit Quantity(Dimensionless const& magnitude);
  Dimensionless magnitude_;

  friend Length            factories::Metres(Dimensionless const&);
  friend Mass              factories::Kilograms(Dimensionless const&);
  friend Time              factories::Seconds(Dimensionless const&);
  friend Current           factories::Amperes(Dimensionless const&);
  friend Temperature       factories::Kelvins(Dimensionless const&);
  friend Amount            factories::Moles(Dimensionless const&);
  friend LuminousIntensity factories::Candelas(Dimensionless const&);
  friend Winding           factories::Cycles(Dimensionless const&);
  friend Angle             factories::Radians(Dimensionless const&);
  friend SolidAngle        factories::Steradians(Dimensionless const&);

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
  friend Quantity<D> operator*(Quantity<D> const&, Dimensionless const&);
  template<typename D>
  friend Quantity<D> operator*(Dimensionless const&, Quantity<D> const&);
  template<typename D>
  friend Quantity<D> operator/(Quantity<D> const&, Dimensionless const&);
  template<typename D>
  friend Inverse<Quantity<D>> operator/(Dimensionless const&,
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

  template<typename D>
  friend Quantity<D> Abs(Quantity<D> const&);

  template<typename D>
  friend SquareRoot<Quantity<D>> Sqrt(Quantity<D> const& x);
  template<typename D>
  friend Angle ArcTan(Quantity<D> const& y, Quantity<D> const& x);

  template<typename D>
  friend std::wstring ToString(Quantity<D> const&, unsigned char const);
};

template<typename D>
void operator+=(Quantity<D>&, Quantity<D> const&);
template<typename D>
inline void operator-=(Quantity<D>&, Quantity<D> const&);
template<typename D>
inline void operator*=(Quantity<D>&, Dimensionless const&);
template<typename D>
inline void operator/=(Quantity<D>&, Dimensionless const&);
}  // namespace quantities
}  // namespace principia

#include "Quantities/Quantities-body.hpp"
