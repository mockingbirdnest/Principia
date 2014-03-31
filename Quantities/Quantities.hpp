// Quantities.h

#pragma once

#include "Dimensionless.hpp"
#include <string>

namespace Principia {
namespace Quantities {
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
namespace TypeGenerators {
template<typename Left, typename Right> struct ProductGenerator;
template<typename Left, typename Right> struct QuotientGenerator;
template<bool> struct Range;
template<typename Q, int Exponent, typename = Range<true>> 
struct PowerGenerator;
}
template<typename Left, typename Right>
using Quotient = typename TypeGenerators::QuotientGenerator<Left,
                                                            Right>::ResultType;
template<typename Left, typename Right>
using Product = typename TypeGenerators::ProductGenerator<Left,
                                                          Right>::ResultType;
template<typename Left, int Exponent>
using Exponentiation =
  typename TypeGenerators::PowerGenerator<Left, Exponent>::ResultType;
template<typename Right>
using Inverse = Quotient<Dimensionless, Right>;

namespace Factories {
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
}

template<typename D>
std::wstring ToString(Quantity<D> const& quantity,
                      unsigned char const precision = 16);

template<typename D>
class Quantity {
 public:
  typedef typename D Dimensions;
  Quantity() = default;
  template<int Exponent>
  Exponentiation<Quantity<D>, Exponent> Pow() const;
 private:
  explicit Quantity(Dimensionless const magnitude) : magnitude_(magnitude) {}
  Dimensionless magnitude_;

  friend Length            Factories::Metres(Dimensionless const&);
  friend Mass              Factories::Kilograms(Dimensionless const&);
  friend Time              Factories::Seconds(Dimensionless const&);
  friend Current           Factories::Amperes(Dimensionless const&);
  friend Temperature       Factories::Kelvins(Dimensionless const&);
  friend Amount            Factories::Moles(Dimensionless const&);
  friend LuminousIntensity Factories::Candelas(Dimensionless const&);
  friend Winding           Factories::Cycles(Dimensionless const&);
  friend Angle             Factories::Radians(Dimensionless const&);
  friend SolidAngle        Factories::Steradians(Dimensionless const&);

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
}
}

#include "Quantities-inl.hpp"
