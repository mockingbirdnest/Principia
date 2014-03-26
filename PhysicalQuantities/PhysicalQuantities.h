// PhysicalQuantities.h

#pragma once

#include "Dimensionless.h"

namespace PhysicalQuantities {
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
template<typename Left, typename Right> struct ProductGenerator;
template<typename Left, typename Right> struct QuotientGenerator;
template<typename Left, typename Right>
using Quotient = typename QuotientGenerator<Left, Right>::ResultType;
template<typename Left, typename Right>
using Product = typename ProductGenerator<Left, Right>::ResultType;
template<typename Right>
using Inverse = Quotient<Dimensionless, Right>;

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

template<typename D>
class Quantity {
 public:
  Quantity() = default;
  typedef typename D Dimensions;
 private:
  explicit Quantity(Dimensionless const magnitude) : magnitude_(magnitude) {}
  Dimensionless magnitude_;

  friend Length            Metres(Dimensionless const&);
  friend Mass              Kilograms(Dimensionless const&);
  friend Time              Seconds(Dimensionless const&);
  friend Current           Amperes(Dimensionless const&);
  friend Temperature       Kelvins(Dimensionless const&);
  friend Amount            Moles(Dimensionless const&);
  friend LuminousIntensity Candelas(Dimensionless const&);
  friend Winding           Cycles(Dimensionless const&);
  friend Angle             Radians(Dimensionless const&);
  friend SolidAngle        Steradians(Dimensionless const&);

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

#include "PhysicalQuantities.ipp"
