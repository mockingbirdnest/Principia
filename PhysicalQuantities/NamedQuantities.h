// NamedQuantities.h

#pragma once

#include "PhysicalQuantities.h"

namespace PhysicalQuantities {
#pragma region General mechanics
  typedef Quotient<Length, Time>    Speed;
  typedef Quotient<Speed, Time>     Acceleration;
  typedef Product<Mass, Speed>      Momentum;
  typedef Quotient<Momentum, Time>  Force;
  typedef Product<Force, Length>    Energy;
  typedef Quotient<Energy, Time>    Power;
  typedef Product<Energy, Time>     Action;
  typedef Inverse<Time>             AngularFrequency;
  typedef Product<Length, Momentum> AngularMomentum;
#pragma endregion
#pragma region Thermodynamics
  typedef Product<Length, Length>       Surface;
  typedef Product<Surface, Length>      Volume;
  typedef Quotient<Force, Surface>      Pressure;
  typedef Quotient<Energy, Temperature> Entropy;
  typedef Quotient<Mass, Volume>        Density;
#pragma endregion
}