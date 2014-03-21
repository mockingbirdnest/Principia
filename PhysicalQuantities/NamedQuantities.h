// NamedQuantities.h

#pragma once

#include "PhysicalQuantities.h"

namespace PhysicalQuantities {
#pragma region General mechanics
  typedef Quotient<Length, Time>   Speed;
  typedef Quotient<Speed, Time>    Acceleration;
  typedef Product<Mass, Speed>     Momentum;
  typedef Quotient<Momentum, Time> Force;
  typedef Product<Force, Length>   Energy;
#pragma endregion
#pragma region Thermodynamics
  typedef Product<Length, Length>       Surface;
  typedef Product<Surface, Length>      Volume;
  typedef Quotient<Force, Surface>      Pressure;
  typedef Quotient<Energy, Temperature> Entropy;
#pragma endregion
}