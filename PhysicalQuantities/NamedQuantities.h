// NamedQuantities.h

#pragma once

#include "PhysicalQuantities.h"

namespace PhysicalQuantities {
#pragma region General mechanics
#pragma region Linear Motion
typedef Quotient<Length, Time>   Speed;
typedef Quotient<Speed, Time>    Acceleration;
typedef Product<Mass, Speed>     Momentum;
typedef Quotient<Momentum, Time> Force;
#pragma endregion Linear Motion
#pragma region Rotational Motion
// There is no good way of strongly typing angles vs. dimensionless scalars
// because angles *are* dimensionless scalars. This is just a convenient name.
typedef DimensionlessScalar                   Angle;
typedef Quotient<Angle, Time>                 AngularFrequency;
typedef Quotient<AngularFrequency, Time>      AngularAcceleration;
typedef Product<Length, Momentum>             AngularMomentum;
typedef Quotient<AngularMomentum, Time>       Torque;
typedef Quotient<Torque, AngularAcceleration> MomentOfInertia;
#pragma endregion Rotational Motion
typedef Product<Force, Length> Energy;
typedef Quotient<Energy, Time> Power;
typedef Product<Energy, Time>  Action;
#pragma endregion General mechanics
#pragma region Thermodynamics
typedef Product<Length, Length>       Surface;
typedef Product<Surface, Length>      Volume;
typedef Quotient<Force, Surface>      Pressure;
typedef Quotient<Energy, Temperature> Entropy;
typedef Quotient<Mass, Volume>        Density;
typedef DimensionlessScalar           SolidAngle;
#pragma endregion Thermodynamics
}