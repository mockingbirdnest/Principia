// NamedQuantities.h

#pragma once

#include "PhysicalQuantities.h"

namespace PhysicalQuantities {
#pragma region Dimensionless quantities
// There is no good way of strongly typing angles vs. dimensionless scalars
// because angles *are* dimensionless scalars. These are just convenient names.
typedef Dimensionless Angle;
typedef Dimensionless SolidAngle;
#pragma endregion
#pragma region General mechanics
#pragma region Linear motion
typedef Quotient<Length, Time>   Speed;
typedef Quotient<Speed, Time>    Acceleration;
typedef Product<Mass, Speed>     Momentum;
typedef Quotient<Momentum, Time> Force;
#pragma endregion
#pragma region Rotational motion
typedef Quotient<Angle, Time>                 AngularFrequency;
typedef Quotient<AngularFrequency, Time>      AngularAcceleration;
typedef Product<Length, Momentum>             AngularMomentum;
typedef Quotient<AngularMomentum, Time>       Torque;
typedef Quotient<Torque, AngularAcceleration> MomentOfInertia;
#pragma endregion
typedef Product<Force, Length> Energy;
typedef Quotient<Energy, Time> Power;
typedef Product<Energy, Time>  Action;
#pragma endregion
#pragma region Thermodynamics
typedef Product<Length, Length>       Area;
typedef Product<Area, Length>         Volume;
typedef Quotient<Force, Area>         Pressure;
typedef Quotient<Energy, Temperature> Entropy;
typedef Quotient<Mass, Volume>        Density;
typedef Quotient<Volume, Mass>        SpecificVolume;
typedef Quotient<Volume, Amount>      MolarVolume;
#pragma endregion
#pragma region Chemistry
typedef Quotient<Amount, Volume> Concentration;
typedef Quotient<Mass, Amount>   MolarMass;
typedef Quotient<Amount, Time>   CatalyticActivity;
#pragma endregion
#pragma region Optics
typedef Quotient<Angle, Length> Wavenumber;
#pragma endregion
#pragma region Spectroscopy
// Nonstandard, the SI considers cycles as dimensionless. This is annoying
// because of the resulting hopeless confusion between frequency and angular
// frequency. We choose to strongly type cycles.
typedef Quotient<Winding, Time>   Frequency;
typedef Quotient<Length, Winding> Wavelength;
typedef Quotient<Winding, Length> SpectroscopicWavenumber;
#pragma endregion
#pragma region Electricity
typedef Product<Current, Time>          Charge;
typedef Quotient<Energy, Charge>        Voltage;
typedef Quotient<Charge, Voltage>       Capacitance;
typedef Quotient<Voltage, Current>      Resistance;
typedef Quotient<Current, Voltage>      Conductance;
typedef Quotient<Energy, Current>       MagneticFlux;
typedef Quotient<MagneticFlux, Area>    MagneticFluxDensity;
typedef Quotient<MagneticFlux, Current> Inductance;
typedef Quotient<Inductance, Length>    Permeability;
typedef Quotient<Capacitance, Length>   Permittivity;
#pragma endregion
#pragma region Radiometry
typedef Quotient<Power, SolidAngle>      RadiantIntensity;
typedef Quotient<RadiantIntensity, Area> Radiance;
typedef Quotient<Power, Wrapping>        RadiantFlux;
typedef Product<RadiantFlux, Time>       RadiantEnergy;
typedef Quotient<RadiantFlux, Area>      Irradiance;
typedef Product<Irradiance, Time>        RadiantExposure;
#pragma endregion
#pragma region Photometry
typedef Quotient<LuminousIntensity, Area>     Luminance;
typedef Quotient<LuminousIntensity, Wrapping> LuminousFlux;
typedef Product<LuminousFlux, Time>           LuminousEnergy;
typedef Quotient<LuminousFlux, Area>          Illuminance;
typedef Product<Illuminance, Time>            LuminousExposure;
typedef Quotient<LuminousFlux, RadiantFlux>   LuminousEfficacy;
#pragma endregion
}