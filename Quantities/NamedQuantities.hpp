#pragma once

#include "Quantities/Quantities.hpp"

namespace principia {
namespace quantities {
#pragma region General mechanics
typedef Quotient<Length, Time>   Speed;
typedef Quotient<Speed, Time>    Acceleration;
typedef Product<Mass, Speed>     Momentum;
typedef Quotient<Momentum, Time> Force;

typedef Product<Force, Length> Energy;
typedef Quotient<Energy, Time> Power;
typedef Product<Energy, Time>  Action;

typedef Quotient<Angle, Time>                 AngularFrequency;
typedef Quotient<AngularFrequency, Time>      AngularAcceleration;
typedef Quotient<Action, Angle>               AngularMomentum;
typedef Quotient<AngularMomentum, Time>       Torque;
typedef Quotient<Torque, AngularAcceleration> MomentOfInertia;
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
#pragma region Fluid dynamics
typedef Product<Pressure, Time> DynamicViscosity;
typedef Quotient<Area, Time>    KinematicViscosity;
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
#pragma region Electromagnetism
typedef Product<Current, Time>          Charge;
typedef Quotient<Energy, Charge>        Voltage;
typedef Quotient<Charge, Voltage>       Capacitance;
typedef Quotient<Voltage, Current>      Resistance;
typedef Quotient<Current, Voltage>      Conductance;
typedef Quotient<Energy, Current>       MagneticFlux;
typedef Quotient<MagneticFlux, Area>    MagneticFluxDensity;
typedef Quotient<MagneticFlux, Current> Inductance;
typedef Quotient<Force, Charge>         ElectricField;

typedef Product<Quotient<Inductance, Length>, SolidAngle>   Permeability;
typedef Quotient<Quotient<Capacitance, Length>, SolidAngle> Permittivity;

typedef Product<ElectricField, Permittivity> ElectricDisplacementField;
typedef Quotient<MagneticFluxDensity, Permeability> MagneticField;
#pragma endregion
#pragma region Radiometry
typedef Quotient<Power, SolidAngle>      RadiantIntensity;
typedef Quotient<RadiantIntensity, Area> Radiance;
typedef Power                            RadiantFlux;
typedef Product<RadiantFlux, Time>       RadiantEnergy;
typedef Quotient<RadiantFlux, Area>      Irradiance;
typedef Product<Irradiance, Time>        RadiantExposure;
#pragma endregion
#pragma region Photometry
typedef Quotient<LuminousIntensity, Area>      Luminance;
typedef Product<LuminousIntensity, SolidAngle> LuminousFlux;
typedef Product<LuminousFlux, Time>            LuminousEnergy;
typedef Quotient<LuminousFlux, Area>           Illuminance;
typedef Product<Illuminance, Time>             LuminousExposure;
typedef Quotient<LuminousFlux, RadiantFlux>    LuminousEfficacy;
#pragma endregion
}  // namespace quantities
}  // namespace principia
