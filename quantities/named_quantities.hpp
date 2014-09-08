#pragma once

#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
#pragma region General mechanics
using Speed        = Quotient<Length, Time>;
using Acceleration = Quotient<Speed, Time>;
using Momentum     = Product<Mass, Speed>;
using Force        = Quotient<Momentum, Time>;
using Stiffness    = Quotient<Force, Length>;

using Energy = Product<Force, Length>;
using Power  = Quotient<Energy, Time>;
using Action = Product<Energy, Time>;

using AngularFrequency    = Quotient<Angle, Time>;
using AngularAcceleration = Quotient<AngularFrequency, Time>;
using AngularMomentum     = Quotient<Action, Angle>;
using Torque              = Quotient<AngularMomentum, Time>;
using MomentOfInertia     = Quotient<Torque, AngularAcceleration>;

using GravitationalParameter = Product<Length, Product<Speed, Speed>>;
#pragma endregion
#pragma region Astrodynamics
using SpecificEnergy          = Quotient<Energy, Mass>;
using SpecificAngularMomentum = Quotient<AngularMomentum, Mass>;
#pragma endregion
#pragma region Thermodynamics
using Area           = Product<Length, Length>;
using Volume         = Product<Area, Length>;
using Pressure       = Quotient<Force, Area>;
using Entropy        = Quotient<Energy, Temperature>;
using Density        = Quotient<Mass, Volume>;
using SpecificVolume = Quotient<Volume, Mass>;
using MolarVolume    = Quotient<Volume, Amount>;
#pragma endregion
#pragma region Fluid dynamics
using DynamicViscosity   = Product<Pressure, Time>;
using KinematicViscosity = Quotient<Area, Time>;
#pragma endregion
#pragma region Chemistry
using Concentration     = Quotient<Amount, Volume>;
using MolarMass         = Quotient<Mass, Amount>;
using CatalyticActivity = Quotient<Amount, Time>;
#pragma endregion
#pragma region Optics
using Wavenumber = Quotient<Angle, Length>;
#pragma endregion
#pragma region Spectroscopy
// Nonstandard, the SI considers cycles as dimensionless. This is annoying
// because of the resulting hopeless confusion between frequency and angular
// frequency. We choose to strongly type cycles.
using Frequency               = Quotient<Winding, Time>;
using Wavelength              = Quotient<Length, Winding>;
using SpectroscopicWavenumber = Quotient<Winding, Length>;
#pragma endregion
#pragma region Electromagnetism
using Charge              = Product<Current, Time>;
using Voltage             = Quotient<Energy, Charge>;
using Capacitance         = Quotient<Charge, Voltage>;
using Resistance          = Quotient<Voltage, Current>;
using Conductance         = Quotient<Current, Voltage>;
using MagneticFlux        = Quotient<Energy, Current>;
using MagneticFluxDensity = Quotient<MagneticFlux, Area>;
using Inductance          = Quotient<MagneticFlux, Current>;
using ElectricField       = Quotient<Force, Charge>;

using Permeability = Product<Quotient<Inductance, Length>, SolidAngle>;
using Permittivity = Quotient<Quotient<Capacitance, Length>, SolidAngle>;

using ElectricDisplacementField = Product<ElectricField, Permittivity>;
using MagneticField             = Quotient<MagneticFluxDensity, Permeability>;
#pragma endregion
#pragma region Radiometry
using RadiantIntensity = Quotient<Power, SolidAngle>;
using Radiance         = Quotient<RadiantIntensity, Area>;
using RadiantFlux      = Power;
using RadiantEnergy    = Product<RadiantFlux, Time>;
using Irradiance       = Quotient<RadiantFlux, Area>;
using RadiantExposure  = Product<Irradiance, Time>;
#pragma endregion
#pragma region Photometry
using Luminance        = Quotient<LuminousIntensity, Area>;
using LuminousFlux     = Product<LuminousIntensity, SolidAngle>;
using LuminousEnergy   = Product<LuminousFlux, Time>;
using Illuminance      = Quotient<LuminousFlux, Area>;
using LuminousExposure = Product<Illuminance, Time>;
using LuminousEfficacy = Quotient<LuminousFlux, RadiantFlux>;
#pragma endregion
}  // namespace quantities
}  // namespace principia
