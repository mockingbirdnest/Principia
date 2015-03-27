#pragma once

#include "quantities/quantities.hpp"

namespace principia {

namespace quantities {

// |Variation<T>| is the type of the time derivative of a |T|-valued function.
template<typename T>
using Variation = Quotient<Difference<T, T>, Time>;

// General mechanics
using Speed        = Variation<Length>;
using Acceleration = Variation<Speed>;
using Momentum     = Product<Mass, Speed>;
using Force        = Variation<Momentum>;
using Stiffness    = Quotient<Force, Length>;

using Energy = Product<Force, Length>;
using Power  = Variation<Energy>;
using Action = Product<Energy, Time>;

using AngularFrequency    = Variation<Angle>;
using AngularAcceleration = Variation<AngularFrequency>;
using AngularMomentum     = Quotient<Action, Angle>;
using Torque              = Variation<AngularMomentum>;
using MomentOfInertia     = Quotient<Torque, AngularAcceleration>;

using GravitationalParameter = Product<Length, Product<Speed, Speed>>;
using Order2ZonalCoefficient = Quotient<Exponentiation<Length, 5>,
                                        Exponentiation<Time, 2>>;

// Astrodynamics
using SpecificEnergy          = Quotient<Energy, Mass>;
using SpecificAngularMomentum = Quotient<AngularMomentum, Mass>;

// Thermodynamics
using Area           = Product<Length, Length>;
using Volume         = Product<Area, Length>;
using Pressure       = Quotient<Force, Area>;
using Entropy        = Quotient<Energy, Temperature>;
using Density        = Quotient<Mass, Volume>;
using SpecificVolume = Quotient<Volume, Mass>;
using MolarVolume    = Quotient<Volume, Amount>;

// Fluid dynamics
using DynamicViscosity   = Product<Pressure, Time>;
using KinematicViscosity = Quotient<Area, Time>;

// Chemistry
using Concentration     = Quotient<Amount, Volume>;
using MolarMass         = Quotient<Mass, Amount>;
using CatalyticActivity = Quotient<Amount, Time>;

// Optics
using Wavenumber = Quotient<Angle, Length>;

// Spectroscopy
// Nonstandard, the SI considers cycles as dimensionless. This is annoying
// because of the resulting hopeless confusion between frequency and angular
// frequency. We choose to strongly type cycles.
using Frequency               = Quotient<Winding, Time>;
using Wavelength              = Quotient<Length, Winding>;
using SpectroscopicWavenumber = Quotient<Winding, Length>;

// Electromagnetism
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

// Radiometry
using RadiantIntensity = Quotient<Power, SolidAngle>;
using Radiance         = Quotient<RadiantIntensity, Area>;
using RadiantFlux      = Power;
using RadiantEnergy    = Product<RadiantFlux, Time>;
using Irradiance       = Quotient<RadiantFlux, Area>;
using RadiantExposure  = Product<Irradiance, Time>;

// Photometry
using Luminance        = Quotient<LuminousIntensity, Area>;
using LuminousFlux     = Product<LuminousIntensity, SolidAngle>;
using LuminousEnergy   = Product<LuminousFlux, Time>;
using Illuminance      = Quotient<LuminousFlux, Area>;
using LuminousExposure = Product<Illuminance, Time>;
using LuminousEfficacy = Quotient<LuminousFlux, RadiantFlux>;
}  // namespace quantities
}  // namespace principia
