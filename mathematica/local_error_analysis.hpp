
#pragma once

#include "integrators/integrators.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"

namespace principia {
namespace mathematica {

using astronomy::ICRFJ2000Equator;
using base::not_null;
using geometry::Instant;
using integrators::FixedStepSizeIntegrator;
using physics::Ephemeris;
using physics::SolarSystem;
using quantities::Time;

class LocalErrorAnalyser {
 public:
  LocalErrorAnalyser(
      not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>> solar_system,
      FixedStepSizeIntegrator<
          Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const&
          integrator,
      Time const& step);

  void WriteDailyErrors(
      std::experimental::filesystem::path path,
      FixedStepSizeIntegrator<
          Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const&
          fine_integrator,
      Time const& fine_step,
      Time const& granularity,
      Time const& duration) const;

 private:
  not_null<std::unique_ptr<Ephemeris<ICRFJ2000Equator>>> ForkEphemeris(
      Ephemeris<ICRFJ2000Equator> const& original,
      Instant const& t,
      FixedStepSizeIntegrator<
          Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const&
          integrator,
      Time const& step) const;

  not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>> const solar_system_;
  FixedStepSizeIntegrator<
      Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const& integrator_;
  Time const step_;
};

}  // namespace mathematica
}  // namespace principia
