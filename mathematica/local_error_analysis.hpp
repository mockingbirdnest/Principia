
#pragma once

#include <filesystem>

#include "integrators/integrators.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"

namespace principia {
namespace mathematica {

using astronomy::ICRS;
using base::not_null;
using geometry::Instant;
using integrators::FixedStepSizeIntegrator;
using physics::Ephemeris;
using physics::SolarSystem;
using quantities::Time;

// A utility to compute the local errors in the numerical integration of a
// |solar_system| with a given |integrator| and |step|.
class LocalErrorAnalyser {
 public:
  LocalErrorAnalyser(
      not_null<std::unique_ptr<SolarSystem<ICRS>>> solar_system,
      FixedStepSizeIntegrator<
          Ephemeris<ICRS>::NewtonianMotionEquation> const&
          integrator,
      Time const& step);

  // Computes the error over |granularity| between the main integration and a
  // fine integration forked off the main one, for |duration| from the solar
  // system epoch.  Writes the errors to a file with the given |path|.
  void WriteLocalErrors(
      std::filesystem::path const& path,
      FixedStepSizeIntegrator<
          Ephemeris<ICRS>::NewtonianMotionEquation> const&
          fine_integrator,
      Time const& fine_step,
      Time const& granularity,
      Time const& duration) const;

 private:
  not_null<std::unique_ptr<Ephemeris<ICRS>>> ForkEphemeris(
      Ephemeris<ICRS> const& original,
      Instant const& t,
      FixedStepSizeIntegrator<
          Ephemeris<ICRS>::NewtonianMotionEquation> const&
          integrator,
      Time const& step) const;

  not_null<std::unique_ptr<SolarSystem<ICRS>>> const solar_system_;
  FixedStepSizeIntegrator<
      Ephemeris<ICRS>::NewtonianMotionEquation> const& integrator_;
  Time const step_;
};

}  // namespace mathematica
}  // namespace principia
