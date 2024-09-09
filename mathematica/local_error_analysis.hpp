#pragma once

#include <filesystem>
#include <memory>

#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "integrators/integrators.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace mathematica {
namespace _local_error_analysis {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::integrators::_integrators;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

// A utility to compute the local errors in the numerical integration of a
// `solar_system` with a given `integrator` and `step`.
template<typename Frame>
class LocalErrorAnalyser {
 public:
  LocalErrorAnalyser(
      not_null<std::unique_ptr<SolarSystem<Frame>>> solar_system,
      FixedStepSizeIntegrator<typename
          Ephemeris<Frame>::NewtonianMotionEquation> const& integrator,
      Time const& step);

  // Computes the error over `granularity` between the main integration and a
  // fine integration forked off the main one, for `duration` from the solar
  // system epoch.  Writes the errors to a file with the given `path`.
  void WriteLocalErrors(
      std::filesystem::path const& path,
      FixedStepSizeIntegrator<typename
          Ephemeris<Frame>::NewtonianMotionEquation> const& fine_integrator,
      Time const& fine_step,
      Time const& granularity,
      Time const& duration) const;

 private:
  not_null<std::unique_ptr<Ephemeris<Frame>>> ForkEphemeris(
      Ephemeris<Frame> const& original,
      Instant const& t,
      FixedStepSizeIntegrator<typename
          Ephemeris<Frame>::NewtonianMotionEquation> const& integrator,
      Time const& step) const;

  static constexpr Length fitting_tolerance_ = 1 * Milli(Metre);

  not_null<std::unique_ptr<SolarSystem<Frame>>> const solar_system_;
  FixedStepSizeIntegrator<
      typename Ephemeris<Frame>::NewtonianMotionEquation> const& integrator_;
  Time const step_;
};

}  // namespace internal

using internal::LocalErrorAnalyser;

}  // namespace _local_error_analysis
}  // namespace mathematica
}  // namespace principia

#include "mathematica/local_error_analysis_body.hpp"
