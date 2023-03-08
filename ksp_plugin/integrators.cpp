#include "ksp_plugin/integrators.hpp"

#include <limits>

#include "geometry/named_quantities.hpp"
#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

namespace principia {
namespace ksp_plugin {
namespace _integrators {
namespace internal {

using namespace principia::geometry::_named_quantities;
using namespace principia::integrators::
    _embedded_explicit_generalized_runge_kutta_nyström_integrator;
using namespace principia::integrators::
    _embedded_explicit_runge_kutta_nyström_integrator;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::integrators::
    _symplectic_runge_kutta_nyström_integrator;
using namespace principia::quantities::_si;

DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters
DefaultDownsamplingParameters() {
  return DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters{
      .max_dense_intervals = 10'000,
      .tolerance = 10 * Metre,
  };
}

DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters
OrbitAnalyserDownsamplingParameters() {
  // Carefully tuned based on MercuryOrbiter test.
  return DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters{
      .max_dense_intervals = 10'000,
      .tolerance = 1 * Milli(Metre),
  };
}

Ephemeris<Barycentric>::AccuracyParameters
DefaultEphemerisAccuracyParameters() {
  return Ephemeris<Barycentric>::AccuracyParameters(
      /*fitting_tolerance=*/1 * Milli(Metre),
      /*geopotential_tolerance*/ 0x1.0p-24);
}

Ephemeris<Barycentric>::FixedStepParameters
DefaultEphemerisFixedStepParameters() {
  return Ephemeris<Barycentric>::FixedStepParameters(
      SymplecticRungeKuttaNyströmIntegrator<
          BlanesMoan2002SRKN14A,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      /*step=*/35 * Minute);
}

Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
DefaultBurnParameters() {
  return Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters(
      EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
          Fine1987RKNG34,
          Ephemeris<Barycentric>::GeneralizedNewtonianMotionEquation>(),
      /*max_steps=*/1000,
      /*length_integration_tolerance=*/1 * Metre,
      /*speed_integration_tolerance=*/1 * Metre / Second);
}

Ephemeris<Barycentric>::FixedStepParameters DefaultHistoryParameters() {
  return Ephemeris<Barycentric>::FixedStepParameters(
      SymmetricLinearMultistepIntegrator<
          Quinlan1999Order8A,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      /*step=*/10 * Second);
}

Ephemeris<Barycentric>::AdaptiveStepParameters
DefaultPsychohistoryParameters() {
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          DormandالمكاوىPrince1986RKN434FM,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*length_integration_tolerance=*/1 * Milli(Metre),
      /*speed_integration_tolerance=*/1 * Milli(Metre) / Second);
}

Ephemeris<Barycentric>::AdaptiveStepParameters DefaultPredictionParameters() {
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          DormandالمكاوىPrince1986RKN434FM,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      /*max_steps=*/1000,
      /*length_integration_tolerance=*/1 * Metre,
      /*speed_integration_tolerance=*/1 * Metre / Second);
}

}  // namespace internal
}  // namespace _integrators
}  // namespace ksp_plugin
}  // namespace principia
