
#include "ksp_plugin/integrators.hpp"

#include <limits>

#include "geometry/named_quantities.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_integrators {

using geometry::Position;
using integrators::BlanesMoan2002SRKN14A;
using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::Quinlan1999Order8A;
using quantities::si::Minute;
using quantities::si::Second;

Ephemeris<Barycentric>::FixedStepParameters DefaultEphemerisParameters() {
  return Ephemeris<Barycentric>::FixedStepParameters(
             BlanesMoan2002SRKN14A<Position<Barycentric>>(),
             /*step=*/35 * Minute);
}

Ephemeris<Barycentric>::FixedStepParameters DefaultHistoryParameters() {
  return Ephemeris<Barycentric>::FixedStepParameters(
             Quinlan1999Order8A<Position<Barycentric>>(),
             /*step=*/10 * Second);
}

Ephemeris<Barycentric>::AdaptiveStepParameters
DefaultPsychohistoryParameters() {
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
             DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
             /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
             /*length_integration_tolerance=*/1 * Milli(Metre),
             /*speed_integration_tolerance=*/1 * Milli(Metre) / Second);
}

Ephemeris<Barycentric>::AdaptiveStepParameters DefaultPredictionParameters() {
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
             DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
             /*max_steps=*/1000,
             /*length_integration_tolerance=*/1 * Metre,
             /*speed_integration_tolerance=*/1 * Metre / Second);
}

}  // namespace internal_integrators
}  // namespace ksp_plugin
}  // namespace principia
