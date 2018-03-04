
#include "mathematica/local_error_analysis.hpp"

#include "astronomy/stabilize_ksp.hpp"
#include "base/file.hpp"
#include "mathematica/mathematica.hpp"

namespace principia {
namespace mathematica {

using base::OFStream;
using base::make_not_null_unique;
using geometry::Position;
using physics::DegreesOfFreedom;
using physics::MassiveBody;
using quantities::Length;
using quantities::si::Day;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;

LocalErrorAnalyser::LocalErrorAnalyser(
    not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>> solar_system,
    FixedStepSizeIntegrator<
        Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const& integrator,
    Time const& step)
    : solar_system_(std::move(solar_system)),
      integrator_(integrator),
      step_(step) {
  // The system might not be defined from Keplerian elements, so we cannot
  // always turn it into a hierarchical system to take its fingerprint.
  // TODO(eggrobin): arbitrary solar system fingerprinting.
  if (solar_system_->names()[0] == "Bop") {
    LOG(INFO) << "All hail retrobop!";
    astronomy::StabilizeKSP(*solar_system_);
  }
}

void LocalErrorAnalyser::WriteDailyErrors(
    std::experimental::filesystem::path path) const {
  auto const reference_ephemeris = solar_system_->MakeEphemeris(
      /*fitting_tolerance=*/1 * Milli(Metre),
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(integrator_, step_));
  reference_ephemeris->Prolong(solar_system_->epoch());
  std::vector<std::vector<Length>> errors;
  for (int day = 1; day < 500; ++day) {
    Instant const t0 = solar_system_->epoch() + (day - 1) * Day;
    std::unique_ptr<Ephemeris<ICRFJ2000Equator>> refined_ephemeris =
        ForkEphemeris(
            *reference_ephemeris,
            t0,
            integrators::BlanesMoan2002SRKN14A<Position<ICRFJ2000Equator>>(),
            1 * Minute);
    Instant const t = solar_system_->epoch() + day * Day;
    reference_ephemeris->Prolong(t);
    refined_ephemeris->Prolong(t);
    LOG_IF(INFO, day % 10 == 0) << "Prolonged to " << day << " days.";

    errors.emplace_back();
    for (auto const& body_name : solar_system_->names()) {
      int const body_index = solar_system_->index(body_name);
      errors.back().push_back(
          (reference_ephemeris
               ->trajectory(reference_ephemeris->bodies()[body_index])
               ->EvaluatePosition(t) -
           refined_ephemeris
               ->trajectory(refined_ephemeris->bodies()[body_index])
               ->EvaluatePosition(t)).Norm());
    }
  }
  OFStream file(path);
  file << Assign("bodyNames", solar_system_->names());
  file << Assign("errors", ExpressIn(Metre, errors));
}

not_null<std::unique_ptr<Ephemeris<ICRFJ2000Equator>>>
LocalErrorAnalyser::ForkEphemeris(
    Ephemeris<ICRFJ2000Equator> const& original,
    Instant const& t,
    FixedStepSizeIntegrator<
        Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const& integrator,
    Time const& step) const {
  std::vector<DegreesOfFreedom<ICRFJ2000Equator>> degrees_of_freedom;
  for (not_null<MassiveBody const*> const body : original.bodies()) {
    degrees_of_freedom.emplace_back(
        original.trajectory(body)->EvaluateDegreesOfFreedom(t));
  }
  return make_not_null_unique<Ephemeris<ICRFJ2000Equator>>(
      solar_system_->MakeAllMassiveBodies(),
      degrees_of_freedom,
      t,
      /*fitting_tolerance=*/1 * Milli(Metre),
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(integrator, step));
}

}  // namespace mathematica
}  // namespace principia
