
#pragma once

#include "mathematica/local_error_analysis.hpp"

#include <filesystem>
#include <vector>

#include "astronomy/stabilize_ksp.hpp"
#include "astronomy/solar_system_fingerprints.hpp"
#include "base/file.hpp"
#include "mathematica/mathematica.hpp"

namespace principia {
namespace mathematica {

using astronomy::KSP191;
using astronomy::KSPStabilizedSystemFingerprints;
using astronomy::KSPStockSystemFingerprints;
using base::OFStream;
using base::make_not_null_unique;
using geometry::Position;
using physics::DegreesOfFreedom;
using physics::MassiveBody;
using quantities::si::Day;
using quantities::si::Minute;

template<typename Frame>
LocalErrorAnalyser<Frame>::LocalErrorAnalyser(
    not_null<std::unique_ptr<SolarSystem<Frame>>> solar_system,
    FixedStepSizeIntegrator<typename
        Ephemeris<Frame>::NewtonianMotionEquation> const& integrator,
    Time const& step)
    : solar_system_(std::move(solar_system)),
      integrator_(integrator),
      step_(step) {
  if (solar_system_->Fingerprint() == KSPStockSystemFingerprints[KSP191]) {
    LOG(INFO) << "All hail retrobop!";
    astronomy::StabilizeKSP(*solar_system_);
    CHECK_EQ(solar_system_->Fingerprint(),
             KSPStabilizedSystemFingerprints[KSP191]);
  }
}

template<typename Frame>
void LocalErrorAnalyser<Frame>::WriteLocalErrors(
    std::filesystem::path const& path,
    FixedStepSizeIntegrator<typename
        Ephemeris<Frame>::NewtonianMotionEquation> const& fine_integrator,
    Time const& fine_step,
    Time const& granularity,
    Time const& duration) const {
  auto const reference_ephemeris = solar_system_->MakeEphemeris(
      /*accuracy_parameters=*/{fitting_tolerance_,
                               /*geopotential_tolerance=*/0x1p-24},
      typename Ephemeris<Frame>::FixedStepParameters(integrator_, step_));
  reference_ephemeris->Prolong(solar_system_->epoch());
  std::vector<std::vector<Length>> errors;
  for (Instant t0 = solar_system_->epoch(),
               t = t0 + granularity;
       t < solar_system_->epoch() + duration;
       t0 = t, t += granularity) {
    std::unique_ptr<Ephemeris<Frame>> refined_ephemeris =
        ForkEphemeris(*reference_ephemeris, t0, fine_integrator, fine_step);
    reference_ephemeris->Prolong(t);
    refined_ephemeris->Prolong(t);
    LOG_EVERY_N(INFO, 10) << "Prolonged to "
                          << (t - solar_system_->epoch()) / Day << " days.";

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
  file << Assign("errors", errors, ExpressIn(Metre));
}

template<typename Frame>
not_null<std::unique_ptr<Ephemeris<Frame>>>
LocalErrorAnalyser<Frame>::ForkEphemeris(
    Ephemeris<Frame> const& original,
    Instant const& t,
    FixedStepSizeIntegrator<typename
        Ephemeris<Frame>::NewtonianMotionEquation> const& integrator,
    Time const& step) const {
  std::vector<DegreesOfFreedom<Frame>> degrees_of_freedom;
  for (not_null<MassiveBody const*> const body : original.bodies()) {
    degrees_of_freedom.emplace_back(
        original.trajectory(body)->EvaluateDegreesOfFreedom(t));
  }
  return make_not_null_unique<Ephemeris<Frame>>(
      solar_system_->MakeAllMassiveBodies(),
      degrees_of_freedom,
      t,
      typename Ephemeris<Frame>::AccuracyParameters(
          fitting_tolerance_, /*geopotential_tolerance=*/0x1p-24),
      typename Ephemeris<Frame>::FixedStepParameters(integrator, step));
}

}  // namespace mathematica
}  // namespace principia
