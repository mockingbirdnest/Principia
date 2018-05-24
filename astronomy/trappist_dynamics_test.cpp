
#include "astronomy/frames.hpp"
#include "base/file.hpp"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/sign.hpp"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/root_finders.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/si.hpp"

namespace principia {

using base::not_null;
using base::OFStream;
using geometry::Instant;
using geometry::Position;
using geometry::Sign;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Quinlan1999Order8A;
using numerics::Bisect;
using physics::Ephemeris;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::RelativeDegreesOfFreedom;
using physics::SolarSystem;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;

namespace astronomy {

using Transits = std::vector<Instant>;
using TransitsByPlanet = std::map<std::string, Transits>;

class TrappistDynamicsTest : public ::testing::Test {
 protected:
  TrappistDynamicsTest()
      : system_(SOLUTION_DIR / "astronomy" / "trappist_gravity_model.proto.txt",
                SOLUTION_DIR / "astronomy" /
                    "trappist_initial_state_jd_2457282_805700000.proto.txt"),
        ephemeris_(system_.MakeEphemeris(
            /*fitting_tolerance=*/5 * Milli(Metre),
            Ephemeris<Trappist>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                                   Position<Trappist>>(),
                /*step=*/0.07 * Day))) {}

  static std::string SanitizedName(MassiveBody const& body) {
    auto sanitized_name = body.name();
    return sanitized_name.erase(sanitized_name.find_first_of("-"), 1);
  }

  static Time MaxError(TransitsByPlanet const& observations,
                       TransitsByPlanet const& computations) {
    Time max_error;
    for (auto const& pair : observations) {
      auto const& name = pair.first;
      auto const& observed_transits = pair.second;
      auto const& computed_transits = computations.at(name);
      for (auto const& observed_transit : observed_transits) {
        auto const next_computed_transit =
            std::lower_bound(computed_transits.begin(),
                             computed_transits.end(),
                             observed_transit);
        Time error;
        if (next_computed_transit == computed_transits.begin()) {
          error = *next_computed_transit - observed_transit;
        } else if (next_computed_transit == computed_transits.end()) {
          error = observed_transit - computed_transits.back();
        } else {
          error =
              std::min(*next_computed_transit - observed_transit,
                       observed_transit - *std::prev(next_computed_transit));
        }
        CHECK_LE(0.0 * Second, error);
        max_error = std::max(max_error, error);
      }
    }
    return max_error;
  }

  SolarSystem<Trappist> const system_;
  not_null<std::unique_ptr<Ephemeris<Trappist>>> ephemeris_;
};

TEST_F(TrappistDynamicsTest, MathematicaPeriods) {
  Instant const a_century_later = system_.epoch() + 100 * JulianYear;
  ephemeris_->Prolong(a_century_later);

  auto const& star = system_.massive_body(*ephemeris_, "Trappist-1A");
  auto const& star_trajectory = ephemeris_->trajectory(star);

  OFStream file(TEMP_DIR / "trappist_periods.generated.wl");
  auto const bodies = ephemeris_->bodies();
  for (auto const& planet : bodies) {
    if (planet != star) {
      auto const& planet_trajectory = ephemeris_->trajectory(planet);
      std::vector<Time> periods;
      for (Instant t = ephemeris_->t_max() - 2000 * Hour;
           t < ephemeris_->t_max();
           t += 1 * Hour) {
        KeplerOrbit<Trappist> const planet_orbit(
            *star,
            *planet,
            planet_trajectory->EvaluateDegreesOfFreedom(t) -
                star_trajectory->EvaluateDegreesOfFreedom(t),
            t);
        periods.push_back(*planet_orbit.elements_at_epoch().period);
      }

      file << mathematica::Assign("period" + SanitizedName(*planet),
                                  periods);
    }
  }
}

TEST_F(TrappistDynamicsTest, MathematicaTransits) {
  Instant const a_century_later = system_.epoch() + 100 * JulianYear;
  ephemeris_->Prolong(a_century_later);

  auto const& star = system_.massive_body(*ephemeris_, "Trappist-1A");
  auto const& star_trajectory = ephemeris_->trajectory(star);

  TransitsByPlanet transits_by_planet;

  OFStream file(TEMP_DIR / "trappist_transits.generated.wl");
  auto const bodies = ephemeris_->bodies();
  for (auto const& planet : bodies) {
    if (planet != star) {
      auto const& planet_trajectory = ephemeris_->trajectory(planet);

      Transits transits;
      std::optional<Instant> last_t;
      std::optional<Sign> last_xy_displacement_derivative_sign;
      for (Instant t = ephemeris_->t_min();
           t < ephemeris_->t_min() + 10 * JulianYear;
           t += 2 * Hour) {
        RelativeDegreesOfFreedom<Trappist> const relative_dof =
            planet_trajectory->EvaluateDegreesOfFreedom(t) -
            star_trajectory->EvaluateDegreesOfFreedom(t);

        auto const xy_displacement_derivative =
            [&planet_trajectory, &star_trajectory](Instant const& t) {
              RelativeDegreesOfFreedom<Trappist> const relative_dof =
                  planet_trajectory->EvaluateDegreesOfFreedom(t) -
                  star_trajectory->EvaluateDegreesOfFreedom(t);
              // TODO(phl): Why don't we have projections?
              auto xy_displacement =
                  relative_dof.displacement().coordinates();
              xy_displacement.z = 0.0 * Metre;
              auto xy_velocity = relative_dof.velocity().coordinates();
              xy_velocity.z = 0.0 * Metre / Second;
              return Dot(xy_displacement, xy_velocity);
            };

        Sign const xy_displacement_derivative_sign(
            xy_displacement_derivative(t));
        if (relative_dof.displacement().coordinates().z > 0.0 * Metre &&
            last_t &&
            xy_displacement_derivative_sign == Sign(1) &&
            last_xy_displacement_derivative_sign == Sign(-1)) {
          Instant const transit =
              Bisect(xy_displacement_derivative, *last_t, t);
          transits_by_planet[planet->name()].push_back(transit);
        }
        last_t = t;
        last_xy_displacement_derivative_sign =
            xy_displacement_derivative_sign;
      }

      file << mathematica::Assign("transit" + SanitizedName(*planet),
                                  transits);
    }
  }
}

}  // namespace astronomy
}  // namespace principia
