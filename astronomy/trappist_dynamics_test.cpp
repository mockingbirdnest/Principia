
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

TransitsByPlanet const observations = {
    {"Trappist-1b",
     {2457322.51531, 2457325.53910, 2457328.55860, 2457331.58160, 2457334.60480,
      2457337.62644, 2457340.64820, 2457345.18028, 2457361.79945, 2457364.82173,
      2457440.36492, 2457452.45228, 2457463.02847, 2457509.86460, 2457512.88731,
      2457568.78880, 2457586.91824, 2457589.93922, 2457599.00640, 2457602.02805,
      2457612.60595, 2457615.62710, 2457624.69094, 2457645.84400, 2457651.88743,
      2457653.39809, 2457654.90908, 2457656.41900, 2457657.93129, 2457659.44144,
      2457660.95205, 2457662.46358, 2457663.97492, 2457665.48509, 2457666.99567,
      2457668.50668, 2457670.01766, 2457671.52876, 2457721.38747, 2457739.51770,
      2457741.02787, 2457742.53918, 2457744.05089, 2457745.56164, 2457747.07208,
      2457748.58446, 2457750.09387, 2457751.60535, 2457753.11623, 2457754.62804,
      2457756.13856, 2457757.64840, 2457759.15953, 2457760.67112, 2457762.18120,
      2457763.69221, 2457765.20298, 2457766.71479, 2457768.22514, 2457769.73704,
      2457771.24778, 2457772.75738, 2457774.26841, 2457775.77995, 2457777.28899,
      2457778.80118, 2457780.31297, 2457781.82231, 2457783.33410, 2457784.84372,
      2457792.39979, 2457793.90955, 2457795.41987, 2457796.93134, 2457798.44211,
      2457799.95320, 2457801.46314, 2457802.97557, 2457804.48638, 2457805.99697,
      2457807.50731, 2457809.01822, 2457810.52781, 2457812.04038, 2457813.55121,
      2457815.06275, 2457816.57335, 2457818.08382, 2457819.59478, 2457821.10550,
      2457824.12730, 2457825.63813, 2457827.14995, 2457828.66042, 2457830.17087,
      2457833.19257, 2457834.70398, 2457836.21440, 2457837.72526, 2457839.23669,
      2457917.80060, 2457923.84629, 2457935.93288, 2457952.55450, 2457955.57554,
      2457967.66254, 2457973.70596}},
    {"Trappist-1c",
     {2457333.66400, 2457362.72605, 2457367.57051, 2457384.52320, 2457452.33470,
      2457454.75672, 2457512.88094, 2457546.78587, 2457551.62888, 2457580.69137,
      2457585.53577, 2457587.95622, 2457600.06684, 2457604.90975, 2457609.75461,
      2457614.59710, 2457626.70610, 2457631.55024, 2457638.81518, 2457650.92395,
      2457653.34553, 2457655.76785, 2457658.18963, 2457660.61168, 2457663.03292,
      2457665.45519, 2457667.87729, 2457670.29869, 2457672.71944, 2457711.46778,
      2457723.57663, 2457740.53361, 2457742.95276, 2457745.37429, 2457747.79699,
      2457750.21773, 2457752.64166, 2457755.05877, 2457757.48313, 2457759.90281,
      2457762.32806, 2457764.74831, 2457767.16994, 2457769.59209, 2457772.01483,
      2457774.43458, 2457776.85815, 2457779.27911, 2457781.70095, 2457784.12338,
      2457791.38801, 2457793.81141, 2457796.23153, 2457798.65366, 2457801.07631,
      2457803.49747, 2457805.91882, 2457808.34123, 2457810.76273, 2457813.18456,
      2457815.60583, 2457818.02821, 2457820.45019, 2457822.87188, 2457825.29388,
      2457827.71513, 2457830.13713, 2457832.55888, 2457834.98120, 2457837.40280,
      2457839.82415}},
    {"Trappist-1d",
     {2457625.59779, 2457641.79360, 2457645.84360, 2457653.94261, 2457657.99220,
      2457662.04284, 2457666.09140, 2457670.14198, 2457726.83975, 2457738.99169,
      2457743.03953, 2457747.08985, 2457751.14022, 2457755.18894, 2457759.24638,
      2457763.28895, 2457767.33866, 2457771.39077, 2457775.44026, 2457779.48843,
      2457783.54023, 2457791.64083, 2457803.79083, 2457807.84032, 2457811.89116,
      2457815.94064, 2457819.99050, 2457824.04185, 2457828.09082, 2457832.14036,
      2457836.19171, 2457961.73760, 2457969.83708, 2457973.88590}},
    {"Trappist-1e",
     {2457312.71300, 2457367.59683, 2457611.57620, 2457623.77950, 2457654.27862,
      2457660.38016, 2457666.48030, 2457672.57930, 2457721.37514, 2457733.57300,
      2457739.67085, 2457745.77160, 2457751.87007, 2457757.96712, 2457764.06700,
      2457770.17109, 2457776.26378, 2457782.36226, 2457794.56159, 2457800.66354,
      2457806.75758, 2457812.85701, 2457818.95510, 2457825.05308, 2457831.15206,
      2457837.24980, 2457934.83095, 2457940.92995}},
    {"Trappist-1f",
     {2457321.52520,
      2457367.57629,
      2457634.57809,
      2457652.98579,
      2457662.18747,
      2457671.39279,
      2457717.41541,
      2457726.61960,
      2457745.03116,
      2457754.23380,
      2457763.44338,
      2457772.64752,
      2457781.85142,
      2457800.27307,
      2457809.47554,
      2457818.68271,
      2457827.88669,
      2457837.10322,
      2457956.80549}},
    {"Trappist-1g",
     {2457294.78600,
      2457356.53410,
      2457615.92400,
      2457640.63730,
      2457652.99481,
      2457665.35151,
      2457739.48441,
      2457751.83993,
      2457764.19098,
      2457776.54900,
      2457801.25000,
      2457813.60684,
      2457825.96112,
      2457838.30655,
      2457924.77090,
      2457961.82621}},
    {"Trappist-1h",
     {2457662.55467,
      2457756.38740,
      2457775.15390,
      2457793.92300,
      2457812.69870,
      2457831.46625,
      2457962.86271}}};

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
