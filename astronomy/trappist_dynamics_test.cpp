
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
     {"JD2457322.51531"_TT, "JD2457325.53910"_TT, "JD2457328.55860"_TT,
      "JD2457331.58160"_TT, "JD2457334.60480"_TT, "JD2457337.62644"_TT,
      "JD2457340.64820"_TT, "JD2457345.18028"_TT, "JD2457361.79945"_TT,
      "JD2457364.82173"_TT, "JD2457440.36492"_TT, "JD2457452.45228"_TT,
      "JD2457463.02847"_TT, "JD2457509.86460"_TT, "JD2457512.88731"_TT,
      "JD2457568.78880"_TT, "JD2457586.91824"_TT, "JD2457589.93922"_TT,
      "JD2457599.00640"_TT, "JD2457602.02805"_TT, "JD2457612.60595"_TT,
      "JD2457615.62710"_TT, "JD2457624.69094"_TT, "JD2457645.84400"_TT,
      "JD2457651.88743"_TT, "JD2457653.39809"_TT, "JD2457654.90908"_TT,
      "JD2457656.41900"_TT, "JD2457657.93129"_TT, "JD2457659.44144"_TT,
      "JD2457660.95205"_TT, "JD2457662.46358"_TT, "JD2457663.97492"_TT,
      "JD2457665.48509"_TT, "JD2457666.99567"_TT, "JD2457668.50668"_TT,
      "JD2457670.01766"_TT, "JD2457671.52876"_TT, "JD2457721.38747"_TT,
      "JD2457739.51770"_TT, "JD2457741.02787"_TT, "JD2457742.53918"_TT,
      "JD2457744.05089"_TT, "JD2457745.56164"_TT, "JD2457747.07208"_TT,
      "JD2457748.58446"_TT, "JD2457750.09387"_TT, "JD2457751.60535"_TT,
      "JD2457753.11623"_TT, "JD2457754.62804"_TT, "JD2457756.13856"_TT,
      "JD2457757.64840"_TT, "JD2457759.15953"_TT, "JD2457760.67112"_TT,
      "JD2457762.18120"_TT, "JD2457763.69221"_TT, "JD2457765.20298"_TT,
      "JD2457766.71479"_TT, "JD2457768.22514"_TT, "JD2457769.73704"_TT,
      "JD2457771.24778"_TT, "JD2457772.75738"_TT, "JD2457774.26841"_TT,
      "JD2457775.77995"_TT, "JD2457777.28899"_TT, "JD2457778.80118"_TT,
      "JD2457780.31297"_TT, "JD2457781.82231"_TT, "JD2457783.33410"_TT,
      "JD2457784.84372"_TT, "JD2457792.39979"_TT, "JD2457793.90955"_TT,
      "JD2457795.41987"_TT, "JD2457796.93134"_TT, "JD2457798.44211"_TT,
      "JD2457799.95320"_TT, "JD2457801.46314"_TT, "JD2457802.97557"_TT,
      "JD2457804.48638"_TT, "JD2457805.99697"_TT, "JD2457807.50731"_TT,
      "JD2457809.01822"_TT, "JD2457810.52781"_TT, "JD2457812.04038"_TT,
      "JD2457813.55121"_TT, "JD2457815.06275"_TT, "JD2457816.57335"_TT,
      "JD2457818.08382"_TT, "JD2457819.59478"_TT, "JD2457821.10550"_TT,
      "JD2457824.12730"_TT, "JD2457825.63813"_TT, "JD2457827.14995"_TT,
      "JD2457828.66042"_TT, "JD2457830.17087"_TT, "JD2457833.19257"_TT,
      "JD2457834.70398"_TT, "JD2457836.21440"_TT, "JD2457837.72526"_TT,
      "JD2457839.23669"_TT, "JD2457917.80060"_TT, "JD2457923.84629"_TT,
      "JD2457935.93288"_TT, "JD2457952.55450"_TT, "JD2457955.57554"_TT,
      "JD2457967.66254"_TT, "JD2457973.70596"_TT}},
    {"Trappist-1c",
     {"JD2457333.66400"_TT, "JD2457362.72605"_TT, "JD2457367.57051"_TT,
      "JD2457384.52320"_TT, "JD2457452.33470"_TT, "JD2457454.75672"_TT,
      "JD2457512.88094"_TT, "JD2457546.78587"_TT, "JD2457551.62888"_TT,
      "JD2457580.69137"_TT, "JD2457585.53577"_TT, "JD2457587.95622"_TT,
      "JD2457600.06684"_TT, "JD2457604.90975"_TT, "JD2457609.75461"_TT,
      "JD2457614.59710"_TT, "JD2457626.70610"_TT, "JD2457631.55024"_TT,
      "JD2457638.81518"_TT, "JD2457650.92395"_TT, "JD2457653.34553"_TT,
      "JD2457655.76785"_TT, "JD2457658.18963"_TT, "JD2457660.61168"_TT,
      "JD2457663.03292"_TT, "JD2457665.45519"_TT, "JD2457667.87729"_TT,
      "JD2457670.29869"_TT, "JD2457672.71944"_TT, "JD2457711.46778"_TT,
      "JD2457723.57663"_TT, "JD2457740.53361"_TT, "JD2457742.95276"_TT,
      "JD2457745.37429"_TT, "JD2457747.79699"_TT, "JD2457750.21773"_TT,
      "JD2457752.64166"_TT, "JD2457755.05877"_TT, "JD2457757.48313"_TT,
      "JD2457759.90281"_TT, "JD2457762.32806"_TT, "JD2457764.74831"_TT,
      "JD2457767.16994"_TT, "JD2457769.59209"_TT, "JD2457772.01483"_TT,
      "JD2457774.43458"_TT, "JD2457776.85815"_TT, "JD2457779.27911"_TT,
      "JD2457781.70095"_TT, "JD2457784.12338"_TT, "JD2457791.38801"_TT,
      "JD2457793.81141"_TT, "JD2457796.23153"_TT, "JD2457798.65366"_TT,
      "JD2457801.07631"_TT, "JD2457803.49747"_TT, "JD2457805.91882"_TT,
      "JD2457808.34123"_TT, "JD2457810.76273"_TT, "JD2457813.18456"_TT,
      "JD2457815.60583"_TT, "JD2457818.02821"_TT, "JD2457820.45019"_TT,
      "JD2457822.87188"_TT, "JD2457825.29388"_TT, "JD2457827.71513"_TT,
      "JD2457830.13713"_TT, "JD2457832.55888"_TT, "JD2457834.98120"_TT,
      "JD2457837.40280"_TT, "JD2457839.82415"_TT}},
    {"Trappist-1d",
     {"JD2457625.59779"_TT, "JD2457641.79360"_TT, "JD2457645.84360"_TT,
      "JD2457653.94261"_TT, "JD2457657.99220"_TT, "JD2457662.04284"_TT,
      "JD2457666.09140"_TT, "JD2457670.14198"_TT, "JD2457726.83975"_TT,
      "JD2457738.99169"_TT, "JD2457743.03953"_TT, "JD2457747.08985"_TT,
      "JD2457751.14022"_TT, "JD2457755.18894"_TT, "JD2457759.24638"_TT,
      "JD2457763.28895"_TT, "JD2457767.33866"_TT, "JD2457771.39077"_TT,
      "JD2457775.44026"_TT, "JD2457779.48843"_TT, "JD2457783.54023"_TT,
      "JD2457791.64083"_TT, "JD2457803.79083"_TT, "JD2457807.84032"_TT,
      "JD2457811.89116"_TT, "JD2457815.94064"_TT, "JD2457819.99050"_TT,
      "JD2457824.04185"_TT, "JD2457828.09082"_TT, "JD2457832.14036"_TT,
      "JD2457836.19171"_TT, "JD2457961.73760"_TT, "JD2457969.83708"_TT,
      "JD2457973.88590"_TT}},
    {"Trappist-1e",
     {"JD2457312.71300"_TT, "JD2457367.59683"_TT, "JD2457611.57620"_TT,
      "JD2457623.77950"_TT, "JD2457654.27862"_TT, "JD2457660.38016"_TT,
      "JD2457666.48030"_TT, "JD2457672.57930"_TT, "JD2457721.37514"_TT,
      "JD2457733.57300"_TT, "JD2457739.67085"_TT, "JD2457745.77160"_TT,
      "JD2457751.87007"_TT, "JD2457757.96712"_TT, "JD2457764.06700"_TT,
      "JD2457770.17109"_TT, "JD2457776.26378"_TT, "JD2457782.36226"_TT,
      "JD2457794.56159"_TT, "JD2457800.66354"_TT, "JD2457806.75758"_TT,
      "JD2457812.85701"_TT, "JD2457818.95510"_TT, "JD2457825.05308"_TT,
      "JD2457831.15206"_TT, "JD2457837.24980"_TT, "JD2457934.83095"_TT,
      "JD2457940.92995"_TT}},
    {"Trappist-1f",
     {"JD2457321.52520"_TT, "JD2457367.57629"_TT, "JD2457634.57809"_TT,
      "JD2457652.98579"_TT, "JD2457662.18747"_TT, "JD2457671.39279"_TT,
      "JD2457717.41541"_TT, "JD2457726.61960"_TT, "JD2457745.03116"_TT,
      "JD2457754.23380"_TT, "JD2457763.44338"_TT, "JD2457772.64752"_TT,
      "JD2457781.85142"_TT, "JD2457800.27307"_TT, "JD2457809.47554"_TT,
      "JD2457818.68271"_TT, "JD2457827.88669"_TT, "JD2457837.10322"_TT,
      "JD2457956.80549"_TT}},
    {"Trappist-1g",
     {"JD2457294.78600"_TT, "JD2457356.53410"_TT, "JD2457615.92400"_TT,
      "JD2457640.63730"_TT, "JD2457652.99481"_TT, "JD2457665.35151"_TT,
      "JD2457739.48441"_TT, "JD2457751.83993"_TT, "JD2457764.19098"_TT,
      "JD2457776.54900"_TT, "JD2457801.25000"_TT, "JD2457813.60684"_TT,
      "JD2457825.96112"_TT, "JD2457838.30655"_TT, "JD2457924.77090"_TT,
      "JD2457961.82621"_TT}},
    {"Trappist-1h",
     {"JD2457662.55467"_TT, "JD2457756.38740"_TT, "JD2457775.15390"_TT,
      "JD2457793.92300"_TT, "JD2457812.69870"_TT, "JD2457831.46625"_TT,
      "JD2457962.86271"_TT}}};

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

  TransitsByPlanet computations;

  OFStream file(TEMP_DIR / "trappist_transits.generated.wl");
  auto const bodies = ephemeris_->bodies();
  for (auto const& planet : bodies) {
    if (planet != star) {
      auto const& planet_trajectory = ephemeris_->trajectory(planet);

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
          computations[planet->name()].push_back(transit);
        }
        last_t = t;
        last_xy_displacement_derivative_sign =
            xy_displacement_derivative_sign;
      }

      file << mathematica::Assign("transit" + SanitizedName(*planet),
                                  computations[planet->name()]);
    }
  }

  LOG(ERROR)<<MaxError(observations, computations);
}

}  // namespace astronomy
}  // namespace principia
