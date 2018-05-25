
#include <random>

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
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::RelativeDegreesOfFreedom;
using physics::SolarSystem;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;

namespace astronomy {

using Transits = std::vector<Instant>;
using TransitsByPlanet = std::map<std::string, Transits>;

std::mt19937_64 engine;

class Genome {
 public:
  explicit Genome(std::vector<KeplerianElements<Trappist>> const& elements);

  void Mutate();

  static Genome OnePointCrossover(Genome const& g1, Genome const& g2);

 private:
  std::vector<KeplerianElements<Trappist>> elements_;
};

class Population {
public:
  Population(Genome const& luca, int const size);

  Genome const& Pick() const;

  void ComputeAllFitnesses(std::function<double(Genome const&)> compute_fitness);

  Population BegetChildren() const;

 private:
  std::vector<Genome> individuals_;
  std::vector<double> fitnesses_;
  std::vector<double> cumulative_fitnesses_;
};

Genome::Genome(std::vector<KeplerianElements<Trappist>> const& elements)
    : elements_(elements) {}


void Genome::Mutate()  {
  for (auto& element : elements_) {
    element.asymptotic_true_anomaly = std::nullopt;
    element.turning_angle = std::nullopt;
    element.semimajor_axis = std::nullopt;
    element.specific_energy = std::nullopt;
    element.characteristic_energy = std::nullopt;
    element.mean_motion = std::nullopt;
    element.hyperbolic_mean_motion = std::nullopt;
    element.hyperbolic_excess_velocity = std::nullopt;
    element.semiminor_axis = std::nullopt;
    element.impact_parameter = std::nullopt;
    element.semilatus_rectum = std::nullopt;
    element.specific_angular_momentum = std::nullopt;
    element.periapsis_distance = std::nullopt;
    element.apoapsis_distance = std::nullopt;
    element.longitude_of_periapsis = std::nullopt;
    element.true_anomaly = std::nullopt;
    element.hyperbolic_mean_anomaly = std::nullopt;
    std::uniform_real_distribution<> angle_distribution(-1.0, 1.0);
    *element.argument_of_periapsis += angle_distribution(engine) * Degree;
    *element.mean_anomaly += angle_distribution(engine) * Degree;
  }
}

Genome Genome::OnePointCrossover(Genome const& g1, Genome const& g2) {
  CHECK_EQ(g1.elements_.size(), g2.elements_.size());
  std::vector<KeplerianElements<Trappist>> new_elements;
  std::uniform_int_distribution<> order_distribution(0, 1);
  std::uniform_int_distribution<> split_distribution(1,
                                                     g1.elements_.size() - 1);
  bool const reverse = order_distribution(engine) == 1;
  int const split = split_distribution(engine);
  if (reverse) {
    for (int i = 0; i < split; ++i) {
      new_elements.push_back(g1.elements_[i]);
    }
    for (int i = split; i < g1.elements_.size(); ++i) {
      new_elements.push_back(g2.elements_[i]);
    }
  } else {
    for (int i = 0; i < split; ++i) {
      new_elements.push_back(g2.elements_[i]);
    }
    for (int i = split; i < g1.elements_.size(); ++i) {
      new_elements.push_back(g1.elements_[i]);
    }
  }
  return Genome(new_elements);
}

Population::Population(Genome const& luca, int const size)
    : individuals_(size, luca) {
  for (auto& individual : individuals_) {
    individual.Mutate();
  }
}

Genome const& Population::Pick() const {
  std::uniform_real_distribution<> fitness_distribution(
      cumulative_fitnesses_.front(), cumulative_fitnesses_.back());
  double const picked_fitness = fitness_distribution(engine);
  auto const picked_it = std::lower_bound(cumulative_fitnesses_.begin(),
                                          cumulative_fitnesses_.end(),
                                          picked_fitness);
  CHECK(picked_it != cumulative_fitnesses_.begin());
  CHECK(picked_it != cumulative_fitnesses_.end());
  return individuals_[std::distance(cumulative_fitnesses_.begin(),
                                    picked_it) - 1];
}

void Population::ComputeAllFitnesses(
    std::function<double(Genome const&)> compute_fitness) {
  fitnesses_.clear();
  cumulative_fitnesses_.clear();
  cumulative_fitnesses_.push_back(0.0);
  for (int i = 0; i < individuals_.size(); ++i) {
    fitnesses_.push_back(compute_fitness(individuals_[i]));
    if (i > 0) {
      cumulative_fitnesses_.push_back(cumulative_fitnesses_[i - 1] +
                                      fitnesses_.back());
    }
  }
}

Population Population::BegetChildren() const {
  return Population();
}

// TODO(phl): Literals are broken in 15.8.0 Preview 1.0 and are off by an
// integral number of days.  Use this function as a stopgap measure and switch
// to literals once MSFT have fixed their bugs.
constexpr Instant JD(double const jd) {
  return Instant{} + (jd - 2451545.0) * Day;
}

TransitsByPlanet const observations = {
    {"Trappist-1b",
     {JD(2457322.51531), JD(2457325.53910), JD(2457328.55860),
      JD(2457331.58160), JD(2457334.60480), JD(2457337.62644),
      JD(2457340.64820), JD(2457345.18028), JD(2457361.79945),
      JD(2457364.82173), JD(2457440.36492), JD(2457452.45228),
      JD(2457463.02847), JD(2457509.86460), JD(2457512.88731),
      JD(2457568.78880), JD(2457586.91824), JD(2457589.93922),
      JD(2457599.00640), JD(2457602.02805), JD(2457612.60595),
      JD(2457615.62710), JD(2457624.69094), JD(2457645.84400),
      JD(2457651.88743), JD(2457653.39809), JD(2457654.90908),
      JD(2457656.41900), JD(2457657.93129), JD(2457659.44144),
      JD(2457660.95205), JD(2457662.46358), JD(2457663.97492),
      JD(2457665.48509), JD(2457666.99567), JD(2457668.50668),
      JD(2457670.01766), JD(2457671.52876), JD(2457721.38747),
      JD(2457739.51770), JD(2457741.02787), JD(2457742.53918),
      JD(2457744.05089), JD(2457745.56164), JD(2457747.07208),
      JD(2457748.58446), JD(2457750.09387), JD(2457751.60535),
      JD(2457753.11623), JD(2457754.62804), JD(2457756.13856),
      JD(2457757.64840), JD(2457759.15953), JD(2457760.67112),
      JD(2457762.18120), JD(2457763.69221), JD(2457765.20298),
      JD(2457766.71479), JD(2457768.22514), JD(2457769.73704),
      JD(2457771.24778), JD(2457772.75738), JD(2457774.26841),
      JD(2457775.77995), JD(2457777.28899), JD(2457778.80118),
      JD(2457780.31297), JD(2457781.82231), JD(2457783.33410),
      JD(2457784.84372), JD(2457792.39979), JD(2457793.90955),
      JD(2457795.41987), JD(2457796.93134), JD(2457798.44211),
      JD(2457799.95320), JD(2457801.46314), JD(2457802.97557),
      JD(2457804.48638), JD(2457805.99697), JD(2457807.50731),
      JD(2457809.01822), JD(2457810.52781), JD(2457812.04038),
      JD(2457813.55121), JD(2457815.06275), JD(2457816.57335),
      JD(2457818.08382), JD(2457819.59478), JD(2457821.10550),
      JD(2457824.12730), JD(2457825.63813), JD(2457827.14995),
      JD(2457828.66042), JD(2457830.17087), JD(2457833.19257),
      JD(2457834.70398), JD(2457836.21440), JD(2457837.72526),
      JD(2457839.23669), JD(2457917.80060), JD(2457923.84629),
      JD(2457935.93288), JD(2457952.55450), JD(2457955.57554),
      JD(2457967.66254), JD(2457973.70596)}},
    {"Trappist-1c",
     {JD(2457333.66400), JD(2457362.72605), JD(2457367.57051),
      JD(2457384.52320), JD(2457452.33470), JD(2457454.75672),
      JD(2457512.88094), JD(2457546.78587), JD(2457551.62888),
      JD(2457580.69137), JD(2457585.53577), JD(2457587.95622),
      JD(2457600.06684), JD(2457604.90975), JD(2457609.75461),
      JD(2457614.59710), JD(2457626.70610), JD(2457631.55024),
      JD(2457638.81518), JD(2457650.92395), JD(2457653.34553),
      JD(2457655.76785), JD(2457658.18963), JD(2457660.61168),
      JD(2457663.03292), JD(2457665.45519), JD(2457667.87729),
      JD(2457670.29869), JD(2457672.71944), JD(2457711.46778),
      JD(2457723.57663), JD(2457740.53361), JD(2457742.95276),
      JD(2457745.37429), JD(2457747.79699), JD(2457750.21773),
      JD(2457752.64166), JD(2457755.05877), JD(2457757.48313),
      JD(2457759.90281), JD(2457762.32806), JD(2457764.74831),
      JD(2457767.16994), JD(2457769.59209), JD(2457772.01483),
      JD(2457774.43458), JD(2457776.85815), JD(2457779.27911),
      JD(2457781.70095), JD(2457784.12338), JD(2457791.38801),
      JD(2457793.81141), JD(2457796.23153), JD(2457798.65366),
      JD(2457801.07631), JD(2457803.49747), JD(2457805.91882),
      JD(2457808.34123), JD(2457810.76273), JD(2457813.18456),
      JD(2457815.60583), JD(2457818.02821), JD(2457820.45019),
      JD(2457822.87188), JD(2457825.29388), JD(2457827.71513),
      JD(2457830.13713), JD(2457832.55888), JD(2457834.98120),
      JD(2457837.40280), JD(2457839.82415)}},
    {"Trappist-1d",
     {JD(2457625.59779), JD(2457641.79360), JD(2457645.84360),
      JD(2457653.94261), JD(2457657.99220), JD(2457662.04284),
      JD(2457666.09140), JD(2457670.14198), JD(2457726.83975),
      JD(2457738.99169), JD(2457743.03953), JD(2457747.08985),
      JD(2457751.14022), JD(2457755.18894), JD(2457759.24638),
      JD(2457763.28895), JD(2457767.33866), JD(2457771.39077),
      JD(2457775.44026), JD(2457779.48843), JD(2457783.54023),
      JD(2457791.64083), JD(2457803.79083), JD(2457807.84032),
      JD(2457811.89116), JD(2457815.94064), JD(2457819.99050),
      JD(2457824.04185), JD(2457828.09082), JD(2457832.14036),
      JD(2457836.19171), JD(2457961.73760), JD(2457969.83708),
      JD(2457973.88590)}},
    {"Trappist-1e",
     {JD(2457312.71300), JD(2457367.59683), JD(2457611.57620),
      JD(2457623.77950), JD(2457654.27862), JD(2457660.38016),
      JD(2457666.48030), JD(2457672.57930), JD(2457721.37514),
      JD(2457733.57300), JD(2457739.67085), JD(2457745.77160),
      JD(2457751.87007), JD(2457757.96712), JD(2457764.06700),
      JD(2457770.17109), JD(2457776.26378), JD(2457782.36226),
      JD(2457794.56159), JD(2457800.66354), JD(2457806.75758),
      JD(2457812.85701), JD(2457818.95510), JD(2457825.05308),
      JD(2457831.15206), JD(2457837.24980), JD(2457934.83095),
      JD(2457940.92995)}},
    {"Trappist-1f",
     {JD(2457321.52520), JD(2457367.57629), JD(2457634.57809),
      JD(2457652.98579), JD(2457662.18747), JD(2457671.39279),
      JD(2457717.41541), JD(2457726.61960), JD(2457745.03116),
      JD(2457754.23380), JD(2457763.44338), JD(2457772.64752),
      JD(2457781.85142), JD(2457800.27307), JD(2457809.47554),
      JD(2457818.68271), JD(2457827.88669), JD(2457837.10322),
      JD(2457956.80549)}},
    {"Trappist-1g",
     {JD(2457294.78600), JD(2457356.53410), JD(2457615.92400),
      JD(2457640.63730), JD(2457652.99481), JD(2457665.35151),
      JD(2457739.48441), JD(2457751.83993), JD(2457764.19098),
      JD(2457776.54900), JD(2457801.25000), JD(2457813.60684),
      JD(2457825.96112), JD(2457838.30655), JD(2457924.77090),
      JD(2457961.82621)}},
    {"Trappist-1h",
     {JD(2457662.55467), JD(2457756.38740), JD(2457775.15390),
      JD(2457793.92300), JD(2457812.69870), JD(2457831.46625),
      JD(2457962.86271)}}};

class TrappistDynamicsTest : public ::testing::Test {
 protected:
  TrappistDynamicsTest()
      : system_(SOLUTION_DIR / "astronomy" / "trappist_gravity_model.proto.txt",
                SOLUTION_DIR / "astronomy" /
                    "trappist_initial_state_jd_2457000_000000000.proto.txt"),
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
        if (error > max_error) {
          max_error = error;
          LOG(ERROR) << name << ": "
                     << (*std::prev(next_computed_transit) - JD(2450000.0)) /
                            Day
                     << " < " << (observed_transit - JD(2450000.0)) / Day
                     << " < "
                     << (*next_computed_transit - JD(2450000.0)) / Day;
        }
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

  LOG(ERROR) << "max error: " << MaxError(observations, computations);
}

}  // namespace astronomy
}  // namespace principia
