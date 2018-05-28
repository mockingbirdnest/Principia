
#include <random>

#include "astronomy/frames.hpp"
#include "base/bundle.hpp"
#include "base/file.hpp"
#include "base/not_null.hpp"
#include "base/status.hpp"
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
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using base::Bundle;
using base::not_null;
using base::OFStream;
using base::Status;
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
using quantities::Square;
using quantities::Sqrt;
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

// The random number generator used by the optimisation.

// The description of the characteristics of an individual, i.e., a
// configuration of the Trappist system.
class Genome {
 public:
  explicit Genome(std::vector<KeplerianElements<Trappist>> const& elements);

  std::vector<KeplerianElements<Trappist>> const& elements() const;

  void Mutate(std::mt19937_64& engine, int generation);

  static Genome OnePointCrossover(Genome const& g1,
                                  Genome const& g2,
                                  std::mt19937_64& engine);
  static Genome TwoPointCrossover(Genome const& g1,
                                  Genome const& g2,
                                  std::mt19937_64& engine);
  static Genome Blend(Genome const& g1,
                      Genome const& g2,
                      std::mt19937_64& engine);

 private:
  std::vector<KeplerianElements<Trappist>> elements_;
};

// A set of genomes which can reproduce based on their fitness.
class Population {
 public:
  Population(Genome const& luca,
             int const size,
             std::function<double(Genome const&)> compute_fitness,
             std::function<std::string(Genome const&)> residual_trace);

  void ComputeAllFitnesses();

  void BegetChildren();

  Genome best_genome() const;

 private:
  Genome const* Pick() const;

  std::function<double(Genome const&)> const compute_fitness_;
  std::function<std::string(Genome const&)> const residual_trace_;
  mutable std::mt19937_64 engine_;
  std::vector<Genome> current_;
  std::vector<Genome> next_;
  std::vector<double> fitnesses_;
  std::vector<double> cumulative_fitnesses_;

  int generation_ = 0;
  double best_fitness_ = 0.0;
  std::optional<Genome> best_genome_;
};

Genome::Genome(std::vector<KeplerianElements<Trappist>> const& elements)
    : elements_(elements) {}

std::vector<KeplerianElements<Trappist>> const& Genome::elements() const {
  return elements_;
}

void Genome::Mutate(std::mt19937_64& engine, int generation)  {
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
    // The standard deviation of the distribution below has a strong effect on
    // the convergence of the algorithm: if it's too small we do not explore the
    // genomic space efficiently and it takes forever to find decent solutions;
    // if it's too large we explore the genomic space haphazardly and suffer
    // from deleterious mutations.
    double multiplicator =
        std::exp2(-2 - (generation % 100) / 15 - 1 * (generation / 100));
    if (generation == -1) multiplicator = 1;
    std::student_t_distribution<> distribution(1);
    *element.argument_of_periapsis +=
        distribution(engine) * 10 * Degree * multiplicator;
    *element.argument_of_periapsis =
        std::fmod(*element.argument_of_periapsis / quantities::si::Radian,
                  2 * π) *
        quantities::si::Radian;
    if (*element.argument_of_periapsis < 0 * quantities::si::Radian) {
      *element.argument_of_periapsis += 2 * π * quantities::si::Radian;
    }
    *element.mean_anomaly += distribution(engine) * 10 * Degree * multiplicator;
    *element.mean_anomaly =
        std::fmod(*element.mean_anomaly / quantities::si::Radian, 2 * π) *
        quantities::si::Radian;
    if (*element.mean_anomaly < 0 * quantities::si::Radian) {
      *element.mean_anomaly += 2 * π * quantities::si::Radian;
    }
    *element.period *= 1 + distribution(engine) * 1e-6 * multiplicator;
    *element.eccentricity += distribution(engine) * 1e-3 * multiplicator;
  }
}

Genome Genome::OnePointCrossover(Genome const& g1,
                                 Genome const& g2,
                                 std::mt19937_64& engine) {
  CHECK_EQ(g1.elements_.size(), g2.elements_.size());
  std::vector<KeplerianElements<Trappist>> new_elements;
  std::uniform_int_distribution<> order_distribution(0, 1);
  std::uniform_int_distribution<> split_distribution(0, g1.elements_.size());
  bool const reverse = order_distribution(engine) == 1;
  int const split = split_distribution(engine);
  if (reverse) {
    for (int i = 0; i < split; ++i) {
      new_elements.push_back(g1.elements_[i]);
    }
    for (int i = split; i < g2.elements_.size(); ++i) {
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

Genome Genome::TwoPointCrossover(Genome const& g1,
                                 Genome const& g2,
                                 std::mt19937_64& engine) {
  CHECK_EQ(g1.elements_.size(), g2.elements_.size());
  std::vector<KeplerianElements<Trappist>> new_elements;
  std::uniform_int_distribution<> order_distribution(0, 1);
  std::uniform_int_distribution<> split_distribution(0, g1.elements_.size());
  bool const reverse = order_distribution(engine) == 1;
  int split1 = split_distribution(engine);
  int split2 = split_distribution(engine);
  if (split2 < split1) {
    std::swap(split1, split2);
  }
  if (reverse) {
    for (int i = 0; i < split1; ++i) {
      new_elements.push_back(g1.elements_[i]);
    }
    for (int i = split1; i < split2; ++i) {
      new_elements.push_back(g2.elements_[i]);
    }
    for (int i = split2; i < g1.elements_.size(); ++i) {
      new_elements.push_back(g1.elements_[i]);
    }
  } else {
    for (int i = 0; i < split1; ++i) {
      new_elements.push_back(g2.elements_[i]);
    }
    for (int i = split1; i < split2; ++i) {
      new_elements.push_back(g1.elements_[i]);
    }
    for (int i = split2; i < g2.elements_.size(); ++i) {
      new_elements.push_back(g2.elements_[i]);
    }
  }
  return Genome(new_elements);
}

Genome Genome::Blend(Genome const& g1,
                     Genome const& g2,
                     std::mt19937_64& engine) {
  CHECK_EQ(g1.elements_.size(), g2.elements_.size());
  std::vector<KeplerianElements<Trappist>> new_elements;
  std::uniform_real_distribution blend_distribution(0.0, 1.0);
  double const blend = blend_distribution(engine);
  for (int i = 0; i < g1.elements_.size(); ++i) {
    KeplerianElements<Trappist> new_element = g1.elements_[i];
    *new_element.argument_of_periapsis =
        *g1.elements_[i].argument_of_periapsis * blend +
        *g2.elements_[i].argument_of_periapsis * (1.0 - blend);
    *new_element.argument_of_periapsis =
        *g1.elements_[i].mean_anomaly * blend +
        *g2.elements_[i].mean_anomaly * (1.0 - blend);
    new_elements.push_back(new_element);
  }
  return Genome(new_elements);
}

Population::Population(Genome const& luca,
                       int const size,
                       std::function<double(Genome const&)> compute_fitness,
                       std::function<std::string(Genome const&)> residual_trace)
    : current_(size, luca),
      next_(size, luca),
      compute_fitness_(std::move(compute_fitness)),
      residual_trace_(std::move(residual_trace)) {
  for (int i = 0; i < current_.size(); ++i) {
    current_[i].Mutate(engine_, -1);
  }
}

void Population::ComputeAllFitnesses() {
  // The fitness computation is expensive, do it in parallel on all genomes.
  {
    Bundle bundle(4);

    fitnesses_.resize(current_.size(), 0.0);
    for (int i = 0; i < current_.size(); ++i) {
      bundle.Add([this, i]() {
        fitnesses_[i] = compute_fitness_(current_[i]);
        return Status();
      });
    }
    bundle.Join();
  }
  
  LOG(ERROR) << "Generation " << generation_;
  double min_fitness = std::numeric_limits<double>::infinity();
  double max_fitness = 0.0;
  Genome const* fittest = nullptr;
  Genome const* least_fit = nullptr;
  cumulative_fitnesses_.clear();
  cumulative_fitnesses_.push_back(0.0);
  for (int i = 0; i < current_.size(); ++i) {
    double const fitness = fitnesses_[i];
    cumulative_fitnesses_.push_back(cumulative_fitnesses_[i] + fitness);
    if (fitness > max_fitness) {
      max_fitness = fitness;
      fittest = &current_[i];
    }
    if (fitness < min_fitness) {
      min_fitness = fitness;
      least_fit = &current_[i];
    }
    if (fitness > best_fitness_) {
      best_fitness_ = fitness;
      LOG(ERROR) << "New best genome:";
      char planet = 'b';
      for (int j = 0; j < current_[i].elements().size(); ++j) {
        LOG(ERROR) << std::string({planet++, ':'});
        if (best_genome_) {
          LOG(ERROR)
              << "old L = "
              << (best_genome_->elements()[j].longitude_of_ascending_node +
                  *best_genome_->elements()[j].argument_of_periapsis +
                  *best_genome_->elements()[j].mean_anomaly) /
                     Degree
              << u8"°";
          LOG(ERROR)
              << u8"   ΔL = "
              << ((current_[i].elements()[j].longitude_of_ascending_node +
                   *current_[i].elements()[j].argument_of_periapsis +
                   *current_[i].elements()[j].mean_anomaly) -
                  (best_genome_->elements()[j].longitude_of_ascending_node +
                   *best_genome_->elements()[j].argument_of_periapsis +
                   *best_genome_->elements()[j].mean_anomaly)) /
                     Degree
              << u8"°";
        }
        LOG(ERROR) << "new L = "
                   << (current_[i].elements()[j].longitude_of_ascending_node +
                       *current_[i].elements()[j].argument_of_periapsis +
                       *current_[i].elements()[j].mean_anomaly) /
                          Degree
                   << u8"°";
        if (best_genome_) {
          LOG(ERROR) << "old e = " << *best_genome_->elements()[j].eccentricity;
          LOG(ERROR) << u8"   Δe = "
                     << *current_[i].elements()[j].eccentricity -
                            *best_genome_->elements()[j].eccentricity;
        }
        LOG(ERROR) << "new e = " << *current_[i].elements()[j].eccentricity;
      }
      best_genome_ = current_[i];
    }
  }
  LOG(ERROR) << "Min: " << min_fitness << " Max: " << max_fitness
             << " Best: " << best_fitness_;
  LOG(ERROR) << "Least fit: " << residual_trace_(*least_fit);
  LOG(ERROR) << "Fittest  : " << residual_trace_(*fittest);
  LOG(ERROR) << "Best     : " << residual_trace_(*best_genome_);
}

void Population::BegetChildren() {
  for (int i = 0; i < next_.size(); ++i) {
    Genome const* const parent1 = Pick();
    Genome const* parent2;
    // Let's have sex like snails: if we find a good partner, fine, otherwise
    // let's go for self-fecundation.
    for (int j = 0; j < 2; ++j) {
      parent2 = Pick();
      if (parent1 != parent2) {
        break;
      }
    }
    next_[i] = Genome::TwoPointCrossover(*parent1, *parent2, engine_);
    next_[i].Mutate(engine_, generation_);
  }
  next_.swap(current_);
  ++generation_;
}

Genome Population::best_genome() const {
  return *best_genome_;
}

Genome const* Population::Pick() const {
  std::uniform_real_distribution<> fitness_distribution(
      cumulative_fitnesses_.front(), cumulative_fitnesses_.back());
  double const picked_fitness = fitness_distribution(engine_);
  auto const picked_it = std::lower_bound(cumulative_fitnesses_.begin(),
                                          cumulative_fitnesses_.end(),
                                          picked_fitness);
  CHECK(picked_it != cumulative_fitnesses_.begin());
  CHECK(picked_it != cumulative_fitnesses_.end());
  int const picked_index =
      std::distance(cumulative_fitnesses_.begin(), picked_it) - 1;
  CHECK_LE(0, picked_index);
  CHECK_LT(picked_index, current_.size());
  return &current_[picked_index];
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
                    "trappist_initial_state_jd_2457010_000000000.proto.txt"),
        ephemeris_(system_.MakeEphemeris(
            /*fitting_tolerance=*/5 * Milli(Metre),
            Ephemeris<Trappist>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                                   Position<Trappist>>(),
                /*step=*/0.07 * Day))) {}

  static Transits ComputeTransits(Ephemeris<Trappist> const& ephemeris,
                                  not_null<MassiveBody const*> const star,
                                  not_null<MassiveBody const*> const planet) {
    Transits transits;
    auto const& star_trajectory = ephemeris.trajectory(star);

    std::optional<Instant> last_t;
    std::optional<Sign> last_xy_displacement_derivative_sign;
      auto const& planet_trajectory = ephemeris.trajectory(planet);
    for (Instant t = ephemeris.t_min();
          t < ephemeris.t_max();
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
        transits.push_back(transit);
      }
      last_t = t;
      last_xy_displacement_derivative_sign =
          xy_displacement_derivative_sign;
    }
    return transits;
  }

  static double ShortDays(Instant const& time) {
    return (time - JD(2450000.0)) / Day;
  }

  static double χ²(TransitsByPlanet const& observations,
                    TransitsByPlanet const& computations,
                    std::vector<Time>* errors) {
    int number_of_transits = 0;
    double sum_of_squared_errors = 0;
    if (errors != nullptr) {
      errors->clear();
    }
    for (auto const& pair : observations) {
      auto const& name = pair.first;
      auto const& observed_transits = pair.second;
      auto const& computed_transits = computations.at(name);
      if (computed_transits.empty()) {
        return std::numeric_limits<double>::infinity();
      }
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
        if (errors != nullptr) {
          errors->push_back(error);
        }
        Time const σ = 0.001 * Day;
        sum_of_squared_errors += quantities::Pow<2>(error / σ);
      }
      number_of_transits += observed_transits.size();
    }
    return sum_of_squared_errors - 1;
  }

  static std::string SanitizedName(MassiveBody const& body) {
    auto sanitized_name = body.name();
    return sanitized_name.erase(sanitized_name.find_first_of("-"), 1);
  }

  constexpr static char star_name[] = "Trappist-1A";
  SolarSystem<Trappist> const system_;
  not_null<std::unique_ptr<Ephemeris<Trappist>>> ephemeris_;
};

constexpr char TrappistDynamicsTest::star_name[];

TEST_F(TrappistDynamicsTest, MathematicaPeriods) {
  Instant const a_century_later = system_.epoch() + 100 * JulianYear;
  ephemeris_->Prolong(a_century_later);

  auto const& star = system_.massive_body(*ephemeris_, star_name);
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

  TransitsByPlanet computations;
  OFStream file(TEMP_DIR / "trappist_transits.generated.wl");

  auto const& star = system_.massive_body(*ephemeris_, star_name);
  auto const bodies = ephemeris_->bodies();
  for (auto const& planet : bodies) {
    if (planet != star) {
      computations[planet->name()] = ComputeTransits(*ephemeris_, star, planet);
      file << mathematica::Assign("transit" + SanitizedName(*planet),
                                  computations[planet->name()]);
    }
  }

  LOG(ERROR) << "min probability: "
             << std::exp(χ²(
                    observations, computations, nullptr));
}

TEST_F(TrappistDynamicsTest, Optimisation) {
  SolarSystem<Trappist> const system(
      SOLUTION_DIR / "astronomy" / "trappist_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "trappist_initial_state_jd_2457010_000000000.proto.txt");

  auto planet_names = system.names();
  planet_names.erase(
      std::find(planet_names.begin(), planet_names.end(), star_name));
  std::vector<KeplerianElements<Trappist>> elements;
  for (auto const& planet_name : planet_names) {
    elements.push_back(SolarSystem<Trappist>::MakeKeplerianElements(
        system.keplerian_initial_state_message(planet_name).elements()));
  }

  auto residual_trace = [&planet_names, &system](Genome const& genome) {
    auto modified_system = system;
    auto const& elements = genome.elements();
    for (int i = 0; i < planet_names.size(); ++i) {
      modified_system.ReplaceElements(planet_names[i], elements[i]);
    }

    auto const ephemeris = modified_system.MakeEphemeris(
            /*fitting_tolerance=*/5 * Milli(Metre),
            Ephemeris<Trappist>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                                   Position<Trappist>>(),
                /*step=*/0.07 * Day));
    ephemeris->Prolong(modified_system.epoch() + 1000 * Day);
    if (!ephemeris->last_severe_integration_status().ok()) {
      return std::string(u8"ἀποκάλυψις");
    }

    TransitsByPlanet computations;
    auto const& star = modified_system.massive_body(*ephemeris, star_name);
    auto const bodies = ephemeris->bodies();
    for (auto const& planet : bodies) {
      if (planet != star) {
        computations[planet->name()] =
            ComputeTransits(*ephemeris, star, planet);
      }
    }

    std::vector<Time> δt;
    double const χ²_statistic = χ²(observations, computations, &δt);
    Time max = 0 * Second;
    Time total = 0 * Second;
    for (auto const residual : δt) {
      max = std::max(max, residual);
      total += residual;
    }
    return u8"χ² = " + std::to_string(χ²_statistic) +
           "; n = " + std::to_string(δt.size()) + u8"; max Δt = " +
           std::to_string(max / Second) + u8" s; avg Δt = " +
           std::to_string(total / δt.size() / Second) + " s";
  };

  auto compute_fitness = [&planet_names, &system](Genome const& genome) {
    auto modified_system = system;
    auto const& elements = genome.elements();
    for (int i = 0; i < planet_names.size(); ++i) {
      modified_system.ReplaceElements(planet_names[i], elements[i]);
    }

    auto const ephemeris = modified_system.MakeEphemeris(
            /*fitting_tolerance=*/5 * Milli(Metre),
            Ephemeris<Trappist>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                                   Position<Trappist>>(),
                /*step=*/0.07 * Day));
    ephemeris->Prolong(modified_system.epoch() + 1000 * Day);
    if (!ephemeris->last_severe_integration_status().ok()) {
      return 0.0;
    }

    TransitsByPlanet computations;
    auto const& star = modified_system.massive_body(*ephemeris, star_name);
    auto const bodies = ephemeris->bodies();
    for (auto const& planet : bodies) {
      if (planet != star) {
        computations[planet->name()] =
            ComputeTransits(*ephemeris, star, planet);
      }
    }

    // This is the place where we cook the sausage.  This function must be steep
    // enough to efficiently separate the wheat from the chaff without leading
    // to monoculture.
    return 1/χ²(observations, computations, nullptr);
  };

  Genome luca(elements);
  Population population(
      luca, 50, std::move(compute_fitness), std::move(residual_trace));
  for (int i = 0; i < 400; ++i) {
    population.ComputeAllFitnesses();
    population.BegetChildren();
  }

  for (int i = 0; i < planet_names.size(); ++i) {
    LOG(ERROR) << planet_names[i] << ": "
               << population.best_genome().elements()[i];
  }
}

}  // namespace astronomy
}  // namespace principia
