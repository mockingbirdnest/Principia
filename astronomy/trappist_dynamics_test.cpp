
#include <algorithm>
#include <iomanip>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/bundle.hpp"
#include "base/file.hpp"
#include "base/not_null.hpp"
#include "base/status.hpp"
#include "geometry/grassmann.hpp"
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
using geometry::Displacement;
using geometry::InnerProduct;
using geometry::Instant;
using geometry::OrientedAngleBetween;
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
using quantities::Abs;
using quantities::Angle;
using quantities::ArcTan;
using quantities::Cos;
using quantities::Derivative;
using quantities::Difference;
using quantities::Mod;
using quantities::Pow;
using quantities::Sin;
using quantities::Square;
using quantities::Sqrt;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;

namespace astronomy {

namespace genetics {

// The description of the characteristics of an individual, i.e., a
// configuration of the Trappist system.  This is merely a wrapper on the
// elements of each planet.
class Genome {
 public:
  explicit Genome(std::vector<KeplerianElements<Sky>> const& elements);

  std::vector<KeplerianElements<Sky>> const& elements() const;

  // Mutate this individual depending on its generation.
  void Mutate(std::mt19937_64& engine, int generation);

  // We tried other forms of crossover (one-point crossover, linear blending)
  // but then didn't produce good results.
  static Genome TwoPointCrossover(Genome const& g1,
                                  Genome const& g2,
                                  std::mt19937_64& engine);

 private:
  std::vector<KeplerianElements<Sky>> elements_;
};

// A set of genomes which can reproduce based on their fitness.
class Population {
 public:
  using ComputeFitness = std::function<double(Genome const&, std::string&)>;

  // Constructs an initial population made of |size| mutated copies of |luca|.
  // If |elitism| is true, the best individual is preserved unchanged in the
  // next generation.
  Population(Genome const& luca,
             int size,
             bool elitism,
             ComputeFitness compute_fitness,
             std::mt19937_64& engine);

  // Compute all the fitnesses for the current population, as well as the
  // cumulative fitnesses used for reproduction.  This is the expensive step.
  void ComputeAllFitnesses();

  // Produce the next generation by crossover and mutation depending on the
  // fitness of each individual.
  void BegetChildren();

  Genome best_genome() const;
  double best_genome_fitness() const;
  std::string best_genome_trace() const;

 private:
  Genome const* Pick() const;

  void TraceNewBestGenome(Genome const& genome) const;

  ComputeFitness const compute_fitness_;
  bool const elitism_;
  std::mt19937_64& engine_;
  std::vector<Genome> current_;
  std::vector<Genome> next_;
  std::vector<double> fitnesses_;
  std::vector<std::string> traces_;
  std::vector<double> cumulative_fitnesses_;

  int generation_ = 0;
  double best_fitness_ = 0.0;
  std::string best_trace_;
  std::optional<Genome> best_genome_;
};

Genome::Genome(std::vector<KeplerianElements<Sky>> const& elements)
    : elements_(elements) {}

std::vector<KeplerianElements<Sky>> const& Genome::elements() const {
  return elements_;
}

void Genome::Mutate(std::mt19937_64& engine,
                    int const generation) {
  std::student_t_distribution<> distribution(1);

  // The scale of the perturbation has a strong effect on the convergence of the
  // algorithm: if it's too small we do not explore the genomic space
  // efficiently and it takes forever to find decent solutions; if it's too
  // large we explore the genomic space haphazardly and suffer from deleterious
  // mutations.  The |multiplicator| is used to decay the perturbation over
  // time.
  double const multiplicator =
      generation == -1 ? 1 : std::exp2(-2 - std::min(generation, 800) / 120);

  for (auto& element : elements_) {
    KeplerianElements<Sky> mutated_element;
    // The elements that are preserved.
    mutated_element.inclination = element.inclination;
    mutated_element.longitude_of_ascending_node =
        element.longitude_of_ascending_node;
    // The elements that are mutated.
    mutated_element.argument_of_periapsis =
        *element.argument_of_periapsis +
        distribution(engine) * 10 * Degree * multiplicator;
    mutated_element.argument_of_periapsis =
        Mod(*mutated_element.argument_of_periapsis, 2 * π * Radian);
    mutated_element.period =
        *element.period +
        distribution(engine) * 5 * Second * Sqrt(multiplicator);
    mutated_element.eccentricity =
        *element.eccentricity +
        distribution(engine) * 1e-3 * Sqrt(multiplicator);
    for (;;) {
      if (*mutated_element.eccentricity < 0.0) {
        *mutated_element.eccentricity = -*mutated_element.eccentricity;
      } else if (*mutated_element.eccentricity > 0.2) {
        *mutated_element.eccentricity = 0.2 - *mutated_element.eccentricity;
      } else {
        break;
      }
    }
    mutated_element.mean_anomaly =
        *element.mean_anomaly +
        distribution(engine) * 10 * Degree * multiplicator;
    mutated_element.mean_anomaly =
        Mod(*mutated_element.mean_anomaly, 2 * π * Radian);

    element = mutated_element;
  }
}

Genome Genome::TwoPointCrossover(Genome const& g1,
                                 Genome const& g2,
                                 std::mt19937_64& engine) {
  CHECK_EQ(g1.elements_.size(), g2.elements_.size());
  std::vector<KeplerianElements<Sky>> new_elements;
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

Population::Population(Genome const& luca,
                       int const size,
                       bool const elitism,
                       ComputeFitness compute_fitness,
                       std::mt19937_64& engine)
    : current_(size, luca),
      next_(size, luca),
      compute_fitness_(std::move(compute_fitness)),
      elitism_(elitism),
      engine_(engine) {
  for (int i = 0; i < current_.size(); ++i) {
    current_[i].Mutate(engine_, /*generation=*/-1);
  }
}

void Population::ComputeAllFitnesses() {
  // The fitness computation is expensive, do it in parallel on all genomes.
  {
    Bundle bundle;

    fitnesses_.resize(current_.size(), 0.0);
    traces_.resize(current_.size(), "");
    int i = 0;
    if (elitism_ && best_genome_.has_value()) {
      // Elitism: fitnesses_[0] is already best_fitness_.
      CHECK_EQ(fitnesses_[i], best_fitness_);
      ++i;
    }
    for (; i < current_.size(); ++i) {
      bundle.Add([this, i]() {
        fitnesses_[i] = compute_fitness_(current_[i], traces_[i]);
        return Status();
      });
    }
    bundle.Join();
  }

  LOG(ERROR) << "------ Generation " << generation_;
  double min_fitness = std::numeric_limits<double>::infinity();
  double max_fitness = 0.0;
  std::string* fittest_info = nullptr;
  std::string* least_fit_info = nullptr;

  // Compute the cumulative fitness: we'll use it to pick who is allowed to
  // reproduce.
  cumulative_fitnesses_.clear();
  cumulative_fitnesses_.push_back(0.0);
  for (int i = 0; i < current_.size(); ++i) {
    double const fitness = fitnesses_[i];
    cumulative_fitnesses_.push_back(cumulative_fitnesses_[i] + fitness);

    if (fitness > max_fitness) {
      max_fitness = fitness;
      fittest_info = &traces_[i];
    }
    if (fitness < min_fitness) {
      min_fitness = fitness;
      least_fit_info = &traces_[i];
    }
    if (fitness > best_fitness_) {
      best_fitness_ = fitness;
      best_trace_ = traces_[i];
      TraceNewBestGenome(current_[i]);
      best_genome_ = current_[i];
    }
  }
  LOG(ERROR) << "Least fit: " << *least_fit_info;
  LOG(ERROR) << "Fittest  : " << *fittest_info;
  if (!elitism_) {
    LOG(ERROR) << "Best     : " << best_trace_;
  }
}

void Population::BegetChildren() {
  int i = 0;
  if (elitism_ && best_genome_.has_value()) {
    // Elitism: always retain the best genome at index 0.
    next_[i] = *best_genome_;
    fitnesses_[i] = best_fitness_;
    traces_[i] = best_trace_;
    ++i;
  }
  for (; i < next_.size(); ++i) {
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

double Population::best_genome_fitness() const {
  return best_fitness_;
}

std::string Population::best_genome_trace() const {
  return best_trace_;
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

void Population::TraceNewBestGenome(Genome const& genome) const {
  LOG(ERROR) << "New best genome:";
  char planet = 'b';
  for (int j = 0; j < genome.elements().size(); ++j) {
    LOG(ERROR) << std::string({planet++, ':'});
    if (best_genome_) {
      LOG(ERROR)
          << "old L = "
          << Mod((best_genome_->elements()[j].longitude_of_ascending_node +
                  *best_genome_->elements()[j].argument_of_periapsis +
                  *best_genome_->elements()[j].mean_anomaly),
                 2 * π * Radian) / Degree
          << u8"°";

      LOG(ERROR) << u8"   ΔL = "
                 << ((genome.elements()[j].longitude_of_ascending_node +
                      *genome.elements()[j].argument_of_periapsis +
                      *genome.elements()[j].mean_anomaly) -
                     (best_genome_->elements()[j].longitude_of_ascending_node +
                      *best_genome_->elements()[j].argument_of_periapsis +
                      *best_genome_->elements()[j].mean_anomaly)) / Degree
                 << u8"°";
    }
    LOG(ERROR) << "new L = "
               << Mod((genome.elements()[j].longitude_of_ascending_node +
                       *genome.elements()[j].argument_of_periapsis +
                       *genome.elements()[j].mean_anomaly),
                      2 * π * Radian) / Degree
               << u8"°";
    if (best_genome_) {
      LOG(ERROR) << "old e = " << *best_genome_->elements()[j].eccentricity;
      LOG(ERROR) << u8"   Δe = "
                 << *genome.elements()[j].eccentricity -
                        *best_genome_->elements()[j].eccentricity;
    }
    LOG(ERROR) << "new e = " << *genome.elements()[j].eccentricity;
    if (best_genome_) {
      LOG(ERROR) << "old T = " << *best_genome_->elements()[j].period / Day
                 << " d";
      LOG(ERROR) << u8"   ΔT = "
                 << (*genome.elements()[j].period -
                     *best_genome_->elements()[j].period) / Second
                 << " s";
    }
    LOG(ERROR) << "new T = " << *genome.elements()[j].period / Day << " d";
  }
}
}  // namespace genetics

// DEMCMC stands for Differential Evolution - Марков Chain Monte-Carlo.
// See for instance "A Markov Chain Monte Carlo version of the genetic algorithm
// Differential Evolution: easy Bayesian computing for real parameter spaces",
// Braak, 2006.
namespace deмcmc {

// These parameters follow "The Laplace resonance in the Kepler-60 planetary
// system", Goździewski et al., 2016.
struct PlanetParameters {
  constexpr static int count = 4;
  Time period;
  double x{};
  double y{};
  Time time_to_first_transit;
};

using SystemParameters = std::array<PlanetParameters, 7>;
using Population = std::vector<SystemParameters>;
using ComputeLogPdf =
    std::function<double(SystemParameters const&, std::string&)>;

PlanetParameters operator+(PlanetParameters const& left,
                           PlanetParameters const& right) {
  PlanetParameters result = left;
  result.period += right.period;
  result.x += right.x;
  result.y += right.y;
  result.time_to_first_transit += right.time_to_first_transit;
  return result;
}

PlanetParameters operator-(PlanetParameters const& left,
                           PlanetParameters const& right) {
  PlanetParameters result = left;
  result.period -= right.period;
  result.x -= right.x;
  result.y -= right.y;
  result.time_to_first_transit -= right.time_to_first_transit;
  return result;
}

PlanetParameters operator*(double const left, PlanetParameters const& right) {
  PlanetParameters result = right;
  result.period *= left;
  result.x *= left;
  result.y *= left;
  result.time_to_first_transit *= left;
  return result;
}

SystemParameters operator+(SystemParameters const& left,
                           SystemParameters const& right) {
  SystemParameters result = left;
  for (int i = 0; i < result.size(); ++i) {
    result[i] = result[i] + right[i];
  }
  return result;
}

SystemParameters operator-(SystemParameters const& left,
                           SystemParameters const& right) {
  SystemParameters result = left;
  for (int i = 0; i < result.size(); ++i) {
    result[i] = result[i] - right[i];
  }
  return result;
}

SystemParameters operator*(double const left, SystemParameters const& right) {
  SystemParameters result = right;
  for (int i = 0; i < result.size(); ++i) {
    result[i] = left * result[i];
  }
  return result;
}

KeplerianElements<Sky> MakeKeplerianElements(
    KeplerianElements<Sky> const& blueprint,
    PlanetParameters const& parameters) {
  KeplerianElements<Sky> elements;
  // The elements that are copied from the blueprint.
  elements.inclination = blueprint.inclination;
  elements.longitude_of_ascending_node = blueprint.longitude_of_ascending_node;
  // The elements that are computed from the parameters.
  elements.eccentricity = Sqrt(Pow<2>(parameters.x) + Pow<2>(parameters.y));
  elements.period = parameters.period;
  elements.argument_of_periapsis = ArcTan(parameters.y, parameters.x);
  elements.mean_anomaly =
      π / 2 * Radian - *elements.argument_of_periapsis -
      (2 * π * Radian) * parameters.time_to_first_transit / *elements.period;
  return elements;
}

PlanetParameters MakePlanetParameters(
    KeplerianElements<Sky> const& elements) {
  PlanetParameters result;
  result.period = *elements.period;
  result.x = *elements.eccentricity * Cos(*elements.argument_of_periapsis);
  result.y = *elements.eccentricity * Sin(*elements.argument_of_periapsis);
  result.time_to_first_transit =
      *elements.period / (2 * π * Radian) *
      (π / 2 * Radian - *elements.argument_of_periapsis -
       *elements.mean_anomaly);
  return result;
}

std::vector<double> EvaluatePopulation(
    Population const& population,
    ComputeLogPdf const& compute_log_pdf,
    std::vector<std::string>& info) {
  std::vector<double> log_pdf(population.size());
  info.resize(population.size());
  Bundle bundle;
  for (int i = 0; i < population.size(); ++i) {
    auto const& parameters = population[i];
    bundle.Add([&compute_log_pdf, i, &log_pdf, &parameters, &info]() {
      log_pdf[i] = compute_log_pdf(parameters, info[i]);
      return Status::OK;
    });
  }
  bundle.Join();
  return log_pdf;
}

Population GenerateTrialStates(Population const& population,
                               double const γ,
                               double const ε,
                               std::mt19937_64& engine) {
  Population trial(population.size());
  std::uniform_int_distribution<> j_distribution(0, population.size() - 2);
  std::uniform_int_distribution<> k_distribution(0, population.size() - 3);
  std::normal_distribution<> perturbation_distribution;
  for (int i = 0; i < population.size(); ++i) {
    // Choose head (k) and tail (j) for perturbation vector.
    int j = j_distribution(engine);
    if (j >= i) {
      ++j;
    }
    int k = k_distribution(engine);
    if (k >= i) {
      ++k;
    }
    if (k >= j) {
      ++k;
    }

    // Choose scale factor.
    double const scale = (1.0 + ε * perturbation_distribution(engine)) * γ;

    trial[i] = population[i] + scale * (population[k] - population[j]);
  }
  return trial;
}

SystemParameters Run(Population& population,
                     int const number_of_generations,
                     int const number_of_generations_between_kicks,
                     int const number_of_burn_in_generations,
                     double const ε,
                     ComputeLogPdf const& compute_log_pdf) {
  CHECK_LE(1, number_of_generations);
  CHECK_LT(std::tuple_size<SystemParameters>::value * PlanetParameters::count,
           population.size());
  std::mt19937_64 engine;
  std::uniform_real_distribution<> distribution(0.0, 1.0);

  static double best_log_pdf = -std::numeric_limits<double>::infinity();
  std::string best_info;
  SystemParameters best_system_parameters;

  std::vector<std::string> infos;
  auto log_pdf = EvaluatePopulation(population, compute_log_pdf, infos);

  // Loop over generations.
  for (int generation = 0; generation < number_of_generations; ++generation) {
    int accepted = 0;

    // Every 10th generation try full-size steps.
    double const
        γ = generation < number_of_burn_in_generations
                ? 0.01
                : generation % number_of_generations_between_kicks == 0
                      ? 1.0
                      : 2.38 /
                            Sqrt(2 * std::tuple_size<SystemParameters>::value *
                                 PlanetParameters::count);

    // Evaluate model for each set of trial parameters.
    auto const trial = GenerateTrialStates(population, γ, ε, engine);
    std::vector<std::string> trial_infos;
    auto const log_pdf_trial =
        EvaluatePopulation(trial, compute_log_pdf, trial_infos);

    // For each member of population.
    for (int i = 0; i < population.size(); ++i) {
      double const log_pdf_ratio = log_pdf_trial[i] - log_pdf[i];
      if (log_pdf_ratio > 0.0 ||
          log_pdf_ratio > std::log(distribution(engine))) {
        population[i] = trial[i];
        log_pdf[i] = log_pdf_trial[i];
        infos[i] = trial_infos[i];
        ++accepted;
      }
    }

    // Traces.
    int const max_index =
        std::max_element(log_pdf.begin(), log_pdf.end()) - log_pdf.begin();
    if (best_log_pdf < log_pdf[max_index]) {
      best_system_parameters = population[max_index];
      best_log_pdf = log_pdf[max_index];
      best_info = infos[max_index];
    }
    LOG(ERROR) << "Generation " << generation << "; Acceptance: " << accepted
               << " / " << population.size();
    LOG(ERROR) << "Max  : " << infos[max_index];
    if (best_info != infos[max_index]) {
      LOG(ERROR) << "Best : " << best_info;
    }
  }
  return best_system_parameters;
}

}  // namespace deмcmc

double ShortDays(Instant const& time) {
    return (time - "JD2450000.0"_TT) / Day;
}

using Transits = std::vector<Instant>;
using TransitsByPlanet = std::map<std::string, Transits>;

template<typename Value>
struct Measured {
  Value estimated_value;
  Difference<Value> standard_uncertainty;
};

using MeasuredTransits = std::vector<Measured<Instant>>;
using MeasuredTransitsByPlanet = std::map<std::string, MeasuredTransits>;

std::map<std::string, Time> nominal_periods = {
    {"Trappist-1b", 1.51087637 * Day},
    {"Trappist-1c", 2.42180746 * Day},
    {"Trappist-1d", 4.049959 * Day},
    {"Trappist-1e", 6.099043 * Day},
    {"Trappist-1f", 9.205585 * Day},
    {"Trappist-1g", 12.354473 * Day},
    {"Trappist-1h", 18.767953 * Day}};

MeasuredTransitsByPlanet const observations = {
    {"Trappist-1b", {{"JD2457322.51531"_TT, 0.00071 * Day},
                     {"JD2457325.53910"_TT, 0.00100 * Day},
                     {"JD2457328.55860"_TT, 0.00130 * Day},
                     {"JD2457331.58160"_TT, 0.00100 * Day},
                     {"JD2457334.60480"_TT, 0.00017 * Day},
                     {"JD2457337.62644"_TT, 0.00092 * Day},
                     {"JD2457340.64820"_TT, 0.00140 * Day},
                     {"JD2457345.18028"_TT, 0.00080 * Day},
                     {"JD2457361.79945"_TT, 0.00028 * Day},
                     {"JD2457364.82173"_TT, 0.00077 * Day},
                     {"JD2457440.36492"_TT, 0.00020 * Day},
                     {"JD2457452.45228"_TT, 0.00014 * Day},
                     {"JD2457463.02847"_TT, 0.00019 * Day},
                     {"JD2457509.86460"_TT, 0.00210 * Day},
                     {"JD2457512.88731"_TT, 0.00029 * Day},
                     {"JD2457568.78880"_TT, 0.00100 * Day},
                     {"JD2457586.91824"_TT, 0.00064 * Day},
                     {"JD2457589.93922"_TT, 0.00092 * Day},
                     {"JD2457599.00640"_TT, 0.00021 * Day},
                     {"JD2457602.02805"_TT, 0.00071 * Day},
                     {"JD2457612.60595"_TT, 0.00085 * Day},
                     {"JD2457615.62710"_TT, 0.00160 * Day},
                     {"JD2457624.69094"_TT, 0.00066 * Day},
                     {"JD2457645.84400"_TT, 0.00110 * Day},
                     {"JD2457651.88743"_TT, 0.00022 * Day},
                     {"JD2457653.39809"_TT, 0.00026 * Day},
                     {"JD2457654.90908"_TT, 0.00084 * Day},
                     {"JD2457656.41900"_TT, 0.00029 * Day},
                     {"JD2457657.93129"_TT, 0.00020 * Day},
                     {"JD2457659.44144"_TT, 0.00017 * Day},
                     {"JD2457660.95205"_TT, 0.00035 * Day},
                     {"JD2457662.46358"_TT, 0.00020 * Day},
                     {"JD2457663.97492"_TT, 0.00070 * Day},
                     {"JD2457665.48509"_TT, 0.00017 * Day},
                     {"JD2457666.99567"_TT, 0.00025 * Day},
                     {"JD2457668.50668"_TT, 0.00030 * Day},
                     {"JD2457670.01766"_TT, 0.00034 * Day},
                     {"JD2457671.52876"_TT, 0.00033 * Day},
                     {"JD2457721.38747"_TT, 0.00035 * Day},
                     {"JD2457739.51770"_TT, 0.00059 * Day},
                     {"JD2457741.02787"_TT, 0.00055 * Day},
                     {"JD2457742.53918"_TT, 0.00058 * Day},
                     {"JD2457744.05089"_TT, 0.00061 * Day},
                     {"JD2457745.56164"_TT, 0.00072 * Day},
                     {"JD2457747.07208"_TT, 0.00085 * Day},
                     {"JD2457748.58446"_TT, 0.00087 * Day},
                     {"JD2457750.09387"_TT, 0.00089 * Day},
                     {"JD2457751.60535"_TT, 0.00082 * Day},
                     {"JD2457753.11623"_TT, 0.00075 * Day},
                     {"JD2457754.62804"_TT, 0.00077 * Day},
                     {"JD2457756.13856"_TT, 0.00060 * Day},
                     {"JD2457757.64840"_TT, 0.00089 * Day},
                     {"JD2457759.15953"_TT, 0.00073 * Day},
                     {"JD2457760.67112"_TT, 0.00082 * Day},
                     {"JD2457762.18120"_TT, 0.00073 * Day},
                     {"JD2457763.69221"_TT, 0.00071 * Day},
                     {"JD2457765.20298"_TT, 0.00077 * Day},
                     {"JD2457766.71479"_TT, 0.00055 * Day},
                     {"JD2457768.22514"_TT, 0.00103 * Day},
                     {"JD2457769.73704"_TT, 0.00064 * Day},
                     {"JD2457771.24778"_TT, 0.00091 * Day},
                     {"JD2457772.75738"_TT, 0.00075 * Day},
                     {"JD2457774.26841"_TT, 0.00080 * Day},
                     {"JD2457775.77995"_TT, 0.00058 * Day},
                     {"JD2457777.28899"_TT, 0.00099 * Day},
                     {"JD2457778.80118"_TT, 0.00062 * Day},
                     {"JD2457780.31297"_TT, 0.00068 * Day},
                     {"JD2457781.82231"_TT, 0.00145 * Day},
                     {"JD2457783.33410"_TT, 0.00071 * Day},
                     {"JD2457784.84372"_TT, 0.00068 * Day},
                     {"JD2457792.39979"_TT, 0.00110 * Day},
                     {"JD2457793.90955"_TT, 0.00064 * Day},
                     {"JD2457795.41987"_TT, 0.00058 * Day},
                     {"JD2457796.93134"_TT, 0.00065 * Day},
                     {"JD2457798.44211"_TT, 0.00061 * Day},
                     {"JD2457799.95320"_TT, 0.00083 * Day},
                     {"JD2457801.46314"_TT, 0.00127 * Day},
                     {"JD2457802.97557"_TT, 0.00016 * Day},
                     {"JD2457804.48638"_TT, 0.00053 * Day},
                     {"JD2457805.99697"_TT, 0.00016 * Day},
                     {"JD2457807.50731"_TT, 0.00017 * Day},
                     {"JD2457809.01822"_TT, 0.00017 * Day},
                     {"JD2457810.52781"_TT, 0.00110 * Day},
                     {"JD2457812.04038"_TT, 0.00020 * Day},
                     {"JD2457813.55121"_TT, 0.00014 * Day},
                     {"JD2457815.06275"_TT, 0.00017 * Day},
                     {"JD2457816.57335"_TT, 0.00011 * Day},
                     {"JD2457818.08382"_TT, 0.00015 * Day},
                     {"JD2457819.59478"_TT, 0.00017 * Day},
                     {"JD2457821.10550"_TT, 0.00020 * Day},
                     {"JD2457824.12730"_TT, 0.00018 * Day},
                     {"JD2457825.63813"_TT, 0.00018 * Day},
                     {"JD2457827.14995"_TT, 0.00012 * Day},
                     {"JD2457828.66042"_TT, 0.00024 * Day},
                     {"JD2457830.17087"_TT, 0.00021 * Day},
                     {"JD2457833.19257"_TT, 0.00018 * Day},
                     {"JD2457834.70398"_TT, 0.00016 * Day},
                     {"JD2457836.21440"_TT, 0.00017 * Day},
                     {"JD2457837.72526"_TT, 0.00014 * Day},
                     {"JD2457839.23669"_TT, 0.00017 * Day},
                     {"JD2457917.80060"_TT, 0.00110 * Day},
                     {"JD2457923.84629"_TT, 0.00045 * Day},
                     {"JD2457935.93288"_TT, 0.00023 * Day},
                     {"JD2457952.55450"_TT, 0.00110 * Day},
                     {"JD2457955.57554"_TT, 0.00069 * Day},
                     {"JD2457967.66254"_TT, 0.00050 * Day},
                     {"JD2457973.70596"_TT, 0.00040 * Day}}},
    {"Trappist-1c", {{"JD2457282.80570"_TT, 0.00140 * Day},
                     {"JD2457333.66400"_TT, 0.00090 * Day},
                     {"JD2457362.72605"_TT, 0.00038 * Day},
                     {"JD2457367.57051"_TT, 0.00033 * Day},
                     {"JD2457384.52320"_TT, 0.00130 * Day},
                     {"JD2457452.33470"_TT, 0.00015 * Day},
                     {"JD2457454.75672"_TT, 0.00066 * Day},
                     {"JD2457512.88094"_TT, 0.00009 * Day},
                     {"JD2457546.78587"_TT, 0.00075 * Day},
                     {"JD2457551.62888"_TT, 0.00066 * Day},
                     {"JD2457580.69137"_TT, 0.00031 * Day},
                     {"JD2457585.53577"_TT, 0.00250 * Day},
                     {"JD2457587.95622"_TT, 0.00054 * Day},
                     {"JD2457600.06684"_TT, 0.00036 * Day},
                     {"JD2457604.90975"_TT, 0.00063 * Day},
                     {"JD2457609.75461"_TT, 0.00072 * Day},
                     {"JD2457614.59710"_TT, 0.00130 * Day},
                     {"JD2457626.70610"_TT, 0.00110 * Day},
                     {"JD2457631.55024"_TT, 0.00056 * Day},
                     {"JD2457638.81518"_TT, 0.00048 * Day},
                     {"JD2457650.92395"_TT, 0.00023 * Day},
                     {"JD2457653.34553"_TT, 0.00024 * Day},
                     {"JD2457655.76785"_TT, 0.00043 * Day},
                     {"JD2457658.18963"_TT, 0.00024 * Day},
                     {"JD2457660.61168"_TT, 0.00051 * Day},
                     {"JD2457663.03292"_TT, 0.00028 * Day},
                     {"JD2457665.45519"_TT, 0.00025 * Day},
                     {"JD2457667.87729"_TT, 0.00031 * Day},
                     {"JD2457670.29869"_TT, 0.00035 * Day},
                     {"JD2457672.71944"_TT, 0.00081 * Day},
                     {"JD2457711.46778"_TT, 0.00064 * Day},
                     {"JD2457723.57663"_TT, 0.00050 * Day},
                     {"JD2457740.53361"_TT, 0.00088 * Day},
                     {"JD2457742.95276"_TT, 0.00115 * Day},
                     {"JD2457745.37429"_TT, 0.00063 * Day},
                     {"JD2457747.79699"_TT, 0.00056 * Day},
                     {"JD2457750.21773"_TT, 0.00096 * Day},
                     {"JD2457752.64166"_TT, 0.00093 * Day},
                     {"JD2457755.05877"_TT, 0.00165 * Day},
                     {"JD2457757.48313"_TT, 0.00066 * Day},
                     {"JD2457759.90281"_TT, 0.00058 * Day},
                     {"JD2457762.32806"_TT, 0.00081 * Day},
                     {"JD2457764.74831"_TT, 0.00072 * Day},
                     {"JD2457767.16994"_TT, 0.00125 * Day},
                     {"JD2457769.59209"_TT, 0.00081 * Day},
                     {"JD2457772.01483"_TT, 0.00100 * Day},
                     {"JD2457774.43458"_TT, 0.00081 * Day},
                     {"JD2457776.85815"_TT, 0.00102 * Day},
                     {"JD2457779.27911"_TT, 0.00089 * Day},
                     {"JD2457781.70095"_TT, 0.00072 * Day},
                     {"JD2457784.12338"_TT, 0.00054 * Day},
                     {"JD2457791.38801"_TT, 0.00064 * Day},
                     {"JD2457793.81141"_TT, 0.00079 * Day},
                     {"JD2457796.23153"_TT, 0.00052 * Day},
                     {"JD2457798.65366"_TT, 0.00082 * Day},
                     {"JD2457801.07631"_TT, 0.00084 * Day},
                     {"JD2457803.49747"_TT, 0.00020 * Day},
                     {"JD2457805.91882"_TT, 0.00017 * Day},
                     {"JD2457808.34123"_TT, 0.00023 * Day},
                     {"JD2457810.76273"_TT, 0.00019 * Day},
                     {"JD2457813.18456"_TT, 0.00024 * Day},
                     {"JD2457815.60583"_TT, 0.00017 * Day},
                     {"JD2457818.02821"_TT, 0.00020 * Day},
                     {"JD2457820.45019"_TT, 0.00022 * Day},
                     {"JD2457822.87188"_TT, 0.00021 * Day},
                     {"JD2457825.29388"_TT, 0.00022 * Day},
                     {"JD2457827.71513"_TT, 0.00022 * Day},
                     {"JD2457830.13713"_TT, 0.00026 * Day},
                     {"JD2457832.55888"_TT, 0.00015 * Day},
                     {"JD2457834.98120"_TT, 0.00025 * Day},
                     {"JD2457837.40280"_TT, 0.00017 * Day},
                     {"JD2457839.82415"_TT, 0.00031 * Day}}},
    {"Trappist-1d", {{"JD2457560.79730"_TT, 0.00230 * Day},
                     {"JD2457625.59779"_TT, 0.00078 * Day},
                     {"JD2457641.79360"_TT, 0.00290 * Day},
                     {"JD2457645.84360"_TT, 0.00210 * Day},
                     {"JD2457653.94261"_TT, 0.00051 * Day},
                     {"JD2457657.99220"_TT, 0.00063 * Day},
                     {"JD2457662.04284"_TT, 0.00051 * Day},
                     {"JD2457666.09140"_TT, 0.00160 * Day},
                     {"JD2457670.14198"_TT, 0.00066 * Day},
                     {"JD2457726.83975"_TT, 0.00029 * Day},
                     {"JD2457738.99169"_TT, 0.00160 * Day},
                     {"JD2457743.03953"_TT, 0.00180 * Day},
                     {"JD2457747.08985"_TT, 0.00145 * Day},
                     {"JD2457751.14022"_TT, 0.00195 * Day},
                     {"JD2457755.18894"_TT, 0.00155 * Day},
                     {"JD2457759.24638"_TT, 0.00225 * Day},
                     {"JD2457763.28895"_TT, 0.00150 * Day},
                     {"JD2457767.33866"_TT, 0.00190 * Day},
                     {"JD2457771.39077"_TT, 0.00260 * Day},
                     {"JD2457775.44026"_TT, 0.00125 * Day},
                     {"JD2457779.48843"_TT, 0.00190 * Day},
                     {"JD2457783.54023"_TT, 0.00240 * Day},
                     {"JD2457791.64083"_TT, 0.00135 * Day},
                     {"JD2457803.79083"_TT, 0.00049 * Day},
                     {"JD2457807.84032"_TT, 0.00030 * Day},
                     {"JD2457811.89116"_TT, 0.00050 * Day},
                     {"JD2457815.94064"_TT, 0.00030 * Day},
                     {"JD2457819.99050"_TT, 0.00050 * Day},
                     {"JD2457824.04185"_TT, 0.00067 * Day},
                     {"JD2457828.09082"_TT, 0.00043 * Day},
                     {"JD2457832.14036"_TT, 0.00037 * Day},
                     {"JD2457836.19171"_TT, 0.00042 * Day},
                     {"JD2457961.73760"_TT, 0.00130 * Day},
                     {"JD2457969.83708"_TT, 0.00068 * Day},
                     {"JD2457973.88590"_TT, 0.00066 * Day}}},
    {"Trappist-1e", {{"JD2457312.71300"_TT, 0.00270 * Day},
                     {"JD2457367.59683"_TT, 0.00037 * Day},
                     {"JD2457611.57620"_TT, 0.00310 * Day},
                     {"JD2457623.77950"_TT, 0.00100 * Day},
                     {"JD2457654.27862"_TT, 0.00049 * Day},
                     {"JD2457660.38016"_TT, 0.00078 * Day},
                     {"JD2457666.48030"_TT, 0.00180 * Day},
                     {"JD2457672.57930"_TT, 0.00260 * Day},
                     {"JD2457721.37514"_TT, 0.00099 * Day},
                     {"JD2457733.57300"_TT, 0.00140 * Day},
                     {"JD2457739.67085"_TT, 0.00135 * Day},
                     {"JD2457745.77160"_TT, 0.00120 * Day},
                     {"JD2457751.87007"_TT, 0.00034 * Day},
                     {"JD2457757.96712"_TT, 0.00160 * Day},
                     {"JD2457764.06700"_TT, 0.00240 * Day},
                     {"JD2457770.17109"_TT, 0.00215 * Day},
                     {"JD2457776.26378"_TT, 0.00160 * Day},
                     {"JD2457782.36226"_TT, 0.00175 * Day},
                     {"JD2457794.56159"_TT, 0.00160 * Day},
                     {"JD2457800.66354"_TT, 0.00170 * Day},
                     {"JD2457806.75758"_TT, 0.00041 * Day},
                     {"JD2457812.85701"_TT, 0.00034 * Day},
                     {"JD2457818.95510"_TT, 0.00030 * Day},
                     {"JD2457825.05308"_TT, 0.00035 * Day},
                     {"JD2457831.15206"_TT, 0.00027 * Day},
                     {"JD2457837.24980"_TT, 0.00025 * Day},
                     {"JD2457934.83095"_TT, 0.00050 * Day},
                     {"JD2457940.92995"_TT, 0.00086 * Day}}},
    {"Trappist-1f", {{"JD2457321.52520"_TT, 0.00200 * Day},
                     {"JD2457367.57629"_TT, 0.00044 * Day},
                     {"JD2457634.57809"_TT, 0.00061 * Day},
                     {"JD2457652.98579"_TT, 0.00032 * Day},
                     {"JD2457662.18747"_TT, 0.00040 * Day},
                     {"JD2457671.39279"_TT, 0.00072 * Day},
                     {"JD2457717.41541"_TT, 0.00091 * Day},
                     {"JD2457726.61960"_TT, 0.00026 * Day},
                     {"JD2457745.03116"_TT, 0.00135 * Day},
                     {"JD2457754.23380"_TT, 0.00155 * Day},
                     {"JD2457763.44338"_TT, 0.00024 * Day},
                     {"JD2457772.64752"_TT, 0.00160 * Day},
                     {"JD2457781.85142"_TT, 0.00180 * Day},
                     {"JD2457800.27307"_TT, 0.00140 * Day},
                     {"JD2457809.47554"_TT, 0.00027 * Day},
                     {"JD2457818.68271"_TT, 0.00032 * Day},
                     {"JD2457827.88669"_TT, 0.00030 * Day},
                     {"JD2457837.10322"_TT, 0.00032 * Day},
                     {"JD2457956.80549"_TT, 0.00054 * Day}}},
    {"Trappist-1g", {{"JD2457294.78600"_TT, 0.00390 * Day},
                     {"JD2457356.53410"_TT, 0.00200 * Day},
                     {"JD2457615.92400"_TT, 0.00170 * Day},
                     {"JD2457640.63730"_TT, 0.00100 * Day},
                     {"JD2457652.99481"_TT, 0.00030 * Day},
                     {"JD2457665.35151"_TT, 0.00028 * Day},
                     {"JD2457739.48441"_TT, 0.00115 * Day},
                     {"JD2457751.83993"_TT, 0.00017 * Day},
                     {"JD2457764.19098"_TT, 0.00155 * Day},
                     {"JD2457776.54900"_TT, 0.00110 * Day},
                     {"JD2457801.25000"_TT, 0.00093 * Day},
                     {"JD2457813.60684"_TT, 0.00023 * Day},
                     {"JD2457825.96112"_TT, 0.00020 * Day},
                     {"JD2457838.30655"_TT, 0.00028 * Day},
                     {"JD2457924.77090"_TT, 0.00140 * Day},
                     {"JD2457961.82621"_TT, 0.00068 * Day}}},
    {"Trappist-1h", {{"JD2457662.55467"_TT, 0.00054 * Day},
                     {"JD2457756.38740"_TT, 0.00130 * Day},
                     {"JD2457775.15390"_TT, 0.00160 * Day},
                     {"JD2457793.92300"_TT, 0.00250 * Day},
                     {"JD2457812.69870"_TT, 0.00450 * Day},
                     {"JD2457831.46625"_TT, 0.00047 * Day},
                     {"JD2457962.86271"_TT, 0.00083 * Day}}}};

class TrappistDynamicsTest : public ::testing::Test {
 protected:
  TrappistDynamicsTest()
      : system_(SOLUTION_DIR / "astronomy" / "trappist_gravity_model.proto.txt",
                SOLUTION_DIR / "astronomy" /
                    "trappist_initial_state_jd_2457000_000000000.proto.txt"),
        ephemeris_(system_.MakeEphemeris(
            /*fitting_tolerance=*/5 * Milli(Metre),
            Ephemeris<Sky>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                                   Position<Sky>>(),
                /*step=*/0.07 * Day))) {}

  static Transits ComputeTransits(Ephemeris<Sky> const& ephemeris,
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
      RelativeDegreesOfFreedom<Sky> const relative_dof =
          planet_trajectory->EvaluateDegreesOfFreedom(t) -
          star_trajectory->EvaluateDegreesOfFreedom(t);

      auto const xy_displacement_derivative =
          [&planet_trajectory, &star_trajectory](Instant const& t) {
            RelativeDegreesOfFreedom<Sky> const relative_dof =
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

  static double Transitsχ²(MeasuredTransitsByPlanet const& observations,
                           TransitsByPlanet const& computations,
                           std::string& info) {
    double sum_of_squared_errors = 0.0;
    Time max_Δt;
    Time total_Δt;
    std::string transit_with_max_Δt;
    int number_of_observations = 0;
    for (auto const& pair : observations) {
      auto const& name = pair.first;
      auto const& observed_transits = pair.second;
      auto const& computed_transits = computations.at(name);
      if (computed_transits.empty()) {
        return std::numeric_limits<double>::infinity();
      }
      Instant const& initial_observed_transit =
          observed_transits.front().estimated_value;
      auto initial_computed_transit =
          std::lower_bound(computed_transits.begin(),
                           computed_transits.end(),
                           initial_observed_transit);
      if (initial_computed_transit == computed_transits.end()) {
        --initial_computed_transit;
      } else if (initial_computed_transit != computed_transits.begin() &&
                 *initial_computed_transit - initial_observed_transit >
                     initial_observed_transit - initial_computed_transit[-1]) {
        --initial_computed_transit;
      }
      int const relevant_computed_transits_size =
          computed_transits.end() - initial_computed_transit;
      for (auto const& observed_transit : observed_transits) {
        int const transit_epoch = std::round(
            (observed_transit.estimated_value - initial_observed_transit) /
            nominal_periods.at(name));
        if (transit_epoch >= relevant_computed_transits_size) {
          // No computed transit corresponds to the observed transit.  Either
          // the planet has escaped, or its period is so low that it does not
          // transit enough over the simulation interval.  In any case,
          // something is very wrong.
          return std::numeric_limits<double>::infinity();
        }
        auto const computed_transit = initial_computed_transit[transit_epoch];
        Time const error =
            Abs(computed_transit - observed_transit.estimated_value);
        if (error > max_Δt) {
          max_Δt = error;
          transit_with_max_Δt =
              name + "[" +
              std::to_string(ShortDays(observed_transit.estimated_value)) + "]";
        }
        total_Δt += error;
        ++number_of_observations;
        sum_of_squared_errors +=
            Pow<2>(error / observed_transit.standard_uncertainty);
      }
    }
    info = u8"max Δt = " + std::to_string(max_Δt / Second) + " s (" +
           transit_with_max_Δt + u8") avg Δt = " +
           std::to_string(total_Δt / number_of_observations / Second) + " s";
    return sum_of_squared_errors;
  }

  static double ProlongAndComputeTransitsχ²(SolarSystem<Sky>& system,
                                            std::string& info) {
    auto const ephemeris = system.MakeEphemeris(
        /*fitting_tolerance=*/5 * Milli(Metre),
        Ephemeris<Sky>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                               Position<Sky>>(),
            /*step=*/0.07 * Day));
    ephemeris->Prolong(system.epoch() + 1000 * Day);

    // For some combinations we get an apocalyse.  In this case the dispersion
    // is infinite.
    if (!ephemeris->last_severe_integration_status().ok()) {
      info = u8"ἀποκάλυψις";
      return std::numeric_limits<double>::infinity();
    }

    TransitsByPlanet computations;
    auto const& star = system.massive_body(*ephemeris, star_name);
    auto const bodies = ephemeris->bodies();
    for (auto const& planet : bodies) {
      if (planet != star) {
        computations[planet->name()] =
            ComputeTransits(*ephemeris, star, planet);
      }
    }
    std::string χ²_info;
    double const χ² = Transitsχ²(observations, computations, χ²_info);
    info = u8"χ² = " + std::to_string(χ²) + " " + χ²_info;
    return χ²;
  }

  static std::string SanitizedName(MassiveBody const& body) {
    auto sanitized_name = body.name();
    return sanitized_name.erase(sanitized_name.find_first_of("-"), 1);
  }

  constexpr static char home_name[] = "Trappist-1e";
  constexpr static char star_name[] = "Trappist-1";
  SolarSystem<Sky> const system_;
  not_null<std::unique_ptr<Ephemeris<Sky>>> ephemeris_;
};

constexpr char TrappistDynamicsTest::home_name[];
constexpr char TrappistDynamicsTest::star_name[];

#if !defined(_DEBUG)
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
        KeplerOrbit<Sky> const planet_orbit(
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

  std::string info;
  double const χ² = Transitsχ²(observations, computations, info);
  CHECK_LT(χ², 482.0);
  CHECK_GT(χ², 470.0);
  LOG(ERROR) << u8"χ²: " << χ² << " " << info;
}

TEST_F(TrappistDynamicsTest, MathematicaAlignments) {
  Instant const a_century_later = system_.epoch() + 100 * JulianYear;
  ephemeris_->Prolong(a_century_later);

  OFStream file(TEMP_DIR / "trappist_alignments.generated.wl");

  auto const& star = system_.massive_body(*ephemeris_, star_name);
  auto const& star_trajectory = ephemeris_->trajectory(star);
  auto const& home = system_.rotating_body(*ephemeris_, home_name);
  auto const& home_trajectory = ephemeris_->trajectory(home);

  auto const bodies = ephemeris_->bodies();
  for (auto const& planet : bodies) {
    if (planet != star && planet != &*home) {
      auto const& planet_trajectory = ephemeris_->trajectory(planet);
      std::vector<double> views;
      for (Instant t = ephemeris_->t_min();
           t < ephemeris_->t_min() + 1000 * Hour;
           t += 10 * Minute) {
        Angle const view =
            OrientedAngleBetween(star_trajectory->EvaluatePosition(t) -
                                     home_trajectory->EvaluatePosition(t),
                                 planet_trajectory->EvaluatePosition(t) -
                                     home_trajectory->EvaluatePosition(t),
                                 home->angular_velocity());
        views.push_back(view / Degree);
      }
      file << mathematica::Assign("alignment" + SanitizedName(*planet), views);
    }
  }
}

TEST_F(TrappistDynamicsTest, PlanetBPlanetDAlignment) {
  Instant const a_century_later = system_.epoch() + 100 * JulianYear;
  ephemeris_->Prolong(a_century_later);

  auto const& star = system_.massive_body(*ephemeris_, star_name);
  auto const& star_trajectory = ephemeris_->trajectory(star);
  auto const& home = system_.rotating_body(*ephemeris_, home_name);
  auto const& home_trajectory = ephemeris_->trajectory(home);
  auto const& planet_b = system_.massive_body(*ephemeris_, "Trappist-1b");
  auto const& planet_b_trajectory = ephemeris_->trajectory(planet_b);
  auto const& planet_d = system_.massive_body(*ephemeris_, "Trappist-1d");
  auto const& planet_d_trajectory = ephemeris_->trajectory(planet_d);

  Angle min_angle = 1 * Radian;

  for (Instant t = ephemeris_->t_min();
        t < ephemeris_->t_min() + 200000 * Hour;
        t += 10 * Minute) {
    Displacement<Sky> const home_star =
        star_trajectory->EvaluatePosition(t) -
        home_trajectory->EvaluatePosition(t);
    Displacement<Sky> const home_planet_b =
        planet_b_trajectory->EvaluatePosition(t) -
        home_trajectory->EvaluatePosition(t);
    Displacement<Sky> const planet_b_star =
        star_trajectory->EvaluatePosition(t) -
        planet_b_trajectory->EvaluatePosition(t);
    Displacement<Sky> const home_planet_d =
        planet_d_trajectory->EvaluatePosition(t) -
        home_trajectory->EvaluatePosition(t);
    Displacement<Sky> const planet_d_star =
        star_trajectory->EvaluatePosition(t) -
        planet_d_trajectory->EvaluatePosition(t);
    Angle const star_b = AngleBetween(home_star, home_planet_b);
    Angle const star_d = AngleBetween(home_star, home_planet_d);
    Angle const angle = Sqrt(star_b * star_b + star_d * star_d);
    bool in_front =
        InnerProduct(home_star, planet_b_star) > 0.0 * Metre * Metre &&
        InnerProduct(home_star, planet_d_star) > 0.0 * Metre * Metre;
    if (in_front && angle < min_angle) {
      LOG(ERROR) << "Time: " << t << ", Angle: " << angle / Degree << u8"°";
      min_angle = angle;
    }
  }
}
#endif

TEST_F(TrappistDynamicsTest, DISABLED_Optimization) {
  SolarSystem<Sky> const system(
      SOLUTION_DIR / "astronomy" / "trappist_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "trappist_preoptimization_initial_state_jd_2457000_000000000"
          ".proto.txt");

  auto planet_names = system.names();
  planet_names.erase(
      std::find(planet_names.begin(), planet_names.end(), star_name));
  std::vector<KeplerianElements<Sky>> elements;
  for (auto const& planet_name : planet_names) {
    elements.push_back(SolarSystem<Sky>::MakeKeplerianElements(
        system.keplerian_initial_state_message(planet_name).elements()));
  }

  auto compute_fitness =
      [&planet_names, &system](genetics::Genome const& genome,
                               std::string& info) {
        auto modified_system = system;
        auto const& elements = genome.elements();
        for (int i = 0; i < planet_names.size(); ++i) {
          modified_system.ReplaceElements(planet_names[i], elements[i]);
        }
        double const χ² = ProlongAndComputeTransitsχ²(modified_system, info);

        // This is the place where we cook the sausage.  This function must be
        // steep enough to efficiently separate the wheat from the chaff without
        // leading to monoculture.
        return 1 / χ²;
      };

  auto compute_log_pdf =
      [&elements, &planet_names, &system](
          deмcmc::SystemParameters const& system_parameters,
          std::string& info) {
        auto modified_system = system;
        for (int i = 0; i < planet_names.size(); ++i) {
          modified_system.ReplaceElements(
              planet_names[i],
              MakeKeplerianElements(elements[i],
                                    system_parameters[i]));
        }
        double const χ² = ProlongAndComputeTransitsχ²(modified_system, info);
        return -χ² / 2.0;
      };

  std::optional<genetics::Genome> great_old_one;
  double great_old_one_fitness = 0.0;
  {
    // First, let's do 5 rounds of evolution with a population of 9 individuals
    // based on |luca|.  The best of all of them is the Great Old One.
    std::mt19937_64 engine;
    genetics::Genome luca(elements);
    for (int i = 0; i < 5; ++i) {
      genetics::Population population(luca,
                                      9,
                                      /*elitism=*/true,
                                      compute_fitness,
                                      engine);
      for (int i = 0; i < 20'000; ++i) {
        population.ComputeAllFitnesses();
        population.BegetChildren();
      }
      LOG(ERROR) << "Great Old One #" << i;
      LOG(ERROR) << population.best_genome_trace();
      for (int i = 0; i < planet_names.size(); ++i) {
        LOG(ERROR) << planet_names[i] << ": "
                   << population.best_genome().elements()[i];
      }
      if (population.best_genome_fitness() > great_old_one_fitness) {
        great_old_one = population.best_genome();
        great_old_one_fitness = population.best_genome_fitness();
      }
    }
  }
  {
    // Next, let's build a population of 50 minor variants of the Great Old One,
    // the Outer Gods.  Use DEMCMC to improve them.  The best of them is the
    // Blind Idiot God.
    std::mt19937_64 engine;
    deмcmc::Population outer_gods;
    for (int i = 0; i < 50; ++i) {
      outer_gods.emplace_back();
      deмcmc::SystemParameters& outer_god = outer_gods.back();
      std::normal_distribution<> angle_distribution(0.0, 0.1);
      std::normal_distribution<> period_distribution(0.0, 1.0);
      std::normal_distribution<> eccentricity_distribution(0.0, 1e-4);
      for (int j = 0; j < outer_god.size(); ++j) {
        auto perturbed_elements = great_old_one->elements()[j];
        *perturbed_elements.period += period_distribution(engine) * Second;
        *perturbed_elements.argument_of_periapsis +=
            angle_distribution(engine) * Degree;
        *perturbed_elements.mean_anomaly += angle_distribution(engine) * Degree;
        *perturbed_elements.eccentricity += eccentricity_distribution(engine);
        outer_god[j] = deмcmc::MakePlanetParameters(perturbed_elements);
      }
    }

    auto const the_blind_idiot_god =
        deмcmc::Run(outer_gods,
                    /*number_of_generations=*/10'000,
                    /*number_of_generations_between_kicks=*/30,
                    /*number_of_burn_in_generations=*/10,
                    /*ε=*/0.05,
                    compute_log_pdf);
    LOG(ERROR) << "The Blind Idiot God";
    for (int i = 0; i < the_blind_idiot_god.size(); ++i) {
      LOG(ERROR) << planet_names[i];
      auto const elements = MakeKeplerianElements(great_old_one->elements()[i],
                                                  the_blind_idiot_god[i]);
      LOG(ERROR) << std::setprecision(
                        std::numeric_limits<double>::max_digits10)
                 << "        eccentricity                : "
                 << *elements.eccentricity;
      LOG(ERROR) << std::setprecision(
                        std::numeric_limits<double>::max_digits10)
                 << "        period                      : \""
                 << *elements.period / Day << " d\"";
      LOG(ERROR) << std::setprecision(
                        std::numeric_limits<double>::max_digits10)
                 << "        argument_of_periapsis       : \""
                 << Mod(*elements.argument_of_periapsis, 2 * π * Radian) /
                        Degree
                 << " deg\"";
      LOG(ERROR) << std::setprecision(
                        std::numeric_limits<double>::max_digits10)
                 << "        mean_anomaly                : \""
                 << Mod(*elements.mean_anomaly, 2 * π * Radian) / Degree
                 << " deg\"";
    }
  }
}

}  // namespace astronomy
}  // namespace principia
