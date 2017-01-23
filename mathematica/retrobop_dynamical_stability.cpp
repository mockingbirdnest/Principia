
#include "mathematica/retrobop_dynamical_stability.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <list>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "base/array.hpp"
#include "base/bundle.hpp"
#include "base/get_line.hpp"
#include "base/hexadecimal.hpp"
#include "ksp_plugin/frames.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/hierarchical_system.hpp"
#include "quantities/astronomy.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using base::Bundle;
using base::not_null;
using base::Status;
using base::GetLine;
using base::HexadecimalDecode;
using base::UniqueBytes;
using geometry::BarycentreCalculator;
using geometry::Instant;
using geometry::Position;
using geometry::Sign;
using geometry::Vector;
using integrators::FixedStepSizeIntegrator;
using ksp_plugin::Barycentric;
using physics::DegreesOfFreedom;
using physics::Ephemeris;
using physics::HierarchicalSystem;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::MasslessBody;
using quantities::Angle;
using quantities::Cos;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Pow;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;
using quantities::si::Radian;
using testing_utilities::AbsoluteError;

namespace mathematica {

namespace {

enum Celestial {
  Sun,
  Moho,
  Eve,
  Gilly,
  Kerbin,
  Mun,
  Minmus,
  Duna,
  Ike,
  Dres,
  Jool,
  Laythe,
  Vall,
  Tylo,
  Bop,
  Pol,
  Eeloo,
};

constexpr std::array<Celestial, 17> celestials = {
    Sun,
    Moho,
    Eve,
    Gilly,
    Kerbin,
    Mun,
    Minmus,
    Duna,
    Ike,
    Dres,
    Jool,
    Laythe,
    Vall,
    Tylo,
    Bop,
    Pol,
    Eeloo,
};

constexpr std::array<char const*, 17> names = {
    "Sun",
    "Moho",
    "Eve",
    "Gilly",
    "Kerbin",
    "Mun",
    "Minmus",
    "Duna",
    "Ike",
    "Dres",
    "Jool",
    "Laythe",
    "Vall",
    "Tylo",
    "Bop",
    "Pol",
    "Eeloo",
};

constexpr char const* system_file = "ksp_fixed_system.proto.hex";

constexpr Instant ksp_epoch;
constexpr Instant a_century_hence = ksp_epoch + 100 * JulianYear;
constexpr Time step = 5 * Minute;

constexpr std::array<Celestial, 6> jool_system =
    {Jool, Laythe, Vall, Tylo, Bop, Pol};
constexpr std::array<Celestial, 5> jool_moons = {Laythe, Vall, Tylo, Bop, Pol};

template<typename Message>
std::unique_ptr<Message> Read(std::ifstream & file) {
  std::string const line = GetLine(file);
  if (line.empty()) {
    return nullptr;
  }
  std::uint8_t const* const hexadecimal =
      reinterpret_cast<std::uint8_t const*>(line.data());
  int const hexadecimal_size = line.size();
  UniqueBytes bytes(hexadecimal_size >> 1);
  HexadecimalDecode({hexadecimal, hexadecimal_size},
                    {bytes.data.get(), bytes.size});
  auto message = std::make_unique<Message>();
  CHECK(
      message->ParseFromArray(bytes.data.get(), static_cast<int>(bytes.size)));
  return std::move(message);
}

HierarchicalSystem<Barycentric>::BarycentricSystem ReadSystem(
    std::string const& file_name) {
  std::ifstream input_file(file_name, std::ios::in);
  CHECK(!input_file.fail());
  HierarchicalSystem<Barycentric>::BarycentricSystem system;
  for (auto body = Read<serialization::MassiveBody>(input_file);
       body != nullptr;
       body = Read<serialization::MassiveBody>(input_file)) {
    auto const degrees_of_freedom = Read<serialization::Pair>(input_file);
    CHECK(degrees_of_freedom != nullptr);
    system.bodies.push_back(MassiveBody::ReadFromMessage(*body));
    system.degrees_of_freedom.push_back(
        DegreesOfFreedom<Barycentric>::ReadFromMessage(*degrees_of_freedom));
  }
  input_file.close();
  return system;
}

Position<Barycentric> EvaluatePosition(Ephemeris<Barycentric> const& ephemeris,
                                       Celestial const celestial,
                                       Instant const& t) {
  return ephemeris.trajectory(ephemeris.bodies()[celestial])
      ->EvaluatePosition(t, /*hint=*/nullptr);
}

DegreesOfFreedom<Barycentric> EvaluateDegreesOfFreedom(
    Ephemeris<Barycentric> const& ephemeris,
    Celestial const celestial,
    Instant const& t) {
  return ephemeris.trajectory(ephemeris.bodies()[celestial])
      ->EvaluateDegreesOfFreedom(t, /*hint=*/nullptr);
}

DegreesOfFreedom<Barycentric> JoolSystemBarycentre(
    Ephemeris<Barycentric> const& ephemeris,
    Instant const& t) {
  BarycentreCalculator<DegreesOfFreedom<Barycentric>, GravitationalParameter>
      jool_system_barycentre;
  for (auto const celestial : jool_system) {
    jool_system_barycentre.Add(
        EvaluateDegreesOfFreedom(ephemeris, celestial, t),
        ephemeris.bodies()[celestial]->gravitational_parameter());
  }
  return jool_system_barycentre.Get();
}

not_null<std::unique_ptr<Ephemeris<Barycentric>>> MakeEphemeris(
    HierarchicalSystem<Barycentric>::BarycentricSystem&& system,  // NOLINT
    FixedStepSizeIntegrator<
        Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator,
    Time const& step) {
  return std::make_unique<Ephemeris<Barycentric>>(
      std::move(system.bodies),
      system.degrees_of_freedom,
      ksp_epoch,
      /*fitting_tolerance=*/1 * Milli(Metre),
      Ephemeris<Barycentric>::FixedStepParameters(integrator, step));
}

template<typename BitGenerator>
Vector<double, Barycentric> RandomUnitVector(BitGenerator& generator) {
  std::uniform_real_distribution<> longitude_distribution(-π, π);
  std::uniform_real_distribution<> z_distribution(-1, 1);
  double const z = z_distribution(generator);
  Angle const longitude = longitude_distribution(generator) * Radian;
  return Vector<double, Barycentric>({Cos(longitude) * Sqrt(1 - Pow<2>(z)),
                                      Sin(longitude) * Sqrt(1 - Pow<2>(z)),
                                      z});
}

template<typename Integrator>
std::list<not_null<std::unique_ptr<Ephemeris<Barycentric>>>>
MakePerturbedEphemerides(int const count,
                         Integrator const& integrator,
                         Time const& step) {
  std::mt19937_64 generator;
  std::list<not_null<std::unique_ptr<Ephemeris<Barycentric>>>> result;
  for (int i = 0; i < count; ++i) {
    auto system = ReadSystem(system_file);
    for (auto const celestial : jool_system) {
      system.degrees_of_freedom[celestial] = {
          system.degrees_of_freedom[celestial].position() +
              RandomUnitVector(generator) * Milli(Metre),
          system.degrees_of_freedom[celestial].velocity()};
    }
    result.emplace_back(MakeEphemeris(std::move(system), integrator, step));
  }
  return result;
}

void FillPositions(
    Ephemeris<Barycentric> const& ephemeris,
    Instant const& initial_time,
    Time const& duration,
    std::vector<std::vector<Vector<double, Barycentric>>>& container) {
  for (int n = 0; n * step < duration; ++n) {
    Instant const t = initial_time + n * step;
    Position<Barycentric> const jool_system_barycentre =
        JoolSystemBarycentre(ephemeris, t).position();
    container.emplace_back();
    for (auto const celestial : jool_system) {
      container.back().emplace_back(
          (EvaluatePosition(ephemeris, celestial, t) - jool_system_barycentre) /
          Metre);
    }
  }
}

void ProduceCenturyPlots(Ephemeris<Barycentric>& ephemeris) {
  for (int i = 1; i < (a_century_hence - ksp_epoch) / JulianYear; ++i) {
    LOG(INFO) << "year " << i;
    ephemeris.Prolong(ksp_epoch + i * JulianYear);
  }
  ephemeris.Prolong(a_century_hence);

  std::map<Celestial, std::vector<double>> extremal_separations_in_m;
  std::map<Celestial, std::vector<double>> times_in_s;
  Instant t = ksp_epoch;
  std::map<Celestial, Length> last_separations;
  std::map<Celestial, Sign> last_separation_changes;

  // Stock elements.
  std::vector<double> bop_eccentricities;
  std::vector<double> bop_inclinations_in_degrees;
  std::vector<double> bop_nodes_in_degrees;
  std::vector<double> bop_arguments_of_periapsis_in_degrees;

  // Elements around the barycentre of Jool, Laythe, and Vall.
  std::vector<double> bop_jacobi_eccentricities;
  std::vector<double> bop_jacobi_nodes_in_degrees;
  std::vector<double> bop_jacobi_inclinations_in_degrees;
  std::vector<double> bop_jacobi_arguments_of_periapsis_in_degrees;

  std::vector<double> tylo_bop_separations_in_m;
  std::vector<double> pol_bop_separations_in_m;

  std::map<Celestial, double> record_separation_in_m;

  for (auto const moon : jool_moons) {
    last_separation_changes.emplace(moon, Sign(+1));
  }
  for (int n = 0; t < a_century_hence; ++n, t = ksp_epoch + n * Hour) {
    auto const jool_position = EvaluatePosition(ephemeris, Jool, t);

    for (auto const moon : jool_moons) {
      Length const separation =
          (jool_position - EvaluatePosition(ephemeris, moon, t)).Norm();
      Sign const separation_change = Sign(separation - last_separations[moon]);
      if (separation_change != last_separation_changes.at(moon)) {
        extremal_separations_in_m[moon].emplace_back(last_separations[moon] /
                                                     Metre);
        times_in_s[moon].emplace_back((t - 1 * Hour - ksp_epoch) / Second);
        record_separation_in_m[moon] =
            std::max(record_separation_in_m[moon],
                     extremal_separations_in_m[moon].back());
        if (extremal_separations_in_m[moon].back() > 2.2e+08 &&
            extremal_separations_in_m[moon].back() ==
                record_separation_in_m[moon]) {
          LOG(WARNING) << "After "
                       << times_in_s[moon].back() * Second / JulianYear
                       << " a, " << names[moon] << " has an apsis at "
                       << extremal_separations_in_m[moon].back() * Metre;
        }
      }
      last_separations[moon] = separation;
      last_separation_changes.at(moon) = separation_change;
    }

    tylo_bop_separations_in_m.emplace_back(
        (EvaluatePosition(ephemeris, Tylo, t) -
         EvaluatePosition(ephemeris, Bop, t)).Norm() / Metre);
    pol_bop_separations_in_m.emplace_back(
        (EvaluatePosition(ephemeris, Pol, t) -
         EvaluatePosition(ephemeris, Bop, t)).Norm() / Metre);

    {
      // KSP's osculating elements.
      auto const bop_elements =
          KeplerOrbit<Barycentric>(
              *ephemeris.bodies()[Jool],
              MasslessBody(),
              EvaluateDegreesOfFreedom(ephemeris, Bop, t) -
                  EvaluateDegreesOfFreedom(ephemeris, Jool, t), t)
              .elements_at_epoch();
      bop_eccentricities.emplace_back(bop_elements.eccentricity);
      bop_inclinations_in_degrees.emplace_back(bop_elements.inclination /
                                               Degree);
      bop_nodes_in_degrees.emplace_back(
          bop_elements.longitude_of_ascending_node / Degree);
      bop_arguments_of_periapsis_in_degrees.emplace_back(
          bop_elements.argument_of_periapsis / Degree);
    }

    {
      BarycentreCalculator<DegreesOfFreedom<Barycentric>,
                           GravitationalParameter>
          innermost_jool_system;
      for (auto const celestial : {Jool, Laythe, Vall, Tylo}) {
        innermost_jool_system.Add(
            EvaluateDegreesOfFreedom(ephemeris, celestial, t),
            ephemeris.bodies()[celestial]->gravitational_parameter());
      }
      auto const bop_jacobi_elements =
          KeplerOrbit<Barycentric>(MassiveBody(innermost_jool_system.weight()),
                                   *ephemeris.bodies()[Bop],
                                   EvaluateDegreesOfFreedom(ephemeris, Bop, t) -
                                       innermost_jool_system.Get(),
                                   t).elements_at_epoch();
      bop_jacobi_eccentricities.emplace_back(bop_jacobi_elements.eccentricity);
      bop_jacobi_inclinations_in_degrees.emplace_back(
          bop_jacobi_elements.inclination / Degree);
      bop_jacobi_nodes_in_degrees.emplace_back(
          bop_jacobi_elements.longitude_of_ascending_node / Degree);
      bop_jacobi_arguments_of_periapsis_in_degrees.emplace_back(
          bop_jacobi_elements.argument_of_periapsis / Degree);
    }
  }

  std::ofstream file;
  file.open("retrobop_century.generated.wl");
  file << mathematica::Assign("laytheTimes", times_in_s[Laythe]);
  file << mathematica::Assign("vallTimes", times_in_s[Vall]);
  file << mathematica::Assign("tyloTimes", times_in_s[Tylo]);
  file << mathematica::Assign("polTimes", times_in_s[Pol]);
  file << mathematica::Assign("bopTimes", times_in_s[Bop]);
  file << mathematica::Assign("laytheSeparations",
                              extremal_separations_in_m[Laythe]);
  file << mathematica::Assign("vallSeparations",
                              extremal_separations_in_m[Vall]);
  file << mathematica::Assign("tyloSeparations",
                              extremal_separations_in_m[Tylo]);
  file << mathematica::Assign("polSeparations", extremal_separations_in_m[Pol]);
  file << mathematica::Assign("bopSeparations", extremal_separations_in_m[Bop]);

  file << mathematica::Assign("bopEccentricities", bop_eccentricities);
  file << mathematica::Assign("bopInclinations", bop_inclinations_in_degrees);
  file << mathematica::Assign("bopNodes", bop_nodes_in_degrees);
  file << mathematica::Assign("bopArguments",
                              bop_arguments_of_periapsis_in_degrees);
  file << mathematica::Assign("bopJacobiEccentricities",
                              bop_jacobi_eccentricities);
  file << mathematica::Assign("bopJacobiInclinations",
                              bop_jacobi_inclinations_in_degrees);
  file << mathematica::Assign("bopJacobiNodes", bop_jacobi_nodes_in_degrees);
  file << mathematica::Assign("bopJacobiArguments",
                              bop_jacobi_arguments_of_periapsis_in_degrees);

  file << mathematica::Assign("tyloBop", tylo_bop_separations_in_m);
  file << mathematica::Assign("polBop", pol_bop_separations_in_m);
  file.close();
}

void ComputeHighestMoonError(Ephemeris<Barycentric> const& left,
                             Ephemeris<Barycentric> const& right,
                             Instant const& t,
                             Length& error,
                             Celestial& most_erroneous_moon) {
  error = Length{};
  auto const left_barycentre = JoolSystemBarycentre(left, t).position();
  auto const right_barycentre = JoolSystemBarycentre(right, t).position();
  for (auto const moon : jool_moons) {
    Length const moon_error =
        AbsoluteError(EvaluatePosition(left, moon, t) - left_barycentre,
                      EvaluatePosition(right, moon, t) - right_barycentre);
    if (moon_error > error) {
      error = moon_error;
      most_erroneous_moon = moon;
    }
  }
}

}  // namespace

void PlotPredictableYears() {
  auto const ephemeris = MakeEphemeris(
      ReadSystem(system_file),
      integrators::QuinlanTremaine1990Order12<Position<Barycentric>>(),
      step);

  for (int i = 1; i <= 5; ++i) {
    ephemeris->Prolong(ksp_epoch + i * JulianYear);
    LOG(INFO) << "Prolonged to year " << i;
  }

  std::vector<std::vector<Vector<double, Barycentric>>>
      barycentric_positions_1_year;
  FillPositions(
      *ephemeris, ksp_epoch, 1 * JulianYear, barycentric_positions_1_year);
  std::vector<std::vector<Vector<double, Barycentric>>>
      barycentric_positions_2_year;
  FillPositions(
      *ephemeris, ksp_epoch, 2 * JulianYear, barycentric_positions_2_year);
  std::vector<std::vector<Vector<double, Barycentric>>>
      barycentric_positions_5_year;
  FillPositions(
      *ephemeris, ksp_epoch, 5 * JulianYear, barycentric_positions_5_year);

  std::ofstream file;
  file.open("retrobop_predictable_years.generated.wl");
  file << mathematica::Assign("barycentricPositions1",
                              barycentric_positions_1_year);
  file << mathematica::Assign("barycentricPositions2",
                              barycentric_positions_2_year);
  file << mathematica::Assign("barycentricPositions5",
                              barycentric_positions_5_year);
  file.close();
}

void PlotCentury() {
  ProduceCenturyPlots(*MakeEphemeris(
      ReadSystem(system_file),
      integrators::QuinlanTremaine1990Order12<Position<Barycentric>>(),
      step));
}

void AnalyseGlobalError() {
  auto const reference_ephemeris = MakeEphemeris(
      ReadSystem(system_file),
      integrators::QuinlanTremaine1990Order12<Position<Barycentric>>(),
      step);
  std::unique_ptr<Ephemeris<Barycentric>> refined_ephemeris = MakeEphemeris(
      ReadSystem(system_file),
      integrators::QuinlanTremaine1990Order12<Position<Barycentric>>(),
      step / 2);
  std::list<not_null<std::unique_ptr<Ephemeris<Barycentric>>>>
      perturbed_ephemerides = MakePerturbedEphemerides(
          100,
          integrators::QuinlanTremaine1990Order12<Position<Barycentric>>(),
          step);

  bool log_radius = true;
  // Errors below this are invisible on plots.
  Length const visible_threshold = 1e6 * Metre;
  // Errors above this mean we are pretty much completely out of phase.
  Length const chaotic_threshold = 1e8 * Metre;

  for (int year = 1;; ++year) {
    Instant const t = ksp_epoch + year * JulianYear;
    Bundle bundle{7};
    if (reference_ephemeris != nullptr) {
      bundle.Add([&reference_ephemeris = *reference_ephemeris, t ]() {
        reference_ephemeris.Prolong(t);
        reference_ephemeris.ForgetBefore(t);
        return Status::OK;
      });
    }
    if (refined_ephemeris != nullptr) {
      bundle.Add([&refined_ephemeris = *refined_ephemeris, t ]() {
        refined_ephemeris.Prolong(t);
        refined_ephemeris.ForgetBefore(t);
        return Status::OK;
      });
    }
    for (auto const& ephemeris : perturbed_ephemerides) {
      bundle.Add([ephemeris = ephemeris.get(), t]() {
        ephemeris->Prolong(t);
        ephemeris->ForgetBefore(t);
        return Status::OK;
      });
    }
    bundle.Join();
    LOG(INFO) << "year " << year;

    if (refined_ephemeris != nullptr) {
      Length numerical_error;
      Celestial most_erroneous_moon;
      ComputeHighestMoonError(*refined_ephemeris,
                              *reference_ephemeris,
                              t,
                              numerical_error,
                              most_erroneous_moon);
      LOG(INFO) << "Numerical error: " << numerical_error << " ("
                << names[most_erroneous_moon] << ")";
      LOG_IF(INFO, numerical_error < visible_threshold) << "invisible on plots";
      if (numerical_error > chaotic_threshold) {
        LOG(INFO) << u8"The wrath of Ляпунов is upon us!";
        refined_ephemeris.reset();
      }
    }

    if (log_radius) {
      Length cluster_radius;
      Celestial most_erroneous_moon;
      for (auto const& ephemeris : perturbed_ephemerides) {
        Length moon_error;
        Celestial moon;
        ComputeHighestMoonError(*ephemeris,
                                *reference_ephemeris,
                                t,
                                moon_error,
                                moon);
        if (moon_error > cluster_radius) {
          cluster_radius = moon_error;
          most_erroneous_moon = moon;
        }
      }
      LOG(INFO) << "Cluster radius: " << cluster_radius << " ("
                << names[most_erroneous_moon] << ")";
      LOG_IF(INFO, cluster_radius < visible_threshold) << "invisible on plots";
      if (cluster_radius > chaotic_threshold) {
        LOG(INFO) << u8"The wrath of Ляпунов is upon us!";
        log_radius = false;
      }
    }

    if (!log_radius && refined_ephemeris == nullptr) {
      return;
    }
  }
}

void StatisticallyAnalyseStability() {
  std::list<not_null<std::unique_ptr<Ephemeris<Barycentric>>>>
      perturbed_ephemerides = MakePerturbedEphemerides(
          100,
          integrators::QuinlanTremaine1990Order12<Position<Barycentric>>(),
          step);

  std::map<not_null<Ephemeris<Barycentric>*>, bool> numerically_unsound;
  for (auto const& ephemeris : perturbed_ephemerides) {
    ephemeris->Prolong(ksp_epoch);
    numerically_unsound[ephemeris.get()] = false;
  }
  int total_breakdowns = 0;
  // If the error between an integration at step and one at step/2 exceeds this
  // over a year, we assume that things have happened that our integrator cannot
  // handle, probably close encounters.
  // TODO(egg): This is a very lousy substitute for a proper estimation of the
  // local forward error.  We probably want to have a way to actually estimate
  // the local error (on every step), and perhaps even the local backward error
  // (though that may be costly if done naïvely).
  Length const yearly_allowed_numerical_error = 1000 * Metre;

  for (int year = 1; year <= 200; ++year) {
    Instant const t = ksp_epoch + year * JulianYear;
    Bundle bundle{7};
    for (auto const& ephemeris : perturbed_ephemerides) {
      bundle.Add([
        &numerically_unsound,
        ephemeris = ephemeris.get(),
        t,
        yearly_allowed_numerical_error
      ]() {
        auto system = ReadSystem(system_file);
        for (auto const celestial : celestials) {
          system.degrees_of_freedom[celestial] =
              EvaluateDegreesOfFreedom(*ephemeris,
                                       celestial,
                                       ephemeris->t_min());
        }
        Ephemeris<Barycentric> refined(
            std::move(system.bodies),
            system.degrees_of_freedom,
            ephemeris->t_min(),
            1 * Milli(Metre),
            Ephemeris<Barycentric>::FixedStepParameters(
                integrators::QuinlanTremaine1990Order12<
                    Position<Barycentric>>(),
                step / 2));
        ephemeris->Prolong(t);
        ephemeris->ForgetBefore(t);
        refined.Prolong(t);
        Length numerical_error;
        Celestial most_erroneous_moon;
        ComputeHighestMoonError(refined,
                                *ephemeris,
                                t,
                                numerical_error,
                                most_erroneous_moon);
        if (numerical_error > yearly_allowed_numerical_error) {
          LOG(INFO) << "high numerical error " << numerical_error << " ("
                    << names[most_erroneous_moon] << ")";
          numerically_unsound[ephemeris] = true;
        }
        return Status::OK;
      });
    }
    bundle.Join();
    LOG(INFO) << "year " << year;

    int yearly_breakdowns = 0;
    for (auto it = perturbed_ephemerides.begin();
         it != perturbed_ephemerides.end();) {
      if (numerically_unsound[it->get()]) {
        goto erase_and_skip_to_next_perturbed_ephemeris;
      } else {
        Ephemeris<Barycentric> const& ephemeris = **it;
        auto const jool_barycentre =
            JoolSystemBarycentre(ephemeris, t).position();
        for (auto const moon : jool_moons) {
          Length const distance =
              (EvaluatePosition(ephemeris, moon, t) - jool_barycentre).Norm();
          if (distance > 3e8 * Metre) {
            LOG(INFO) << names[moon] << " escape, " << distance
                      << " from Jool.";
            ++yearly_breakdowns;
            ++total_breakdowns;
            goto erase_and_skip_to_next_perturbed_ephemeris;
          }
        }
      }
      ++it;
      continue;
    erase_and_skip_to_next_perturbed_ephemeris:
      perturbed_ephemerides.erase(it++);
    }
    LOG(INFO) << "cluster size is " << perturbed_ephemerides.size();
    LOG_IF(INFO, yearly_breakdowns > 0) << yearly_breakdowns << " breakdowns";
    LOG_IF(INFO, total_breakdowns > 0) << total_breakdowns << " thus far";
  }
}

}  // namespace mathematica
}  // namespace principia
