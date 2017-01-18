
#include "mathematica/retrobop_dynamical_stability.hpp"

#include <array>
#include <fstream>
#include <memory>

#include "base/array.hpp"
#include "base/get_line.hpp"
#include "base/hexadecimal.hpp"
#include "ksp_plugin/frames.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/hierarchical_system.hpp"
#include "quantities/astronomy.hpp"

namespace principia {

using base::GetLine;
using base::HexadecimalDecode;
using base::UniqueBytes;
using geometry::BarycentreCalculator;
using geometry::Instant;
using geometry::Position;
using geometry::Sign;
using geometry::Vector;
using ksp_plugin::Barycentric;
using physics::DegreesOfFreedom;
using physics::Ephemeris;
using physics::HierarchicalSystem;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::MasslessBody;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;

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

constexpr Instant ksp_epoch;
constexpr Instant a_century_hence = ksp_epoch + 1000 * JulianYear;
constexpr Time step = 8 * Minute;

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

void FillPositions(
    Ephemeris<Barycentric> const& ephemeris,
    Instant const& initial_time,
    Time const& duration,
    std::vector<std::vector<Vector<double, Barycentric>>>& container) {
  for (Instant t = initial_time; t < initial_time + duration;
       t += step) {
    auto const position = [t, &ephemeris](Celestial celestial) {
      return ephemeris.trajectory(ephemeris.bodies()[celestial])
          ->EvaluatePosition(t, /*hint=*/nullptr);
    };
    BarycentreCalculator<Position<Barycentric>, GravitationalParameter>
        jool_system_barycentre;
    for (auto const celestial : jool_system) {
      jool_system_barycentre.Add(
          position(celestial),
          ephemeris.bodies()[celestial]->gravitational_parameter());
    }
    container.emplace_back();
    for (auto const celestial : jool_system) {
      container.back().emplace_back(
          (position(celestial) - jool_system_barycentre.Get()) / Metre);
    }
  }
}

}  // namespace

void SimulateFixedSystem(bool produce_file) {
  std::ifstream input_file("ksp_fixed_system.proto.hex", std::ios::in);
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
    LOG(ERROR) << system.bodies.back()->gravitational_parameter();
    LOG(ERROR) << system.degrees_of_freedom.back();
  }
  input_file.close();
  Ephemeris<Barycentric> ephemeris(
      std::move(system.bodies),
      system.degrees_of_freedom,
      ksp_epoch,
      /*fitting_tolerance=*/1 * Milli(Metre),
      Ephemeris<Barycentric>::FixedStepParameters(
          integrators::QuinlanTremaine1990Order12<Position<Barycentric>>(),
          step));
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
    auto const position = [t, &ephemeris](Celestial const celestial) {
      return ephemeris.trajectory(ephemeris.bodies()[celestial])
          ->EvaluatePosition(t, /*hint=*/nullptr);
    };
    auto const degrees_of_freedom = [t, &ephemeris](Celestial const celestial) {
      return ephemeris.trajectory(ephemeris.bodies()[celestial])
          ->EvaluateDegreesOfFreedom(t, /*hint=*/nullptr);
    };
    auto const jool_position = position(Jool);

    for (auto const moon : jool_moons) {
      Length const separation = (jool_position - position(moon)).Norm();
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
        (position(Tylo) - position(Bop)).Norm() / Metre);
    pol_bop_separations_in_m.emplace_back(
        (position(Pol) - position(Bop)).Norm() / Metre);

    {
      // KSP's osculating elements.
      auto const bop_elements =
          KeplerOrbit<Barycentric>(
              *ephemeris.bodies()[Jool],
              MasslessBody(),
              degrees_of_freedom(Bop) - degrees_of_freedom(Jool),
              t).elements_at_epoch();
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
            degrees_of_freedom(celestial),
            ephemeris.bodies()[celestial]->gravitational_parameter());
      }
      auto const bop_jacobi_elements =
          KeplerOrbit<Barycentric>(
              MassiveBody(innermost_jool_system.weight()),
              *ephemeris.bodies()[Bop],
              degrees_of_freedom(Bop) - innermost_jool_system.Get(),
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

  std::vector<std::vector<Vector<double, Barycentric>>>
      barycentric_positions_1_year;
  FillPositions(ephemeris,
                ksp_epoch,
                1 * JulianYear,
                barycentric_positions_1_year);
  std::vector<std::vector<Vector<double, Barycentric>>>
      barycentric_positions_2_year;
  FillPositions(ephemeris,
                ksp_epoch,
                2 * JulianYear,
                barycentric_positions_2_year);

  if (!produce_file) {
    return;
  }

  std::ofstream file;
  file.open("ksp_system.generated.wl");
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

  file << mathematica::Assign("barycentricPositions1",
                              barycentric_positions_1_year);
  file << mathematica::Assign("barycentricPositions2",
                              barycentric_positions_2_year);
  file.close();
}

}  // namespace mathematica
}  // namespace principia
