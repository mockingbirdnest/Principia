#include "tools/generate_configuration.hpp"

#include <filesystem>
#include <iomanip>
#include <limits>
#include <string>

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "base/fingerprint2011.hpp"
#include "base/serialization.hpp"
#include "glog/logging.h"
#include "numerics/elementary_functions.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/solar_system.hpp"
#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/parser.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/astronomy.pb.h"

namespace principia {
namespace tools {
namespace _generate_configuration {
namespace internal {

using namespace principia::astronomy::_epoch;
using namespace principia::astronomy::_frames;
using namespace principia::base::_fingerprint2011;
using namespace principia::base::_serialization;
using namespace principia::numerics::_elementary_functions;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_constants;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_parser;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

namespace {
constexpr char cfg[] = "cfg";
constexpr char proto_txt[] = "proto.txt";
}  // namespace

std::string NormalizeLength(std::string const& s) {
  // If the string contains an R, it's expressed using astronomical radii.
  // Convert to kilometers.
  if (s.find('R') == std::string::npos) {
    return s;
  } else {
    Length const length = ParseQuantity<Length>(s);
    std::ostringstream stream;
    stream << std::scientific
           << std::setprecision(std::numeric_limits<double>::max_digits10)
           << length / Kilo(Metre) << " km";
    return stream.str();
  }
}

void GenerateConfiguration(std::string const& game_epoch,
                           std::string const& gravity_model_stem,
                           std::string const& initial_state_stem,
                           std::string const& numerics_blueprint_stem,
                           std::string const& needs) {
  std::filesystem::path const directory =
      SOLUTION_DIR / "astronomy";
  SolarSystem<ICRS> solar_system(
      (directory / gravity_model_stem).replace_extension(proto_txt),
      (directory / initial_state_stem).replace_extension(proto_txt),
      /*ignore_frame=*/true);

  std::ofstream gravity_model_cfg(
      (directory / gravity_model_stem).replace_extension(cfg));
  CHECK(gravity_model_cfg.good());
  gravity_model_cfg << "principia_gravity_model:NEEDS[" << needs << "] {\n";

  // Find the star.  This is needed to construct Kepler orbits below.
  std::optional<serialization::GravityModel::Body> star;
  for (std::string const& name : solar_system.names()) {
    serialization::GravityModel::Body const& body =
        solar_system.gravity_model_message(name);
    if (solar_system.has_keplerian_initial_state_message(name)) {
      bool const is_star =
          !solar_system.keplerian_initial_state_message(name).has_parent();
      if (is_star) {
        star = body;
      }
    }
  }

  for (std::string const& name : solar_system.names()) {
    serialization::GravityModel::Body const& body =
        solar_system.gravity_model_message(name);
    LOG(INFO) << "Fingerprint " << std::setw(16) << std::hex << std::uppercase
              << Fingerprint2011(SerializeAsBytes(body).get())
              << " for " << name;
    gravity_model_cfg << "  body {\n";
    gravity_model_cfg << "    name                    = "
                      << (star.has_value() && name == star->name() ? "Sun"
                                                                   : name)
                      << "\n";
    if (body.has_gravitational_parameter()) {
      gravity_model_cfg << "    gravitational_parameter = "
                        << body.gravitational_parameter() << "\n";
    } else {
      // If the mass was provided, convert to a gravitational parameter.
      Mass const mass = ParseQuantity<Mass>(body.mass());
      GravitationalParameter const gravitational_parameter =
          GravitationalConstant * mass;
      gravity_model_cfg << "    gravitational_parameter = "
                        << std::scientific
                        << std::setprecision(
                               std::numeric_limits<double>::max_digits10)
                        << gravitational_parameter /
                               (Pow<3>(Kilo(Metre)) / Pow<2>(Second))
                        << " km^3/s^2\n";
    }
    if (body.has_reference_instant()) {
      gravity_model_cfg << "    reference_instant       = "
                        << body.reference_instant() << "\n";
    }
    // The fields min_radius, mean_radius and max_radius come from the game and
    // are not copied from the proto to the configuration.
    if (body.has_axis_right_ascension()) {
      gravity_model_cfg << "    axis_right_ascension    = "
                        << body.axis_right_ascension() << "\n";
    }
    if (body.has_axis_declination()) {
      gravity_model_cfg << "    axis_declination        = "
                        << body.axis_declination() << "\n";
    }
    if (body.has_reference_angle()) {
      gravity_model_cfg << "    reference_angle         = "
                        << body.reference_angle() << "\n";
    }
    if (body.has_angular_frequency()) {
      gravity_model_cfg << "    angular_frequency       = "
                        << body.angular_frequency() << "\n";
    }
    if (body.has_reference_radius()) {
      gravity_model_cfg << "    reference_radius        = "
                        << NormalizeLength(body.reference_radius()) << "\n";
    }
    switch (body.oblateness_case()) {
      case serialization::GravityModel::Body::kJ2:
        gravity_model_cfg << "    j2                      = "
                          << std::scientific
                          << std::setprecision(
                                 std::numeric_limits<double>::max_digits10)
                          << body.j2() << "\n";
        break;
      case serialization::GravityModel::Body::kGeopotential:
        for (auto const& row : body.geopotential().row()) {
          gravity_model_cfg << "    geopotential_row {\n"
                            << "      degree = " << row.degree() << "\n";
          for (auto const& column : row.column()) {
            gravity_model_cfg << "      geopotential_column {\n"
                              << "        order = " << column.order() << "\n"
                              << std::scientific
                              << std::setprecision(
                                     std::numeric_limits<double>::max_digits10);
            if (column.has_j()) {
              gravity_model_cfg << "        j     = " << column.j() << "\n";
            }
            if (column.has_cos()) {
              gravity_model_cfg << "        cos   = " << column.cos() << "\n";
            }
            gravity_model_cfg << "        sin   = " << column.sin() << "\n"
                              << "      }\n";
          }
          gravity_model_cfg << "    }\n";
        }
        break;
      case serialization::GravityModel::Body::OBLATENESS_NOT_SET:
        break;
    }
    gravity_model_cfg << "  }\n";
  }
  gravity_model_cfg << "}\n";

  std::ofstream initial_state_cfg(
      (directory / initial_state_stem).replace_extension(cfg));
  CHECK(initial_state_cfg.good());
  initial_state_cfg << "principia_initial_state:NEEDS[" << needs << "] {\n";
  initial_state_cfg << "  game_epoch = " << game_epoch << "\n";
  initial_state_cfg << "  solar_system_epoch = "
                    << solar_system.epoch_literal() << "\n";

  if (solar_system.has_keplerian_initial_state_message(
          solar_system.names().front())) {
    // If the configuration is given in keplerian elements, convert it to
    // cartesian coordinates here.  This avoids depending on conversions done
    // by the game.
    auto hierarchical_system = solar_system.MakeHierarchicalSystem();
    auto const barycentric_system =
        hierarchical_system->ConsumeBarycentricSystem();

    auto displacement = [](DegreesOfFreedom<ICRS> const& dof) {
      return (dof.position() - ICRS::origin).coordinates();
    };
    auto velocity = [](DegreesOfFreedom<ICRS> const& dof) {
      return dof.velocity().coordinates();
    };

    for (int i = 0; i < barycentric_system.bodies.size(); ++i) {
      auto const& body = barycentric_system.bodies[i];
      auto const& dof = barycentric_system.degrees_of_freedom[i];
      initial_state_cfg << "  body {\n";
      initial_state_cfg << "    name = "
                        << (star.has_value() && body->name() == star->name()
                                ? "Sun"
                                : body->name())
                        << "\n";
      initial_state_cfg << "    x    = " << displacement(dof).x << "\n";
      initial_state_cfg << "    y    = " << displacement(dof).y << "\n";
      initial_state_cfg << "    z    = " << displacement(dof).z << "\n";
      initial_state_cfg << "    vx   = " << velocity(dof).x << "\n";
      initial_state_cfg << "    vy   = " << velocity(dof).y << "\n";
      initial_state_cfg << "    vz   = " << velocity(dof).z << "\n";
      initial_state_cfg << "  }\n";
    }
  } else {
    for (std::string const& name : solar_system.names()) {
      serialization::InitialState::Cartesian::Body const& body =
          solar_system.cartesian_initial_state_message(name);
      initial_state_cfg << "  body {\n";
      initial_state_cfg << "    name = " << name << "\n";
      initial_state_cfg << "    x    = " << body.x() << "\n";
      initial_state_cfg << "    y    = " << body.y() << "\n";
      initial_state_cfg << "    z    = " << body.z() << "\n";
      initial_state_cfg << "    vx   = " << body.vx() << "\n";
      initial_state_cfg << "    vy   = " << body.vy() << "\n";
      initial_state_cfg << "    vz   = " << body.vz() << "\n";
      initial_state_cfg << "  }\n";
    }
  }
  initial_state_cfg << "}\n";

  // Parse the numerics blueprint file here, it doesn't belong in class
  // SolarSystem.
  std::filesystem::path numerics_blueprint_filename =
      (directory / numerics_blueprint_stem).replace_extension(proto_txt);
  serialization::SolarSystemFile numerics_blueprint;
  std::ifstream numerics_blueprint_ifstream(numerics_blueprint_filename);
  CHECK(numerics_blueprint_ifstream.good());
  google::protobuf::io::IstreamInputStream numerics_blueprint_zcs(
                                               &numerics_blueprint_ifstream);
  CHECK(google::protobuf::TextFormat::Parse(&numerics_blueprint_zcs,
                                            &numerics_blueprint));
  CHECK(numerics_blueprint.has_numerics_blueprint());
  auto const& ephemeris = numerics_blueprint.numerics_blueprint().ephemeris();
  auto const& history = numerics_blueprint.numerics_blueprint().history();
  auto const& psychohistory =
      numerics_blueprint.numerics_blueprint().psychohistory();

  std::ofstream numerics_blueprint_cfg(
      (directory / numerics_blueprint_stem).replace_extension(cfg));
  CHECK(numerics_blueprint_cfg.good());
  numerics_blueprint_cfg <<
      "principia_numerics_blueprint:NEEDS[" << needs << "] {\n";

  if (numerics_blueprint.numerics_blueprint().has_downsampling()) {
    auto const& downsampling =
        numerics_blueprint.numerics_blueprint().downsampling();
    numerics_blueprint_cfg << "  downsampling {\n";
    numerics_blueprint_cfg << "    max_dense_intervals = "
                           << downsampling.max_dense_intervals() << "\n";
    numerics_blueprint_cfg << "    tolerance = "
                           << downsampling.tolerance() << "\n";
    numerics_blueprint_cfg << "  }\n";
  }

  numerics_blueprint_cfg << "  ephemeris {\n";
  numerics_blueprint_cfg
      << "    fixed_step_size_integrator = "
      << serialization::FixedStepSizeIntegrator::Kind_Name(
            ephemeris.integrator()) << "\n";
  numerics_blueprint_cfg << "    integration_step_size = " << ephemeris.step()
                         << "\n";
  numerics_blueprint_cfg << "    fitting_tolerance = "
      << ephemeris.fitting_tolerance() << "\n";
  numerics_blueprint_cfg << "    geopotential_tolerance = "
      << ephemeris.geopotential_tolerance() << "\n";
  numerics_blueprint_cfg << "  }\n";

  numerics_blueprint_cfg << "  history {\n";
  numerics_blueprint_cfg
      << "    fixed_step_size_integrator = "
      << serialization::FixedStepSizeIntegrator::Kind_Name(
            history.integrator()) << "\n";
  numerics_blueprint_cfg << "    integration_step_size = " << history.step()
                         << "\n";
  numerics_blueprint_cfg << "  }\n";

  numerics_blueprint_cfg << "  psychohistory {\n";
  numerics_blueprint_cfg
      << "    adaptive_step_size_integrator = "
      << serialization::AdaptiveStepSizeIntegrator::Kind_Name(
             psychohistory.integrator()) << "\n";
  numerics_blueprint_cfg << "    length_integration_tolerance = "
                         << psychohistory.length_integration_tolerance()
                         << "\n";
  numerics_blueprint_cfg << "    speed_integration_tolerance = "
                         << psychohistory.speed_integration_tolerance()
                         << "\n";
  numerics_blueprint_cfg << "  }\n";

  numerics_blueprint_cfg << "}\n";
}

}  // namespace internal
}  // namespace _generate_configuration
}  // namespace tools
}  // namespace principia
