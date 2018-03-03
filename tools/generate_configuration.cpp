
#include "tools/generate_configuration.hpp"

#include <experimental/filesystem>
#include <iomanip>
#include <limits>
#include <string>

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "physics/solar_system.hpp"
#include "quantities/parser.hpp"
#include "quantities/si.hpp"
#include "serialization/astronomy.pb.h"

namespace principia {

using astronomy::ICRFJ2000Equator;
using astronomy::J2000;
using astronomy::JulianDate;
using physics::SolarSystem;
using quantities::si::Second;

namespace {
constexpr char cfg[] = "cfg";
constexpr char proto_txt[] = "proto.txt";
}  // namespace

namespace tools {
namespace internal_generate_configuration {

void GenerateConfiguration(std::string const& game_epoch,
                           std::string const& gravity_model_stem,
                           std::string const& initial_state_stem,
                           std::string const& numerics_blueprint_stem) {
  std::experimental::filesystem::path const directory =
      SOLUTION_DIR / "astronomy";
  SolarSystem<ICRFJ2000Equator> solar_system(
      (directory / gravity_model_stem).replace_extension(proto_txt),
      (directory / initial_state_stem).replace_extension(proto_txt));

  std::ofstream gravity_model_cfg(
      (directory / gravity_model_stem).replace_extension(cfg));
  CHECK(gravity_model_cfg.good());
  gravity_model_cfg << "principia_gravity_model:NEEDS[RealSolarSystem] {\n";
  for (std::string const& name : solar_system.names()) {
    serialization::GravityModel::Body const& body =
        solar_system.gravity_model_message(name);
    gravity_model_cfg << "  body {\n";
    gravity_model_cfg << "    name                    = "
                      << name << "\n";
    gravity_model_cfg << "    gravitational_parameter = "
                      << body.gravitational_parameter() << "\n";
    if (body.has_reference_instant()) {
      gravity_model_cfg << "    reference_instant       = "
                        << body.reference_instant() << "\n";
    }
    if (body.has_mean_radius()) {
      gravity_model_cfg << "    mean_radius             = "
                        << body.mean_radius() << "\n";
    }
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
    if (body.has_j2()) {
      gravity_model_cfg << "    j2                      = "
                        << std::scientific
                        << std::setprecision(
                               std::numeric_limits<double>::max_digits10)
                        << body.j2() << "\n";
    }
    if (body.has_reference_radius()) {
      gravity_model_cfg << "    reference_radius        = "
                        << body.reference_radius() << "\n";
    }
    gravity_model_cfg << "  }\n";
  }
  gravity_model_cfg << "}\n";

  std::ofstream initial_state_cfg(
      (directory / initial_state_stem).replace_extension(cfg));
  CHECK(initial_state_cfg.good());
  initial_state_cfg << "principia_initial_state:NEEDS[RealSolarSystem] {\n";
  initial_state_cfg << "  game_epoch = " << game_epoch << "\n";
  initial_state_cfg << "  solar_system_epoch = "
                    << solar_system.epoch_literal() << "\n";
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
  initial_state_cfg << "}\n";

  // Parse the numerics blueprint file here, it doesn't belong in class
  // SolarSystem.
  std::experimental::filesystem::path numerics_blueprint_filename =
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
      "principia_numerics_blueprint {\n";

  numerics_blueprint_cfg << "  ephemeris {\n";
  numerics_blueprint_cfg
      << "    fixed_step_size_integrator = "
      << serialization::FixedStepSizeIntegrator::Kind_Name(
            ephemeris.integrator()) << "\n";
  numerics_blueprint_cfg << "    integration_step_size = " << ephemeris.step()
                         << "\n";
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

}  // namespace internal_generate_configuration
}  // namespace tools
}  // namespace principia
