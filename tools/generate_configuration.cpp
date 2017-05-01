
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

void GenerateConfiguration(Instant const& game_epoch,
                           std::string const& gravity_model_stem,
                           std::string const& initial_state_stem) {
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
      gravity_model_cfg << "    reference_instant       = " << std::fixed
                        << std::setprecision(
                               std::numeric_limits<long double>::digits10 + 1)
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
  initial_state_cfg << "  game_epoch = "
                    << game_epoch - J2000 << "\n";
  initial_state_cfg << "  solar_system_epoch = "
                    << solar_system.epoch() - J2000 << "\n";
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
}

}  // namespace internal_generate_configuration
}  // namespace tools
}  // namespace principia
