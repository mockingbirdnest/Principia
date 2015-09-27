#include <string>

#include <experimental/filesystem>

#include "astronomy/frames.hpp"
#include "glog/logging.h"
#include "physics/solar_system.hpp"
#include "serialization/astronomy.pb.h"

namespace principia {

using physics::SolarSystem;

namespace {
constexpr char kCfg[] = "cfg";
constexpr char kProtoTxt[] = "proto.txt";
}  // namespace

namespace astronomy {

void GenerateConfiguration(std::string const& gravity_model_stem,
                           std::string const& initial_state_stem) {
  std::experimental::filesystem::path const directory =
      SOLUTION_DIR / "astronomy";
  SolarSystem<ICRFJ2000Equator> solar_system;
  solar_system.Initialize(
      (directory / gravity_model_stem).replace_extension(kProtoTxt),
      (directory / initial_state_stem).replace_extension(kProtoTxt));

  std::ofstream gravity_model_cfg(
      (directory / gravity_model_stem).replace_extension(kCfg));
  gravity_model_cfg << "principia_gravity_model:NEEDS["
                    << solar_system.gravity_model_needs()
                    << "] {\n";
  for (std::string const& name : solar_system.names()) {
    serialization::GravityModel::Body const& body =
        solar_system.gravity_model_message(name);
    gravity_model_cfg << "  body {\n";
    gravity_model_cfg << "    name                    = "
                      << name << "\n";
    gravity_model_cfg << "    gravitational_parameter = " 
                      << body.gravitational_parameter() << "\n";
    if (body.has_axis_right_ascension()) {
      gravity_model_cfg << "    axis_right_ascension    = "
                        << body.axis_right_ascension() << "\n";
    }
    if (body.has_axis_declination()) {
      gravity_model_cfg << "    axis_declination        = "
                        << body.axis_declination() << "\n";
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
      (directory / initial_state_stem).replace_extension(kCfg));

}

}  // namespace astronomy
}  // namespace principia

int main(int argc, char const* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::LogToStderr();
  if (argc < 3) {
    std::cerr << "Usage: generate_configuration "
                 "gravity_model_stem initial_state_stem";
    return 1;
  }
  std::string const gravity_model_stem = argv[1];
  std::string const initial_state_stem = argv[2];
  principia::astronomy::GenerateConfiguration(gravity_model_stem,
                                              initial_state_stem);
  return 0;
}
