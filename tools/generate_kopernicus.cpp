
#include "tools/generate_kopernicus.hpp"

#include <filesystem>
#include <map>

#include "astronomy/frames.hpp"
#include "base/map_util.hpp"
#include "physics/solar_system.hpp"
#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/parser.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using base::FindOrDie;
using physics::SolarSystem;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::Mass;
using quantities::ParseQuantity;
using quantities::constants::GravitationalConstant;

namespace tools {

namespace {

constexpr char cfg[] = "cfg";
constexpr char proto_txt[] = "proto.txt";

std::map<std::string, std::string> const body_description_map = {
    {"Trappist-1",
     "An ultra-cool red dwarf star of spectral type M8 V, 40 light-years from "
     "Earth; also known as 2MUCD 12171 and 2MASS J23062928-0502285."},
    {"Trappist-1b", ""},
    {"Trappist-1c", ""},
    {"Trappist-1d", ""},
    {"Trappist-1e", ""},
    {"Trappist-1f", ""},
    {"Trappist-1g", ""},
    {"Trappist-1h", ""}};
std::map<std::string, std::string> const body_name_map = {
    {"Trappist-1", "Sun"},
    {"Trappist-1b", "Beta"},
    {"Trappist-1c", "Charlie"},
    {"Trappist-1d", "Delta"},
    {"Trappist-1e", "Echo"},
    {"Trappist-1f", "Foxtrot"},
    {"Trappist-1g", "Golf"},
    {"Trappist-1h", "Hotel"}};

}  // namespace

void GenerateKopernicusForSlippist1(
    std::string const& gravity_model_stem,
    std::string const& initial_state_stem) {
  std::filesystem::path const directory =
      SOLUTION_DIR / "astronomy";
  SolarSystem<ICRFJ2000Equator> solar_system(
      (directory / gravity_model_stem).replace_extension(proto_txt),
      (directory / initial_state_stem).replace_extension(proto_txt),
      /*ignore_frame=*/true);

  std::ofstream kopernicus_cfg(
      (directory / (gravity_model_stem + "_slippist1")).replace_extension(cfg));
  CHECK(kopernicus_cfg.good());
  kopernicus_cfg << "@Kopernicus:FOR[Principia]:AFTER[aSLIPPIST-1] {\n";
  for (std::string const& name : solar_system.names()) {
    serialization::GravityModel::Body const& body =
        solar_system.gravity_model_message(name);
    bool const is_star =
        !solar_system.keplerian_initial_state_message(name).has_parent();
    kopernicus_cfg << "@Body[" << FindOrDie(body_name_map, name) << "] {\n";
    if (is_star) {
      kopernicus_cfg << "  @cbNameLater =" << name << "\n";
    }
    kopernicus_cfg << "  @Properties {\n";
    if (is_star) {
      kopernicus_cfg << "    !mass = delete\n";
    } else {
      kopernicus_cfg << "    !geeASL = delete\n";
    }
    kopernicus_cfg << "    %gravParameter = "
                   << GravitationalConstant * ParseQuantity<Mass>(body.mass())
                   << "\n";
    kopernicus_cfg << "    %radius = " << body.mean_radius() << "\n";
    kopernicus_cfg << "    %description = "
                   << FindOrDie(body_description_map, name) << "\n";
    kopernicus_cfg << "    %rotationPeriod = "
                   << 2 * π /
                          ParseQuantity<AngularFrequency>(
                              body.angular_frequency())
                   << "\n";
    if (!is_star) {
      kopernicus_cfg << "    %initialRotation = " << body.reference_angle()
                     << "\n";
      kopernicus_cfg << "    %tidallyLocked = false\n";
    }
    kopernicus_cfg << "  }\n";
  }
  kopernicus_cfg << "}\n";
}

}  // namespace tools
}  // namespace principia
