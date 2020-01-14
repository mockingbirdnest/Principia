
#include "tools/generate_kopernicus.hpp"

#include <filesystem>
#include <map>
#include <string>

#include "astronomy/frames.hpp"
#include "base/map_util.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/solar_system.hpp"
#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/parser.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using astronomy::Sky;
using base::FindOrDie;
using physics::KeplerOrbit;
using physics::SolarSystem;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::DebugString;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Mass;
using quantities::Mod;
using quantities::ParseQuantity;
using quantities::SIUnit;
using quantities::constants::GravitationalConstant;
using quantities::si::Degree;
using quantities::si::Metre;
using quantities::si::Radian;

namespace tools {

namespace {

constexpr char cfg[] = "cfg";
constexpr char proto_txt[] = "proto.txt";
constexpr char kerbin[] = "Kerbin";

std::map<std::string, Angle> const body_angle = {
    {"Trappist-1b", -90 * Degree},
    {"Trappist-1c", -90 * Degree},
    {"Trappist-1d", -90 * Degree},
    {"Trappist-1e", -90 * Degree},
    {"Trappist-1f", 180 * Degree},
    {"Trappist-1g", -90 * Degree},
    {"Trappist-1h", -90 * Degree}};
std::map<std::string, std::string> const body_description_map = {
    {"Trappist-1",
     "An ultra-cool red dwarf star of spectral type M8 V, 40 light-years from "
     "Earth; also known as 2MUCD 12171 and 2MASS J23062928-0502285."},
    {"Trappist-1b",
     "Requires envelope of volatiles in the form of a thick atmosphere.  Above "
     "the runaway greenhouse limit."},
    {"Trappist-1c", "Likely has a rocky interior."},
    {"Trappist-1d",
     "Requires envelope of volatiles in the form of a thick atmosphere, oceans "
     "or ice."},
    {"Trappist-1e", "Likely has a rocky interior."},
    {"Trappist-1f",
     "Requires envelope of volatiles in the form of oceans or ice."},
    {"Trappist-1g",
     "Requires envelope of volatiles in the form of oceans or ice."},
    {"Trappist-1h",
     "Requires envelope of volatiles in the form of oceans or ice."}};
std::map<std::string, std::string> const body_name_map = {
    {"Trappist-1", "Sun"},
    {"Trappist-1b", "Bravo"},
    {"Trappist-1c", "Charlie"},
    {"Trappist-1d", "Delta"},
    {"Trappist-1e", "Kerbin"},
    {"Trappist-1f", "Foxtrot"},
    {"Trappist-1g", "Golf"},
    {"Trappist-1h", "Hotel"}};

}  // namespace

void GenerateKopernicusForSlippist1(
    std::string const& gravity_model_stem,
    std::string const& initial_state_stem) {
  std::filesystem::path const directory =
      SOLUTION_DIR / "astronomy";
  SolarSystem<Sky> solar_system(
      (directory / gravity_model_stem).replace_extension(proto_txt),
      (directory / initial_state_stem).replace_extension(proto_txt),
      /*ignore_frame=*/true);

  std::ofstream kopernicus_cfg(
      (directory / (gravity_model_stem + "_slippist1")).replace_extension(cfg));
  CHECK(kopernicus_cfg.good());

  // Find the star.  This is needed to construct Kepler orbits below.
  std::optional<serialization::GravityModel::Body> star;
  for (std::string const& name : solar_system.names()) {
    serialization::GravityModel::Body const& body =
        solar_system.gravity_model_message(name);
    bool const is_star =
        !solar_system.keplerian_initial_state_message(name).has_parent();
    if (is_star) {
      star = body;
    }
  }

  kopernicus_cfg << "@Kopernicus:AFTER[aSLIPPIST-1] {\n";
  for (std::string const& name : solar_system.names()) {
    serialization::GravityModel::Body const& body =
        solar_system.gravity_model_message(name);
    serialization::InitialState::Keplerian::Body::Elements const& elements =
        solar_system.keplerian_initial_state_message(name).elements();
    bool const is_star =
        !solar_system.keplerian_initial_state_message(name).has_parent();
    bool const is_kerbin =
        FindOrDie(body_name_map, name) == kerbin;
    kopernicus_cfg << "  @Body[" << FindOrDie(body_name_map, name) << "] {\n";
    if (is_kerbin) {
      kopernicus_cfg << "    %cbNameLater = " << name << "\n";
    } else if (!is_star) {
      kopernicus_cfg << "    @name = " << name << "\n";
    }
    kopernicus_cfg << "    @Properties {\n";
    if (is_star) {
      kopernicus_cfg << "      !mass = delete\n";
    } else {
      kopernicus_cfg << "      !geeASL = delete\n";
    }
    kopernicus_cfg << "      @displayName = " << name << "\n";
    kopernicus_cfg << "      %gravParameter = "
                   << DebugString(ParseQuantity<GravitationalParameter>(
                                      body.gravitational_parameter()) /
                                  SIUnit<GravitationalParameter>())
                   << "\n";
    kopernicus_cfg << "      %radius = "
                   << DebugString(ParseQuantity<Length>(body.mean_radius()) /
                                  Metre)
                   << "\n";
    kopernicus_cfg << "      %description = "
                   << FindOrDie(body_description_map, name) << "\n";
    if (!is_star) {
      kopernicus_cfg << "      %tidallyLocked = false\n";
    }
    kopernicus_cfg << "    }\n";
    if (is_star) {
      kopernicus_cfg << "    @ScaledVersion {\n";
      kopernicus_cfg << "      @Light {\n";
      for (char const* const curve :
           {"ScaledIntensityCurve", "IntensityCurve", "IVAIntensityCurve"}) {
        kopernicus_cfg << "        @" << curve << " {\n";
        kopernicus_cfg << "          @key,*[0, ] *= 10\n";
        kopernicus_cfg << "        }\n";
      }
      kopernicus_cfg << "      }\n";
      kopernicus_cfg << "    }\n";
    } else {
      CHECK(star.has_value());
      auto const keplerian_elements =
          solar_system.MakeKeplerianElements(elements);
      KeplerOrbit<Sky> const kepler_orbit(*solar_system.MakeMassiveBody(*star),
                                          *solar_system.MakeMassiveBody(body),
                                          keplerian_elements,
                                          solar_system.epoch());
      kopernicus_cfg << "    @Orbit {\n";
      kopernicus_cfg << "      %semiMajorAxis = "
                     << DebugString(
                            *kepler_orbit.elements_at_epoch().semimajor_axis /
                            Metre)
                     << "\n";
      kopernicus_cfg << "      %eccentricity = "
                     << DebugString(elements.eccentricity()) << "\n";
      kopernicus_cfg << "      %longitudeOfAscendingNode = "
                     << DebugString(
                            ParseQuantity<Angle>(
                                elements.longitude_of_ascending_node()) /
                            Degree)
                     << "\n";
      kopernicus_cfg << "      %argumentOfPeriapsis = "
                     << DebugString(ParseQuantity<Angle>(
                                        elements.argument_of_periapsis()) /
                                    Degree)
                     << "\n";
      kopernicus_cfg << "      %meanAnomalyAtEpoch = "
                     << DebugString(ParseQuantity<Angle>(
                                        elements.mean_anomaly()) / Radian)
                     << "\n";
      kopernicus_cfg << "    }\n";
      kopernicus_cfg << "    @Atmosphere {\n";
      kopernicus_cfg << "      @altitude *= 2\n";
      for (char const* const curve :
           {"temperatureCurve", "temperatureSunMultCurve", "pressureCurve"}) {
        kopernicus_cfg << "      @" << curve << " {\n";
        kopernicus_cfg << "        @key,*[0, ] *= 2\n";
        kopernicus_cfg << "      }\n";
      }
      kopernicus_cfg << "    }\n";
    }
    kopernicus_cfg << "  }\n";
  }
  kopernicus_cfg << "}\n";

  kopernicus_cfg << "@principia_gravity_model:FOR[Principia] {\n";
  for (std::string const& name : solar_system.names()) {
    serialization::InitialState::Keplerian::Body::Elements const& elements =
        solar_system.keplerian_initial_state_message(name).elements();
    bool const is_star =
        !solar_system.keplerian_initial_state_message(name).has_parent();
    kopernicus_cfg << "  @body[" << name << "] {\n";
    if (!is_star) {
      kopernicus_cfg << "      @reference_angle = "
                     << Mod(FindOrDie(body_angle, name) +
                                ParseQuantity<Angle>(
                                    elements.argument_of_periapsis()) +
                                ParseQuantity<Angle>(elements.mean_anomaly()),
                            2 * π * Radian)
                     << "\n";
    }
    kopernicus_cfg << "  }\n";
  }
  kopernicus_cfg << "}\n";
  for (std::string const& name : solar_system.names()) {
    bool const is_star =
        !solar_system.keplerian_initial_state_message(name).has_parent();
    if (!is_star) {
      kopernicus_cfg << "@Scatterer_atmosphere:HAS[@Atmo["
                     << body_name_map.at(name) << "]]:AFTER[aSLIPPIST-1] {\n";
      kopernicus_cfg << "  @Atmo[" << body_name_map.at(name) << "] {\n";
      kopernicus_cfg << "    @name = " << name << "\n";
      kopernicus_cfg << "    @configPoints {\n";
      kopernicus_cfg << "      @Item,* {\n";
      kopernicus_cfg << "        @altitude *= 2\n";
      kopernicus_cfg << "      }\n";
      kopernicus_cfg << "    }\n";
      kopernicus_cfg << "  }\n";
      kopernicus_cfg << "}\n";
      kopernicus_cfg << "@Scatterer_ocean:HAS[@Ocean[" << body_name_map.at(name)
                     << "]]:AFTER[aSLIPPIST-1] {\n";
      kopernicus_cfg << "  @Ocean[" << body_name_map.at(name) << "] {\n";
      kopernicus_cfg << "    @name = " << name << "\n";
      kopernicus_cfg << "  }\n";
      kopernicus_cfg << "}\n";
    }
  }
  kopernicus_cfg << "@Scatterer_planetsList:AFTER[aSLIPPIST-1] {\n";
  kopernicus_cfg << "  @scattererCelestialBodies {\n";
  for (std::string const& name : solar_system.names()) {
    bool const is_star =
        !solar_system.keplerian_initial_state_message(name).has_parent();
    bool const is_kerbin = body_name_map.at(name) == kerbin;
    std::string const slippist_name =
        is_kerbin ? "Echo" : body_name_map.at(name);
    if (!is_star) {
      kopernicus_cfg << "    @Item[" << slippist_name << "] {\n";
      kopernicus_cfg << "      @celestialBodyName = " << name << "\n";
      kopernicus_cfg << "      @transformName = " << name << "\n";
      kopernicus_cfg << "    }\n";
    }
  }
  kopernicus_cfg << "  }\n";
  kopernicus_cfg << "}\n";
  kopernicus_cfg << "@EVE_CLOUDS:AFTER[aSLIPPIST-1] {\n";
  for (std::string const& name : solar_system.names()) {
    bool const is_star =
        !solar_system.keplerian_initial_state_message(name).has_parent();
    if (!is_star) {
      kopernicus_cfg << "  @OBJECT:HAS[#body[" << body_name_map.at(name)
                     << "]] {\n";
      kopernicus_cfg << "    @body = " << name << "\n";
      kopernicus_cfg << "    @altitude *= 2\n";
      kopernicus_cfg << "  }\n";
    }
  }
  kopernicus_cfg << "}\n";
  kopernicus_cfg << "@EVE_SHADOWS:AFTER[aSLIPPIST-1] {\n";
  for (std::string const& name : solar_system.names()) {
    bool const is_star =
        !solar_system.keplerian_initial_state_message(name).has_parent();
    if (!is_star) {
      kopernicus_cfg << "  @OBJECT:HAS[#body[" << body_name_map.at(name)
                     << "]] {\n";
      kopernicus_cfg << "    @body = " << name << "\n";
      kopernicus_cfg << "    !caster,* = delete\n";
      auto const elements = SolarSystem<Sky>::MakeKeplerianElements(
          solar_system.keplerian_initial_state_message(name).elements());
      for (std::string const& caster_name : solar_system.names()) {
        bool const caster_is_star =
            !solar_system.keplerian_initial_state_message(caster_name)
                 .has_parent();
        if (!caster_is_star) {
          auto const caster_elements = SolarSystem<Sky>::MakeKeplerianElements(
              solar_system.keplerian_initial_state_message(caster_name)
                  .elements());
          if (caster_elements.period < elements.period) {
            kopernicus_cfg << "    caster = " << caster_name << "\n";
          }
        }
      }
      kopernicus_cfg << "  }\n";
    }
  }
  kopernicus_cfg << "}\n";
}

}  // namespace tools
}  // namespace principia
