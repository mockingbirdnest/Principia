
#pragma once

#include <experimental/filesystem>
#include <set>
#include <string>
#include <vector>

#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace testing_utilities {

namespace {

void AdjustAccuracy(SolarSystemFactory::Accuracy const accuracy,
                    SolarSystem<ICRFJ2000Equator>* const solar_system) {
  std::set<std::string> existing;
  std::set<std::string> oblate;
  switch (accuracy) {
    case SolarSystemFactory::Accuracy::all_bodies_and_oblateness:
      oblate.insert(SolarSystemFactory::name(SolarSystemFactory::jupiter));
      oblate.insert(SolarSystemFactory::name(SolarSystemFactory::saturn));
      oblate.insert(SolarSystemFactory::name(SolarSystemFactory::neptune));
      oblate.insert(SolarSystemFactory::name(SolarSystemFactory::uranus));
    case SolarSystemFactory::Accuracy::minor_and_major_bodies:
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::titania));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::oberon));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::rhea));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::iapetus));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::charon));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::ariel));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::umbriel));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::dione));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::tethys));
    case SolarSystemFactory::Accuracy::major_bodies_only:
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::sun));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::jupiter));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::saturn));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::neptune));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::uranus));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::earth));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::venus));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::mars));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::mercury));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::ganymede));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::titan));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::callisto));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::io));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::moon));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::europa));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::triton));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::eris));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::pluto));
  }
  std::vector<std::string> bodies_to_remove;
  std::vector<std::string> bodies_to_spherify;
  for (std::string const& solar_system_name : solar_system->names()) {
    if (existing.count(solar_system_name) == 0) {
      bodies_to_remove.push_back(solar_system_name);
    } else if (oblate.count(solar_system_name) == 0) {
      bodies_to_spherify.push_back(solar_system_name);
    }
  }
  for (std::string const& body_to_remove : bodies_to_remove) {
    solar_system->RemoveMassiveBody(body_to_remove);
  }
  for (std::string const& body_to_spherify : bodies_to_spherify) {
    solar_system->RemoveOblateness(body_to_spherify);
  }
}

}  // namespace

inline not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>>
SolarSystemFactory::AtСпутник1Launch(Accuracy const accuracy) {
  auto solar_system =
      base::make_not_null_unique<SolarSystem<ICRFJ2000Equator>>();
  solar_system->Initialize(
      SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "initial_state_jd_2436116_311504629.proto.txt");
  AdjustAccuracy(accuracy, solar_system.get());
  return solar_system;
}

inline not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>>
SolarSystemFactory::AtСпутник2Launch(Accuracy const accuracy) {
  auto solar_system =
      base::make_not_null_unique<SolarSystem<ICRFJ2000Equator>>();
  solar_system->Initialize(
      SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "initial_state_jd_2436145_604166667.proto.txt");
  AdjustAccuracy(accuracy, solar_system.get());
  return solar_system;
}

inline int SolarSystemFactory::parent(int const index) {
  switch (index) {
    case sun:
      LOG(FATAL) << FUNCTION_SIGNATURE << "The Sun has no parent";
      base::noreturn();
    case jupiter:
    case saturn:
    case neptune:
    case uranus:
    case earth:
    case venus:
    case mars:
    case mercury:
    case eris:
    case pluto:
      return sun;
    case ganymede:
    case callisto:
    case io:
    case europa:
      return jupiter;
    case titan:
    case rhea:
    case iapetus:
    case dione:
    case tethys:
      return saturn;
    case moon:
      return earth;
    case triton:
      return neptune;
    case titania:
    case oberon:
    case ariel:
    case umbriel:
      return uranus;
    case charon:
      return pluto;
    default:
      LOG(FATAL) << FUNCTION_SIGNATURE << "Undefined index";
      base::noreturn();
  }
}

inline std::string SolarSystemFactory::name(int const index) {
#define BODY_NAME(name) case k##name: return #name
  switch (index) {
    BODY_NAME(Sun);
    BODY_NAME(Jupiter);
    BODY_NAME(Saturn);
    BODY_NAME(Neptune);
    BODY_NAME(Uranus);
    BODY_NAME(Earth);
    BODY_NAME(Venus);
    BODY_NAME(Mars);
    BODY_NAME(Mercury);
    BODY_NAME(Ganymede);
    BODY_NAME(Titan);
    BODY_NAME(Callisto);
    BODY_NAME(Io);
    BODY_NAME(Moon);
    BODY_NAME(Europa);
    BODY_NAME(Triton);
    BODY_NAME(Eris);
    BODY_NAME(Pluto);
    BODY_NAME(Titania);
    BODY_NAME(Oberon);
    BODY_NAME(Rhea);
    BODY_NAME(Iapetus);
    BODY_NAME(Charon);
    BODY_NAME(Ariel);
    BODY_NAME(Umbriel);
    BODY_NAME(Dione);
    BODY_NAME(Tethys);
    default:
      LOG(FATAL) << FUNCTION_SIGNATURE << "Undefined index";
      base::noreturn();
  }
#undef BODY_NAME
}

}  // namespace testing_utilities
}  // namespace principia
