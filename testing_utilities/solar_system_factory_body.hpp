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
    case SolarSystemFactory::Accuracy::kAllBodiesAndOblateness:
      oblate.insert(SolarSystemFactory::name(SolarSystemFactory::kJupiter));
      oblate.insert(SolarSystemFactory::name(SolarSystemFactory::kSaturn));
      oblate.insert(SolarSystemFactory::name(SolarSystemFactory::kNeptune));
      oblate.insert(SolarSystemFactory::name(SolarSystemFactory::kUranus));
    case SolarSystemFactory::Accuracy::kMinorAndMajorBodies:
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kTitania));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kOberon));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kRhea));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kIapetus));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kCharon));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kAriel));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kUmbriel));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kDione));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kTethys));
    case SolarSystemFactory::Accuracy::kMajorBodiesOnly:
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kSun));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kJupiter));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kSaturn));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kNeptune));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kUranus));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kEarth));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kVenus));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kMars));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kMercury));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kGanymede));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kTitan));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kCallisto));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kIo));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kMoon));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kEuropa));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kTriton));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kEris));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::kPluto));
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
    case kSun:
      LOG(FATAL) << FUNCTION_SIGNATURE << "The Sun has no parent";
      base::noreturn();
    case kJupiter:
    case kSaturn:
    case kNeptune:
    case kUranus:
    case kEarth:
    case kVenus:
    case kMars:
    case kMercury:
    case kEris:
    case kPluto:
      return kSun;
    case kGanymede:
    case kCallisto:
    case kIo:
    case kEuropa:
      return kJupiter;
    case kTitan:
    case kRhea:
    case kIapetus:
    case kDione:
    case kTethys:
      return kSaturn;
    case kMoon:
      return kEarth;
    case kTriton:
      return kNeptune;
    case kTitania:
    case kOberon:
    case kAriel:
    case kUmbriel:
      return kUranus;
    case kCharon:
      return kPluto;
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
