
#pragma once

#include <cctype>
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
    case SolarSystemFactory::Accuracy::AllBodiesAndOblateness:
      oblate.insert(SolarSystemFactory::name(SolarSystemFactory::Jupiter));
      oblate.insert(SolarSystemFactory::name(SolarSystemFactory::Saturn));
      oblate.insert(SolarSystemFactory::name(SolarSystemFactory::Neptune));
      oblate.insert(SolarSystemFactory::name(SolarSystemFactory::Uranus));
    case SolarSystemFactory::Accuracy::MinorAndMajorBodies:
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Titania));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Oberon));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Rhea));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Iapetus));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Charon));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Ariel));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Umbriel));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Dione));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Tethys));
    case SolarSystemFactory::Accuracy::MajorBodiesOnly:
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Sun));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Jupiter));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Saturn));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Neptune));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Uranus));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Earth));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Venus));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Mars));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Mercury));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Ganymede));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Titan));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Callisto));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Io));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Moon));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Europa));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Triton));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Eris));
      existing.insert(SolarSystemFactory::name(SolarSystemFactory::Pluto));
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
    case Sun:
      LOG(FATAL) << FUNCTION_SIGNATURE << "The Sun has no parent";
      base::noreturn();
    case Jupiter:
    case Saturn:
    case Neptune:
    case Uranus:
    case Earth:
    case Venus:
    case Mars:
    case Mercury:
    case Eris:
    case Pluto:
      return Sun;
    case Ganymede:
    case Callisto:
    case Io:
    case Europa:
      return Jupiter;
    case Titan:
    case Rhea:
    case Iapetus:
    case Dione:
    case Tethys:
      return Saturn;
    case Moon:
      return Earth;
    case Triton:
      return Neptune;
    case Titania:
    case Oberon:
    case Ariel:
    case Umbriel:
      return Uranus;
    case Charon:
      return Pluto;
    default:
      LOG(FATAL) << FUNCTION_SIGNATURE << "Undefined index";
      base::noreturn();
  }
}

inline std::string SolarSystemFactory::name(int const index) {
#define BODY_NAME(name) case name: return #name
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
