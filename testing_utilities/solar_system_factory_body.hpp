#pragma once

#include <set>
#include <string>

#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace testing_utilities {

namespace {

#define ADD_BODY(x, names)                                                     \
do {                                                                           \
  std::string const x_string = #x;                                             \
  names.insert(x_string.substr(1));                                            \
} while (true)

void AdjustAccuracy(SolarSystemFactory::Accuracy const accuracy,
                    SolarSystem<ICRFJ2000Equator>* const solar_system) {
  std::set<std::string> existing_names;
  std::set<std::string> oblate_names;
  switch (accuracy) {
    case SolarSystemFactory::Accuracy::kAllBodiesAndOblateness:
      ADD_BODY(kJupiter, oblate_names);
      ADD_BODY(kSaturn, oblate_names);
      ADD_BODY(kNeptune, oblate_names);
      ADD_BODY(kUranus, oblate_names);
    case SolarSystemFactory::Accuracy::kMinorAndMajorBodies:
      ADD_BODY(kTitania, existing_names);
      ADD_BODY(kOberon, existing_names);
      ADD_BODY(kRhea, existing_names);
      ADD_BODY(kIapetus, existing_names);
      ADD_BODY(kCharon, existing_names);
      ADD_BODY(kAriel, existing_names);
      ADD_BODY(kUmbriel, existing_names);
      ADD_BODY(kDione, existing_names);
      ADD_BODY(kTethys, existing_names);
    case SolarSystemFactory::Accuracy::kMajorBodiesOnly:
      ADD_BODY(kSun, existing_names);
      ADD_BODY(kJupiter, existing_names);
      ADD_BODY(kSaturn, existing_names);
      ADD_BODY(kNeptune, existing_names);
      ADD_BODY(kUranus, existing_names);
      ADD_BODY(kEarth, existing_names);
      ADD_BODY(kVenus, existing_names);
      ADD_BODY(kMars, existing_names);
      ADD_BODY(kMercury, existing_names);
      ADD_BODY(kGanymede, existing_names);
      ADD_BODY(kTitan, existing_names);
      ADD_BODY(kCallisto, existing_names);
      ADD_BODY(kIo, existing_names);
      ADD_BODY(kMoon, existing_names);
      ADD_BODY(kEuropa, existing_names);
      ADD_BODY(kTriton, existing_names);
      ADD_BODY(kEris, existing_names);
      ADD_BODY(kPluto, existing_names);
  }
  for (std::string const& solar_system_name : solar_system->names()) {
    if (existing_names.count(solar_system_name) == 0) {
      solar_system->RemoveMassiveBody(solar_system_name);
    } else if (oblate_names.count(solar_system_name) == 0) {
      solar_system->RemoveOblateness(solar_system_name);
    }
  }
}

}  // namespace

not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>>
SolarSystemFactory::AtСпутник1Launch(Accuracy const accuracy) {
  auto solar_system =
      base::make_not_null_unique<SolarSystem<ICRFJ2000Equator>>();
  solar_system->Initialize(
      SOLUTION_DIR "astronomy\\gravity_model.proto.txt",
      SOLUTION_DIR "astronomy\\initial_state_jd_2436116_311504629.proto.txt");
  AdjustAccuracy(accuracy, solar_system.get());
  return solar_system;
}

not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>>
SolarSystemFactory::AtСпутник2Launch(Accuracy const accuracy) {
  auto solar_system =
      base::make_not_null_unique<SolarSystem<ICRFJ2000Equator>>();
  solar_system->Initialize(
      SOLUTION_DIR "astronomy\\gravity_model.proto.txt",
      SOLUTION_DIR "astronomy\\initial_state_jd_2436145_604166667.proto.txt");
  AdjustAccuracy(accuracy, solar_system.get());
  return solar_system;
}

int SolarSystemFactory::parent(int const index) {
}

std::string SolarSystemFactory::name(int const index) {
}

}  // namespace testing_utilities
}  // namespace principia
