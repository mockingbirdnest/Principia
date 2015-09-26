#pragma once

#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace testing_utilities {

namespace {

void AdjustAccuracy(SolarSystemFactory::Accuracy const accuracy,
                    SolarSystem<ICRFJ2000Equator>* const solar_system) {
  switch (accuracy) {
    case SolarSystemFactory::Accuracy::kMajorBodiesOnly:
      break;
    case SolarSystemFactory::Accuracy::kMinorAndMajorBodies:
      break;
    case SolarSystemFactory::Accuracy::kAllBodiesAndOblateness:
      break;
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
}

not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>>
SolarSystemFactory::AtСпутник2Launch(Accuracy const accuracy) {
  auto solar_system =
      base::make_not_null_unique<SolarSystem<ICRFJ2000Equator>>();
  solar_system->Initialize(
      SOLUTION_DIR "astronomy\\gravity_model.proto.txt",
      SOLUTION_DIR "astronomy\\initial_state_jd_2436145_604166667.proto.txt");
}

int SolarSystemFactory::parent(int const index) {
}

std::string SolarSystemFactory::name(int const index) {
}

}  // namespace testing_utilities
}  // namespace principia
