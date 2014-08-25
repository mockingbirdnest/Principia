
#include "clr_benchmarks_adapter/n_body_system.hpp"

#include "benchmarks/n_body_system.hpp"
#include "testing_utilities/solar_system.hpp"

using principia::testing_utilities::SolarSystem;

namespace principia {
namespace clr_benchmarks_adapter {

void NBodySystemCLRBenchmark::SimulateSolarSystem() {
  std::unique_ptr<SolarSystem> solar_system = SolarSystem::AtСпутникLaunch();
  benchmarks::SimulateSolarSystem(solar_system.get());
}

}  // namespace clr_benchmarks_adapter
}  // namespace principia
