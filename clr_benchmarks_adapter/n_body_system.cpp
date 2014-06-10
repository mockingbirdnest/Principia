
#include "clr_benchmarks_adapter/n_body_system.hpp"

#include "benchmarks/n_body_system.hpp"

namespace principia {
namespace clr_benchmarks_adapter {

void NBodySystemCLRBenchmark::SimulateSolarSystem() {
  benchmarks::SimulateSolarSystem();
}

}  // namespace clr_benchmarks_adapter
}  // namespace principia
