
#include <string>

#include "glog/logging.h"
#include "mathematica/integrator_plots.hpp"
#include "mathematica/ksp_dynamical_stability.hpp"

int main(int argc, char const* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::LogToStderr();
  CHECK_LE(argc, 2);
  if (argc == 1 || std::string(argv[1]) == "integrator_plots") {
    principia::mathematica::GenerateSimpleHarmonicMotionWorkErrorGraphs();
    // Circular.
    principia::mathematica::GenerateKeplerProblemWorkErrorGraphs(0.0);
    // Pluto.
    principia::mathematica::GenerateKeplerProblemWorkErrorGraphs(0.25);
    // 67P.
    principia::mathematica::GenerateKeplerProblemWorkErrorGraphs(0.64);
    // 1P.
    principia::mathematica::GenerateKeplerProblemWorkErrorGraphs(0.967);
  } else if (std::string(argv[1]) == "ksp_dynamical_stability") {
    principia::mathematica::SimulateStockSystem();
    principia::mathematica::SimulateFixedSystem();
  } else {
    LOG(FATAL) << "unexpected argument " << argv[1];
  }
  return 0;
}
