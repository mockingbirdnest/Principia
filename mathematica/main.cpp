
#include <string>

#include "glog/logging.h"
#include "mathematica/integrator_plots.hpp"
#include "mathematica/retrobop_dynamical_stability.hpp"

int main(int argc, char const* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::LogToStderr();
  if (argc == 1 || std::string(argv[1]) == "integrator_plots") {
    CHECK_LE(argc, 2);
    principia::mathematica::GenerateSimpleHarmonicMotionWorkErrorGraphs();
    // Circular.
    principia::mathematica::GenerateKeplerProblemWorkErrorGraphs(0.0);
    // Pluto.
    principia::mathematica::GenerateKeplerProblemWorkErrorGraphs(0.25);
    // 67P.
    principia::mathematica::GenerateKeplerProblemWorkErrorGraphs(0.64);
    // 1P.
    principia::mathematica::GenerateKeplerProblemWorkErrorGraphs(0.967);
  } else if (std::string(argv[1]) == "retrobop_dynamical_stability") {
    bool produce_file = false;
    if (argc == 3) {
      if (std::string(argv[2]) == "--produce_file") {
        produce_file = true;
      } else {
        LOG(FATAL) << "unexpected argument " << argv[3];
      }
    }
    principia::mathematica::SimulateFixedSystem(produce_file);
  } else {
    LOG(FATAL) << "unexpected argument " << argv[1];
  }
  return 0;
}
