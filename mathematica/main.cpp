
#include "glog/logging.h"
#include "mathematica/integrator_plots.hpp"

int main(int argc, char const* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::LogToStderr();
  principia::mathematica::GenerateSimpleHarmonicMotionWorkErrorGraphs();
  principia::mathematica::GenerateKeplerProblemWorkErrorGraphs();
  return 0;
}
