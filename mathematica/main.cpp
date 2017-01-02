
#include "glog/logging.h"
#include "mathematica/integrator_plots.hpp"

int main(int argc, char const* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::LogToStderr();/*
  principia::mathematica::GenerateSimpleHarmonicMotionWorkErrorGraphs();
  // Circular.
  principia::mathematica::GenerateKeplerProblemWorkErrorGraphs(0.0);*/
  // Pluto.
  principia::mathematica::GenerateKeplerProblemWorkErrorGraphs(0.25);/*
  // 67P.
  principia::mathematica::GenerateKeplerProblemWorkErrorGraphs(0.64);
  // 1P.
  principia::mathematica::GenerateKeplerProblemWorkErrorGraphs(0.967);*/
  return 0;
}
