
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
    CHECK_EQ(argc, 3);
    std::string const target = argv[2];
    if (target == "plot_predictable_years") {
      principia::mathematica::PlotPredictableYears();
    } else if (target == "plot_century") {
      principia::mathematica::PlotCentury();
    } else if (target == "analyse_global_error") {
      principia::mathematica::AnalyseGlobalError();
    } else if (target == "analyse_local_error") {
      principia::mathematica::AnalyseLocalError();
    } else if (target == "statistically_analyse_stability") {
      principia::mathematica::StatisticallyAnalyseStability();
    } else {
      LOG(FATAL) << "unexpected target " << target;
    }
  } else {
    LOG(FATAL) << "unexpected argument " << argv[1];
  }
  return 0;
}
