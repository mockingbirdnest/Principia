
#include <experimental/filesystem>
#include <string>

#include "astronomy/frames.hpp"
#include "glog/logging.h"
#include "integrators/integrators.hpp"
#include "physics/ephemeris.hpp"
#include "mathematica/integrator_plots.hpp"
#include "mathematica/retrobop_dynamical_stability.hpp"
#include "physics/solar_system.hpp"
#include "quantities/parser.hpp"

using ::principia::astronomy::ICRFJ2000Equator;
using ::principia::quantities::ParseQuantity;
using ::principia::integrators::ParseFixedStepSizeIntegrator;
using ::principia::physics::Ephemeris;
using ::principia::physics::SolarSystem;
using ::principia::quantities::Time;

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
    CHECK_GE(argc, 3);
    std::string const target = argv[2];
    if (target == "plot_predictable_years") {
      principia::mathematica::PlotPredictableYears();
    } else if (target == "plot_century") {
      principia::mathematica::PlotCentury();
    } else if (target == "analyse_global_error") {
      principia::mathematica::AnalyseGlobalError();
    } else if (target == "statistically_analyse_stability") {
      principia::mathematica::StatisticallyAnalyseStability();
    } else {
      LOG(FATAL) << "unexpected target " << target;
    }
  } else if (std::string(argv[1]) == "local_error_analysis") {
    CHECK_EQ(argc, 6) << "Usage: " << argv[0]
                      << " <gravity model path> <initial state path> "
                         "<integrator> <time step>";
    auto const gravity_model_path =
        std::experimental::filesystem::path(argv[2]);
    auto const initial_state_path =
        std::experimental::filesystem::path(argv[3]);
    auto solar_system = std::make_unique<SolarSystem<ICRFJ2000Equator>>(
        gravity_model_path, initial_state_path, /*ignore_frame=*/true);
    auto const& integrator = ParseFixedStepSizeIntegrator<
        Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation>(argv[4]);
    auto const time_step = ParseQuantity<Time>(argv[5]);
  } else {
    LOG(FATAL) << "unexpected argument " << argv[1];
  }
  return 0;
}
