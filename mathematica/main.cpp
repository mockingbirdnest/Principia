
#include <cstdio>
#include <experimental/filesystem>
#include <experimental/optional>
#include <map>
#include <string>

#include "astronomy/frames.hpp"
#include "base/map_util.hpp"
#include "glog/logging.h"
#include "integrators/integrators.hpp"
#include "physics/ephemeris.hpp"
#include "mathematica/integrator_plots.hpp"
#include "mathematica/local_error_analysis.hpp"
#include "mathematica/retrobop_dynamical_stability.hpp"
#include "physics/solar_system.hpp"
#include "quantities/parser.hpp"

using ::principia::astronomy::ICRFJ2000Equator;
using ::principia::base::Contains;
using ::principia::base::make_not_null_unique;
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
    CHECK_EQ(argc, 3);
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
    std::map<std::string, std::experimental::optional<std::string>> flags;
    for (int i = 2; i < argc; ++i) {
      std::string const flag(argv[i]);
      std::size_t const name_begin = flag.find_first_not_of("-/");
      std::size_t const name_end = flag.find('=');
      if (name_begin == std::string::npos) {
        LOG(FATAL) << "Invalid flag syntax '" << flag << "', expected\n"
                   << "(--|/)<flag_name>[=value]";
      }
      std::string const flag_name =
          flag.substr(name_begin, name_end - name_begin);
      if (name_end != std::string::npos) {
        flags.emplace(flag_name, flag.substr(name_end + 1, std::string::npos));
      }
      flags.emplace(flag_name, std::experimental::nullopt);
    }
    if (flags.empty() || Contains(flags, "help") || Contains(flags, "?")) {
      std::printf(
          "Usage:\n"
          "%s %s\n"
          "    --gravity_model=<path>\n"
          "    --initial_state=<path>\n"
          "    --integrator=<fixed_step_size_integrator >\n"
          "    --time_step=<quantity(time)>\n"
          "    [--output_directory=<path>] default: .\n"
          "    [--fine_integrator=<fixed_step_size_integrator >] "
          "        default: BLANES_MOAN_2002_SRKN_14A\n"
          "    [--fine_step=<quantity(time)>] default: 1 min\n"
          "    [--granularity=<quantity(time)>] default: 1 d\n"
          "    [--duration=<quantity(time)>] default: 500 d\n",
          argv[0],
          argv[1]);
      return 0;
    }
    auto const gravity_model_path =
        std::experimental::filesystem::path(*flags["gravity_model"]);
    auto const initial_state_path =
        std::experimental::filesystem::path(*flags["initial_state"]);
    auto solar_system = make_not_null_unique<SolarSystem<ICRFJ2000Equator>>(
        gravity_model_path, initial_state_path, /*ignore_frame=*/true);
    auto const& integrator = ParseFixedStepSizeIntegrator<
        Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation>(
        *flags["integrator"]);
    auto const time_step = ParseQuantity<Time>(*flags["time_step"]);
    auto const out =
        std::experimental::filesystem::path(
            flags["output_directory"].value_or(".")) /
        (std::string("local_error_analysis[") + solar_system->names()[0] + "," +
         solar_system->epoch_literal() + "," + *flags["integrator"] + "," +
         DebugString(time_step) + "].wl");
    principia::mathematica::LocalErrorAnalyser analyser(
        std::move(solar_system), integrator, time_step);
    analyser.WriteDailyErrors(
        out,
        ParseFixedStepSizeIntegrator<
            Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation>(
            flags["integrator"].value_or("BLANES_MOAN_2002_SRKN_14A")),
        ParseQuantity<Time>(flags["fine_step"].value_or("1 min")),
        ParseQuantity<Time>(flags["granularity"].value_or("1 d")),
        ParseQuantity<Time>(flags["duration"].value_or("500 d")));
  } else {
    LOG(FATAL) << "unexpected argument " << argv[1];
  }
  return 0;
}
