#include <map>
#include <string>
#include <vector>

#include "base/map_util.hpp"
#include "base/not_null.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "integrators/integrators.hpp"
#include "mathematica/integrator_plots.hpp"
#include "mathematica/mathematica.hpp"
#include "mathematica/local_error_analysis.hpp"
#include "mathematica/retrobop_dynamical_stability.hpp"
#include "quantities/parser.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using base::Contains;
using base::make_not_null_unique;
using quantities::ParseQuantity;
using integrators::ParseFixedStepSizeIntegrator;

namespace mathematica {

class ErrorAnalysisTest : public ::testing::Test {
 protected:
};

TEST_F(ErrorAnalysisTest, DISABLED_IntegratorPlots) {
  GenerateSimpleHarmonicMotionWorkErrorGraphs();
  // Circular.
  GenerateKeplerProblemWorkErrorGraphs(0.0);
  // Pluto.
  GenerateKeplerProblemWorkErrorGraphs(0.25);
  // 67P.
  GenerateKeplerProblemWorkErrorGraphs(0.64);
  // 1P.
  GenerateKeplerProblemWorkErrorGraphs(0.967);
}

TEST_F(ErrorAnalysisTest, RetrobopDynamicalStability_PlotPredictableYears) {
  PlotPredictableYears();
}

TEST_F(ErrorAnalysisTest, DISABLED_RetrobopDynamicalStability_PlotCentury) {
  PlotCentury();
}

TEST_F(ErrorAnalysisTest,
       DISABLED_RetrobopDynamicalStability_AnalyseGlobalError) {
  AnalyseGlobalError();
}

TEST_F(ErrorAnalysisTest,
       DISABLED_RetrobopDynamicalStability_StatisticallyAnalyseStability) {
  google::LogToStderr();
  StatisticallyAnalyseStability();
}

TEST_F(ErrorAnalysisTest, LocalErrorAnalysis) {
  google::LogToStderr();
  std::vector<std::string> argv = ::testing::internal::GetArgvs();
  std::map<std::string, std::optional<std::string>> flags;
  for (int i = 2; i < argv.size(); ++i) {
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
    flags.emplace(flag_name, std::nullopt);
  }
  if (flags.empty() || Contains(flags, "help") || Contains(flags, "?")) {
    // Example:
    // .\Release\x64\mathematica_tests.exe \
    //   --gtest_filter=ErrorAnalysisTest.LocalErrorAnalysis \
    //   --gravity_model=.\astronomy\kerbol_gravity_model.proto.txt \
    //   --initial_state=.\astronomy\kerbol_initial_state_0_0.proto.txt \
    //   --time_step=1h --integrator=QUINLAN_1999_ORDER_8A
    std::printf(
        "Usage:\n"
        "mathematica_tests --gtest_filter=%s.%s \n"
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
        testing::UnitTest::GetInstance()->current_test_info()->test_case_name(),
        testing::UnitTest::GetInstance()->current_test_info()->name());
    return;
  }
  auto const gravity_model_path =
      std::filesystem::path(*flags["gravity_model"]);
  auto const initial_state_path =
      std::filesystem::path(*flags["initial_state"]);
  auto solar_system = make_not_null_unique<SolarSystem<ICRS>>(
      gravity_model_path, initial_state_path, /*ignore_frame=*/true);
  auto const& integrator =
      ParseFixedStepSizeIntegrator<Ephemeris<ICRS>::NewtonianMotionEquation>(
          *flags["integrator"]);
  auto const time_step = ParseQuantity<Time>(*flags["time_step"]);
  auto const out =
      std::filesystem::path(flags["output_directory"].value_or(".")) /
      (std::string("local_error_analysis[") + solar_system->names()[0] + "," +
       solar_system->epoch_literal() + "," + *flags["integrator"] + "," +
       DebugString(time_step) + "].wl");
  LocalErrorAnalyser analyser(std::move(solar_system), integrator, time_step);
  analyser.WriteLocalErrors(
      out,
      ParseFixedStepSizeIntegrator<Ephemeris<ICRS>::NewtonianMotionEquation>(
          flags["fine_integrator"].value_or("BLANES_MOAN_2002_SRKN_14A")),
      ParseQuantity<Time>(flags["fine_step"].value_or("1 min")),
      ParseQuantity<Time>(flags["granularity"].value_or("1 d")),
      ParseQuantity<Time>(flags["duration"].value_or("500 d")));
}

}  // namespace mathematica
}  // namespace principia
