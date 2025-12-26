
#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <print>
#include <ranges>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_CLANG.
#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "mathematica/logger.hpp"
#include "mathematica/mathematica.hpp"
#include "nanobenchmarks/flag_parsing.hpp"  // 🧙 For std::vector-valued flags.
#include "nanobenchmarks/latency_distribution_table.hpp"
#include "nanobenchmarks/microarchitectures.hpp"
#include "nanobenchmarks/nanobenchmark.hpp"
#include "nanobenchmarks/performance_settings_controller.hpp"
#include "testing_utilities/statistics.hpp"

ABSL_FLAG(std::size_t,
          loop_iterations,
          100,
          "Number of iterations of the measured loop");
ABSL_FLAG(std::size_t,
          samples,
          1'000'000,
          "Number of measurements to perform for each benchmarked function");
ABSL_FLAG(std::string,
          benchmark_filter,
          ".*",
          "Regular expression matching the names of functions to benchmark");
ABSL_FLAG(std::string,
          log_to_mathematica,
          "",
          "File to which to log the measurements");
ABSL_FLAG(double, input, 5, "Input for the benchmarked functions");
ABSL_FLAG(std::vector<double>,
          quantiles,
          (std::vector<double>{0.001, 0.01, 0.05, 0.1, 0.25, 0.5}),
          "Quantiles to report");

namespace principia {
namespace nanobenchmarks {
namespace _main {
namespace {

using namespace principia::geometry::_frame;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::mathematica::_logger;
using namespace principia::mathematica::_mathematica;
using namespace principia::nanobenchmarks::_nanobenchmark;
using namespace principia::nanobenchmarks::_latency_distribution_table;
using namespace principia::nanobenchmarks::_microarchitectures;
using namespace principia::nanobenchmarks::_performance_settings_controller;
using namespace principia::testing_utilities::_statistics;

using World = Frame<serialization::Frame::TestTag,
                    Inertial,
                    Handedness::Right,
                    serialization::Frame::TEST>;

std::size_t FormattedWidth(std::string const& s) {
  // Two columns per code unit is wide enough, since field width is at most 2
  // per extended grapheme cluster.
  std::size_t const wide = 2 * s.size();
  // There is no vformatted_size, so we actually format.
  std::size_t const formatted_size =
      std::vformat("{:" + std::to_string(wide) + "}", std::make_format_args(s))
          .size();
  // The actual width is the field width we allocated, minus the padding spaces
  // added by formatting.
  return wide - (formatted_size - s.size());
}

template<typename Value, typename Argument>
void CalibrateAndRun(std::regex const& filter, Logger* const logger) {
  auto const nanobenchmarks =
      NanobenchmarkRegistry<Value, Argument>::NanobenchmarksMatching(filter);
  auto const& reference_cycle_counts = ReferenceCycleCounts();

  // Would like to use std::views::concat, but not this year.
  auto const nanobenchmark_widths =
      nanobenchmarks |
      std::views::transform(&Nanobenchmark<Value, Argument>::name) |
      std::views::transform(&FormattedWidth);
  auto const reference_cycle_counts_widths =
      reference_cycle_counts | std::views::keys |
      std::views::transform(&Nanobenchmark<>::name) |
      std::views::transform(&FormattedWidth);

  std::size_t name_width = std::ranges::max(reference_cycle_counts_widths);
  if (!std::ranges::empty(nanobenchmark_widths)) {
    name_width = std::max(name_width, std::ranges::max(nanobenchmark_widths));
  }

  std::map<Nanobenchmark<> const*, LatencyDistributionTable>
      reference_measurements;
  std::vprint_unicode(stdout,
                      "{:<" + std::to_string(name_width + 2) + "}{:8}{}\n",
                      std::make_format_args(
                          "RAW TSC:", "", LatencyDistributionTable::Heading()));
  for (auto const& [nanobenchmark, _] : reference_cycle_counts) {
    auto const result = nanobenchmark->Run(logger);
    reference_measurements.emplace(nanobenchmark, result);
    std::vprint_unicode(
        stdout,
        "{:>" + std::to_string(name_width + 2) + "}{:8}{}\n",
        std::make_format_args(nanobenchmark->name(),
                              "",
                              static_cast<std::string const&>(result.Row())));
  }
  std::vector<double> tsc;
  std::vector<double> expected_cycles;
  for (auto const& [nanobenchmark, cycles] : reference_cycle_counts) {
    tsc.push_back(reference_measurements[nanobenchmark].min());
    expected_cycles.push_back(cycles);
  }
  double const a = Slope(tsc, expected_cycles);
  double const b = Mean(expected_cycles) - a * Mean(tsc);
  auto benchmark_cycles =
      [a, b, logger](Nanobenchmark<Value, Argument> const* nanobenchmark) {
        return a * nanobenchmark->Run(logger) + b;
      };
  std::println("Slope: {:0.6f} cycle/TSC    Overhead: {:0.6f} TSC", a, -b / a);
  std::println(
      "Correlation coefficient: {:0.6f}",
      PearsonProductMomentCorrelationCoefficient(tsc, expected_cycles));
  std::vprint_unicode(
      stdout,
      "{:<" + std::to_string(name_width + 2) + "}{:>8}{}\n",
      std::make_format_args(
          "Cycles:", "expected", LatencyDistributionTable::Heading()));

  for (auto const& [nanobenchmark, cycles] : reference_cycle_counts) {
    std::vprint_unicode(
        stdout,
        "R {:>" + std::to_string(name_width) + "}{:>8}{}\n",
        std::make_format_args(
            nanobenchmark->name(),
            static_cast<std::string const&>(std::to_string(cycles)),
            static_cast<std::string const&>(
                benchmark_cycles(nanobenchmark).Row())));
  }

  for (auto const* nanobenchmark : nanobenchmarks) {
    std::vprint_unicode(
        stdout,
        "  {:>" + std::to_string(name_width) + "}        {}\n",
        std::make_format_args(nanobenchmark->name(),
                              static_cast<std::string const&>(
                                  benchmark_cycles(nanobenchmark).Row())));
  }
}

void Main() {
  std::regex const filter(absl::GetFlag(FLAGS_benchmark_filter));
  auto controller = PerformanceSettingsController::New();
  std::unique_ptr<Logger> logger;
  std::string const& filename = absl::GetFlag(FLAGS_log_to_mathematica);
  if (!filename.empty()) {
    logger = std::make_unique<Logger>(filename, /*make_unique=*/false);
  }
  std::println("{} {}",
               principia::base::_cpuid::CPUVendorIdentificationString(),
               principia::base::_cpuid::ProcessorBrandString());
  std::println("Features: {}", principia::base::_cpuid::CPUFeatures());

  CalibrateAndRun<double, double>(filter, logger.get());
}

}  // namespace
}  // namespace _main
}  // namespace nanobenchmarks
}  // namespace principia

int __cdecl main(int const argc, char** const argv) {
  absl::ParseCommandLine(argc, argv);
  principia::nanobenchmarks::_main::Main();
}
