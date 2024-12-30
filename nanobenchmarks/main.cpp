
#include <algorithm>
#include <iostream>
#include <iomanip>
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
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_CLANG.
#include "mathematica/logger.hpp"
#include "mathematica/mathematica.hpp"
#include "nanobenchmarks/function_registry.hpp"
#include "nanobenchmarks/microarchitectures.hpp"
#include "nanobenchmarks/performance_settings_controller.hpp"
#include "testing_utilities/statistics.hpp"


#if PRINCIPIA_COMPILER_MSVC
#include <intrin.h>
#endif

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

// Adding support for flag types only works using ADL (or by being in
// marshalling.h), so we do this, which is UB.
namespace std {

bool AbslParseFlag(absl::string_view const text,
                   std::vector<double>* const flag,
                   std::string* const error) {
  flag->clear();
  for (absl::string_view const element : absl::StrSplit(text, ',')) {
    if (!absl::ParseFlag(element, &flag->emplace_back(), error)) {
      return false;
    }
  }
  return true;
}

std::string AbslUnparseFlag(std::vector<double> const& flag) {
  return absl::StrJoin(flag, ",");
}

}  // namespace std

ABSL_FLAG(std::vector<double>,
          quantiles,
          (std::vector<double>{0.001, 0.01, 0.05, 0.1, 0.25, 0.5}),
          "Quantiles to report");

namespace principia {
namespace nanobenchmarks {
namespace _main {
namespace {

using namespace principia::mathematica::_logger;
using namespace principia::mathematica::_mathematica;
using namespace principia::nanobenchmarks::_function_registry;
using namespace principia::nanobenchmarks::_microarchitectures;
using namespace principia::nanobenchmarks::_performance_settings_controller;
using namespace principia::testing_utilities::_statistics;

struct LatencyDistributionTable {
  double min;
  std::vector<double> quantiles;

  static std::string const& Heading() {
    static std::string const& result = [] {
      std::stringstream& out = *new std::stringstream();
      std::print(out, "{:>8}", "min");
      for (double const q : absl::GetFlag(FLAGS_quantiles)) {
        if (q < 1e-3) {
          std::print(out, "{:>7}‱", 10'000 * q);
        } else if (q < 1e-2) {
          std::print(out, "{:>7}‰", 1000 * q);
        } else {
          std::print(out, "{:>7}%", 100 * q);
        }
      }
      return out.str();
    }();
    return result;
  }

  std::string Row() const {
    std::stringstream out;
    std::print(out, "{:8.2f}", min);
    for (double const quantile : quantiles) {
      std::print(out, "{:+8.2f}", quantile - min);
    }
    return out.str();
  }
};

LatencyDistributionTable operator*(double const a,
                                   LatencyDistributionTable const& x) {
  LatencyDistributionTable result{a * x.min};
  for (double const quantile : x.quantiles) {
      result.quantiles.push_back(a * quantile);
  }
  return result;
}

LatencyDistributionTable operator+(LatencyDistributionTable const& x,
                                   double const b) {
  LatencyDistributionTable result{x.min + b};
  for (double const quantile : x.quantiles) {
      result.quantiles.push_back(quantile + b);
  }
  return result;
}

// We disable inlining on this function so that the overhead is independent of
// the callsite, and so that we actually call the benchmarked function via a
// function pointer, instead of inlining it.
__declspec(noinline) LatencyDistributionTable
    Benchmark(BenchmarkedFunction const f, Logger* logger) {
  std::size_t const sample_count = absl::GetFlag(FLAGS_samples);
  std::size_t const loop_iterations = absl::GetFlag(FLAGS_loop_iterations);
  static std::vector<double>& samples = *new std::vector<double>(
      sample_count, std::numeric_limits<double>::quiet_NaN());
  int registers[4]{};
  int leaf = 0;
  for (int j = 0; j < sample_count; ++j) {
    double const input = absl::GetFlag(FLAGS_input);
    double x = input;
    // The CPUID barriers prevent out-of-order execution; see [Pao10].
    #if PRINCIPIA_COMPILER_MSVC
    __cpuid(registers, leaf);
    #else
    asm volatile("cpuid");
    #endif
    auto const tsc_start = __rdtsc();
    for (int i = 0; i < loop_iterations; ++i) {
      x = f(x);
      x += input - x;
    }
    unsigned int tsc_aux;
    // The use of RDTSCP rather than RDTSC here follows [Pao10].  See the Intel®
    // 64 and IA-32 Architectures Software Developer’s Manual:
    // The RDTSCP  instruction is not a serializing instruction, but it does
    // wait until all previous instructions have executed and all previous loads
    // are globally visible.  But it does not wait for previous stores to be
    // globally visible, and subsequent instructions may begin execution before
    // the read operation is performed.
    auto const tsc_stop = __rdtscp(&tsc_aux);
    #if PRINCIPIA_COMPILER_MSVC
    __cpuid(registers, leaf);
    #else
    asm volatile("cpuid");
    #endif
    double const δtsc = tsc_stop - tsc_start;
    samples[j] = δtsc / loop_iterations;
  }
  if (logger != nullptr) {
    logger->Append(
        "samples",
        std::tuple{FunctionRegistry::names_by_function().at(f),
                   samples});
  }
  std::ranges::sort(samples);
  LatencyDistributionTable result{samples[0]};
  for (double const q : absl::GetFlag(FLAGS_quantiles)) {
    result.quantiles.push_back(samples[(sample_count - 1) * q]);
  }
  return result;
}

std::size_t FormattedWidth(std::string const& s) {
  // Two columns per code unit is wide enough, since field width is at most 2
  // per extended grapheme cluster.
  std::size_t const wide = 2 * s.size();
    // There is no vformatted_size, so we actually format.
    std::size_t const formatted_size =
      std::vformat("{:" + std::to_string(wide) + "}",
                     std::make_format_args(s))
            .size();
  // The actual width is the field width we allocated, minus the padding spaces
  // added by formatting.
  return wide - (formatted_size - s.size());
}

void Main() {
  std::regex const name_matcher(absl::GetFlag(FLAGS_benchmark_filter));
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
  auto name_widths =
      FunctionRegistry::functions_by_name() |
      std::views::filter([&](auto const& pair) {
        auto const& [name, f] = pair;
        return std::regex_match(name, name_matcher) ||
               ReferenceCycleCounts().contains(f);
      }) |
      std::views::keys | std::views::transform(&FormattedWidth);
  std::size_t const name_width = *std::ranges::max_element(name_widths);
  std::map<BenchmarkedFunction, LatencyDistributionTable>
      reference_measurements;
  std::vprint_unicode(
      "{:<" + std::to_string(name_width + 2) + "}{:8}{}\n",
                      std::make_format_args(
                          "RAW TSC:", "", LatencyDistributionTable::Heading()));
  for (auto const& [function, _] : ReferenceCycleCounts()) {
    auto const result = Benchmark(function, logger.get());
    reference_measurements.emplace(function, result);
    std::vprint_unicode(
        "{:>" + std::to_string(name_width + 2) + "}{:8}{}\n",
        std::make_format_args(
                            FunctionRegistry::names_by_function().at(function),
                            "",
            static_cast<std::string const&>(result.Row())));
  }
  std::vector<double> tsc;
  std::vector<double> expected_cycles;
  for (auto const& [f, cycles] : ReferenceCycleCounts()) {
    tsc.push_back(reference_measurements[f].min);
    expected_cycles.push_back(cycles);
  }
  double const a = Slope(tsc, expected_cycles);
  double const b = Mean(expected_cycles) - a * Mean(tsc);
  auto benchmark_cycles = [&](BenchmarkedFunction const f) {
    return a * Benchmark(f, logger.get()) + b;
  };
  std::println("Slope: {:0.6f} cycle/TSC", a);
  std::println(
      "Correlation coefficient: {:0.6f}",
      PearsonProductMomentCorrelationCoefficient(tsc, expected_cycles));
  std::vprint_unicode(
      "{:<" + std::to_string(name_width + 2) + "}{:>8}{}\n",
      std::make_format_args(
          "Cycles:", "expected", LatencyDistributionTable::Heading()));
  for (auto const& [name, f] :
       FunctionRegistry::functions_by_name()) {
    if (!std::regex_match(name, name_matcher) &&
        !ReferenceCycleCounts().contains(f)) {
      continue;
    }
    std::vprint_unicode(
        "{} {:>" + std::to_string(name_width) + "}{:>8}{}\n",
        std::make_format_args(
            ReferenceCycleCounts().contains(f) ? "R" : " ",
            name,
            static_cast<std::string const&>(
                ReferenceCycleCounts().contains(f)
                    ? std::to_string(ReferenceCycleCounts().at(f))
                    : ""),
            static_cast<std::string const&>(benchmark_cycles(f).Row())));
  }
}

}  // namespace
}  // namespace _main
}  // namespace nanobenchmarks
}  // namespace principia

int __cdecl main(int const argc, char** const argv) {
  absl::ParseCommandLine(argc, argv);
  principia::nanobenchmarks::_main::Main();
}
