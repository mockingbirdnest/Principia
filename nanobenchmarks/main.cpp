
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <print>
#include <ranges>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <intrin.h>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "base/cpuid.hpp"
#include "nanobenchmarks/function_registry.hpp"
#include "nanobenchmarks/performance_settings_controller.hpp"
#include "numerics/cbrt.hpp"
#include "testing_utilities/statistics.hpp"

using namespace principia::nanobenchmarks::_function_registry;
using namespace principia::nanobenchmarks::_performance_settings_controller;
using namespace principia::testing_utilities::_statistics;

BENCHMARK_EXTERN_C_FUNCTION(identity);
BENCHMARK_EXTERN_C_FUNCTION(sqrtps_xmm0_xmm0);
BENCHMARK_EXTERN_C_FUNCTION(sqrtsd_xmm0_xmm0);
BENCHMARK_EXTERN_C_FUNCTION(mulsd_xmm0_xmm0);
BENCHMARK_EXTERN_C_FUNCTION(mulsd_xmm0_xmm0_4x);

struct LatencyDistributionTable {
  double min;
  std::vector<double> quantiles;
  static std::vector<double>& quantile_definitions;

  static std::string heading() {
    std::stringstream out;
    std::print(out, "{:>8}", "min");
    for (auto const& n : quantile_definitions) {
      if (n > 1000) {
        std::print(out, "{:>7}‱", 10'000.0 / n);
      } else if (n > 100) {
        std::print(out, "{:>7}‰", 1000.0 / n);
      } else {
        std::print(out, "{:>7}%", 100.0 / n);
      }
    }
    return out.str();
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

std::vector<double>& LatencyDistributionTable::quantile_definitions =
    *new std::vector<double>();

LatencyDistributionTable operator*(double a, LatencyDistributionTable x) {
  LatencyDistributionTable result{a * x.min};
  for (double const quantile : x.quantiles) {
      result.quantiles.push_back(a * quantile);
  }
  return result;
}

LatencyDistributionTable operator+(LatencyDistributionTable x, double b) {
  LatencyDistributionTable result{x.min + b};
  for (double const quantile : x.quantiles) {
      result.quantiles.push_back(quantile + b);
  }
  return result;
}

__declspec(noinline) LatencyDistributionTable benchmark(BenchmarkedFunction f) {
  constexpr int k = 1'000'000;
  static double* durations = new double[k];
  int registers[4]{};
  int leaf = 0;
  for (int j = 0; j < k; ++j) {
    constexpr int n = 100;
    __cpuid(registers, leaf);
    auto const tsc = __rdtsc();
    double x = 5 + tsc % 2 + registers[0] % 2;
    for (int i = 0; i < n; ++i) {
      x = f(x);
      x += 5 - x;
    }
    __cpuid(registers, x);
    double const δtsc = __rdtsc() - tsc;
    durations[j] = δtsc / n;
  }
  std::sort(durations, durations + k);
  LatencyDistributionTable result {
      durations[0]};
  for (int const q : LatencyDistributionTable::quantile_definitions) {
    result.quantiles.push_back(durations[k / q]);
  }
  return result;
}

BENCHMARK_FUNCTION_WITH_NAME(
    "Cbrt 3²ᴄZ5¹ Faithful",
    principia::numerics::_cbrt::internal::method_3²ᴄZ5¹::Cbrt<
        principia::numerics::_cbrt::internal::Rounding::Faithful>);
BENCHMARK_FUNCTION_WITH_NAME(
    "Cbrt 3²ᴄZ5¹ Correct",
    principia::numerics::_cbrt::internal::method_3²ᴄZ5¹::Cbrt<
        principia::numerics::_cbrt::internal::Rounding::Correct>);
BENCHMARK_FUNCTION_WITH_NAME(
    "Cbrt 5²Z4¹FMA Faithful",
    principia::numerics::_cbrt::internal::method_5²Z4¹FMA::Cbrt<
        principia::numerics::_cbrt::internal::Rounding::Faithful>);
BENCHMARK_FUNCTION_WITH_NAME(
    "Cbrt 5²Z4¹FMA Correct",
    principia::numerics::_cbrt::internal::method_5²Z4¹FMA::Cbrt<
        principia::numerics::_cbrt::internal::Rounding::Correct>);
BENCHMARK_FUNCTION_WITH_NAME("Cbrt",
    principia::numerics::_cbrt::Cbrt);

int __cdecl main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);
  auto controller = PerformanceSettingsController::Make();
  LatencyDistributionTable::quantile_definitions = {1000, 100, 20, 10, 4, 2};
  std::println("{} {}",
               principia::base::_cpuid::CPUVendorIdentificationString(),
               principia::base::_cpuid::ProcessorBrandString());
  std::println("Features: {}", principia::base::_cpuid::CPUFeatures());
  auto name_widths =
      std::views::keys(FunctionRegistry::singleton.functions_by_name()) |
      std::views::transform(&std::string::size);
  std::size_t name_width = *std::ranges::max_element(name_widths);
  std::vector reference_functions{
      std::pair{&identity, 0},
      std::pair{&mulsd_xmm0_xmm0, 4},
      std::pair{&mulsd_xmm0_xmm0_4x, 4 * 4},
      std::pair{&sqrtps_xmm0_xmm0, 12},
  };
  std::map<BenchmarkedFunction, LatencyDistributionTable>
      reference_measurements;
  std::vprint_unicode(
      "{:<" + std::to_string(name_width + 2) + "}{}\n",
      std::make_format_args("RAW TSC:", LatencyDistributionTable::heading()));
  for (auto const& [function, _] : reference_functions) {
    auto const result = benchmark(function);
    reference_measurements.emplace(function, result);
    std::vprint_unicode(
        "{:>" + std::to_string(name_width + 2) + "}{}\n",
        std::make_format_args(
            FunctionRegistry::singleton.names_by_function().at(function),
            result.Row()));
  }
  std::vector<double> tsc;
  std::vector<double> expected_cycles;
  for (auto const& [f, cycles] : reference_functions) {
    tsc.push_back(reference_measurements[f].min);
    expected_cycles.push_back(cycles);
  }
  double const a = Slope(tsc, expected_cycles);
  double const b = Mean(expected_cycles) - a * Mean(tsc);
  std::println("Slope: {:0.6f} cycle/TSC", a);
  std::println(
      "Correlation coefficient: {:0.6f}",
      PearsonProductMomentCorrelationCoefficient(tsc, expected_cycles));
  std::vprint_unicode(
      "{:<" + std::to_string(name_width + 2) + "}{}\n",
      std::make_format_args("Cycles:", LatencyDistributionTable::heading()));
  auto bm_cycles = [&](BenchmarkedFunction f) {
    return a * benchmark(f) + b;
  };
  for (auto const& [name, f] : FunctionRegistry::singleton.functions_by_name()) {
    auto const result = benchmark(f);
    std::vprint_unicode(
        "{}{:>" + std::to_string(name_width + 1) + "}{}\n",
        std::make_format_args(
            std::ranges::contains(std::views::keys(reference_functions), f)
                ? "R"
                : " ",
            name,
            (a * result + b).Row()));
  }
}