
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <utility>
#include <vector>
#include <map>
#include <string>
#include <ranges>

#include <intrin.h>
#include <emmintrin.h>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "base/cpuid.hpp"
#include "numerics/cbrt.hpp"
#include "testing_utilities/statistics.hpp"
#include "nanobenchmarks/disable_perf_boost_mode.hpp"

using namespace principia::nanobenchmarks::_disable_perf_boost_mode;

#define BENCHMARK_CALLING_CONVENTION

using BenchmarkedFunction = double(BENCHMARK_CALLING_CONVENTION*)(double);

static std::map<std::string, BenchmarkedFunction>&
    function_registry =
    *new std::map<std::string, BenchmarkedFunction>();

#define BENCHMARK_FUNCTION(f) \
  static bool registered_##f = function_registry.emplace(#f, &(f)).second

#define EXPAND(x) x

#define BENCHMARK_FUNCTION_WITH_NAME(name, ...) \
  BENCHMARK_FUNCTION_WITH_NAME_INTERNAL(__LINE__, name, __VA_ARGS__)
#define BENCHMARK_FUNCTION_WITH_NAME_INTERNAL(line, name, ...) \
  BENCHMARK_FUNCTION_WITH_NAME_INTERNAL2(line, name, __VA_ARGS__)
#define BENCHMARK_FUNCTION_WITH_NAME_INTERNAL2(line, name, ...) \
  static bool registered_##line =                              \
      function_registry.emplace(name, &(__VA_ARGS__)).second

#define BenchmarkedFunction(f)                    \
  double BENCHMARK_CALLING_CONVENTION f(double x); \
  BENCHMARK_FUNCTION(f);                           \
  double BENCHMARK_CALLING_CONVENTION f(double x)

#define BENCHMARK_EXTERN_C_FUNCTION(f) \
extern "C" double BENCHMARK_CALLING_CONVENTION f(double); \
BENCHMARK_FUNCTION(f)

BENCHMARK_EXTERN_C_FUNCTION(identity);
BENCHMARK_EXTERN_C_FUNCTION(sqrtps_xmm0_xmm0);
BENCHMARK_EXTERN_C_FUNCTION(sqrtsd_xmm0_xmm0);
BENCHMARK_EXTERN_C_FUNCTION(mulsd_xmm0_xmm0);
BENCHMARK_EXTERN_C_FUNCTION(mulsd_xmm0_xmm0_4x);

BenchmarkedFunction(twice) {
  return 2 * x;
}

BenchmarkedFunction(thrice) {
  return 3 * x;
}

BenchmarkedFunction(inc) {
  return x + 1;
}

BenchmarkedFunction(add_4_times) {
  return x * x * x * x * x;
}

BenchmarkedFunction(add_16_times) {

  return x + x + x + x +
         x + x + x + x +
         x + x + x + x +
         x + x + x + x + x;}

BenchmarkedFunction(square_root) {
  __m128d x_0 = _mm_set_sd(x);
  return _mm_cvtsd_f64(_mm_sqrt_sd(x_0, x_0));
}

BenchmarkedFunction(sqrt_sqrt) {
  __m128d x_0 = _mm_set_sd(x);
  x_0 = _mm_sqrt_sd(x_0, x_0);
  return _mm_cvtsd_f64(_mm_sqrt_sd(x_0, x_0));
}

BenchmarkedFunction(square_root_division) {
  __m128d x_0 = _mm_set_sd(x);
  return _mm_cvtsd_f64(_mm_div_sd(x_0, _mm_sqrt_sd(x_0, x_0)));
}

struct distribution {
  double min;
  std::vector<double> quantiles;
  static std::vector<double>& quantile_definitions;

  static std::ostream& __cdecl heading(std::ostream& out) {
    out << std::setw(8) << "min";
    for (auto const& n : quantile_definitions) {
      if (n > 100) {
        std::print(out, "{:>7}‰", 1000.0 / n);
      } else {
        std::print(out, "{:>7}%", 100.0 / n);
      }
    }
    return out;
  }
};

std::vector<double>& distribution::quantile_definitions =
    *new std::vector<double>();

std::ostream& operator<<(std::ostream& out, distribution const& x) {
  std::print(out, "{:8.2f}", x.min);
  for (double const quantile : x.quantiles) {
      std::print(out, "{:+8.2f}", quantile - x.min);
  }
  return out;
}

distribution operator*(double a, distribution x) {
  distribution result{a * x.min};
  for (double const quantile : x.quantiles) {
      result.quantiles.push_back(a * quantile);
  }
  return result;
}

distribution operator+(distribution x, double b) {
  distribution result{x.min + b};
  for (double const quantile : x.quantiles) {
      result.quantiles.push_back(quantile + b);
  }
  return result;
}

__declspec(noinline) distribution benchmark(BenchmarkedFunction f) {
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
  distribution result {
      durations[0]};
  for (int const q : distribution::quantile_definitions) {
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
  PerfBoostModeDisabler perf_boost_mode_disabled;
  distribution::quantile_definitions = {1000, 100, 10, 4, 2};
  std::cout << principia::base::_cpuid::CPUVendorIdentificationString() << " "
            << principia::base::_cpuid::ProcessorBrandString() << "\nFeatures:"
            << principia::base::_cpuid::CPUFeatures() << "\n";
  auto name_widths = std::views::keys(function_registry) |
                    std::views::transform(&std::string::size);
  std::size_t name_width = *std::ranges::max_element(name_widths);
  std::vector reference_functions{
      std::pair{&identity, 0},
      std::pair{&mulsd_xmm0_xmm0, 4},
      std::pair{&mulsd_xmm0_xmm0_4x, 4 * 4},
      std::pair{&sqrtps_xmm0_xmm0, 12},
  };
  std::map<BenchmarkedFunction, distribution>
      reference_measurements;
  std::cout << std::setw(name_width + 2) << "RAW TSC:" << distribution::heading
            << "\n";
  for (auto const& [function, _] : reference_functions) {
    auto const result = benchmark(function);
    reference_measurements.emplace(function, result);
    auto const& name =
        std::ranges::find(function_registry, function, [&](auto pair) {
          return pair.second;
        })->first;
    std::vprint_unicode(std::cout,
                        " {:>" + std::to_string(name_width + 1) + "}",
                        std::make_format_args(name));
    std::cout << result << "\n";
  }
  std::vector<double> tsc;
  std::vector<double> expected_cycles;
  for (auto const& [f, cycles] : reference_functions) {
    tsc.push_back(reference_measurements[f].min);
    expected_cycles.push_back(cycles);
  }
  double const a =
      principia::testing_utilities::_statistics::Slope(tsc, expected_cycles);
  double const b =
      principia::testing_utilities::_statistics::Mean(expected_cycles) -
      a * principia::testing_utilities::_statistics::Mean(tsc);
  std::cout << "Slope: " << std::setprecision(6) << a << " cycle/TSC\n";
  std::cout << "Correlation coefficient: "
            << principia::testing_utilities::_statistics::
                   PearsonProductMomentCorrelationCoefficient(tsc,
                                                              expected_cycles)
            << "\n";
  std::cout << std::setw(name_width + 2) << "Cycles:" << distribution::heading
            << "\n";
  auto bm_cycles = [&](BenchmarkedFunction f) {
    return a * benchmark(f) + b;
  };
  for (auto const& [name, f] : function_registry) {
    auto const result = benchmark(f);
    std::cout << (std::ranges::contains(std::views::keys(reference_functions),
                                        f)
                      ? "R"
                      : " ");
    std::vprint_unicode(std::cout,
                        "{:>" + std::to_string(name_width + 1) + "}",
                        std::make_format_args(name));
    std::cout << a * result + b << "\n";
  }
}