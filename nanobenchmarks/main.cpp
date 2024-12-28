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

#include "base/cpuid.hpp"
#include "numerics/cbrt.hpp"
#include "testing_utilities/statistics.hpp"

#define BENCHMARK_CALLING_CONVENTION

static std::map<std::string, double(BENCHMARK_CALLING_CONVENTION*)(double)>&
    function_registry =
        *new std::map<std::string,
                      double(BENCHMARK_CALLING_CONVENTION*)(double)>();

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

#define BENCHMARKED_FUNCTION(f)                    \
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

BENCHMARKED_FUNCTION(twice) {
  return 2 * x;
}

BENCHMARKED_FUNCTION(thrice) {
  return 3 * x;
}

BENCHMARKED_FUNCTION(inc) {
  return x + 1;
}

BENCHMARKED_FUNCTION(add_4_times) {
  return x * x * x * x * x;
}

BENCHMARKED_FUNCTION(add_16_times) {

  return x + x + x + x +
         x + x + x + x +
         x + x + x + x +
         x + x + x + x + x;}

BENCHMARKED_FUNCTION(square_root) {
  __m128d x_0 = _mm_set_sd(x);
  return _mm_cvtsd_f64(_mm_sqrt_sd(x_0, x_0));
}

BENCHMARKED_FUNCTION(sqrt_sqrt) {
  __m128d x_0 = _mm_set_sd(x);
  x_0 = _mm_sqrt_sd(x_0, x_0);
  return _mm_cvtsd_f64(_mm_sqrt_sd(x_0, x_0));
}

BENCHMARKED_FUNCTION(square_root_division) {
  __m128d x_0 = _mm_set_sd(x);
  return _mm_cvtsd_f64(_mm_div_sd(x_0, _mm_sqrt_sd(x_0, x_0)));
}

struct distribution {
  double min;
  double percentile;
  double decile;
  double quartile;
  double median;

  static std::ostream& __cdecl heading(std::ostream& out) {
    return out << std::setw(8) << "min" << std::setw(8) << "1%" << std::setw(8)
               << "10%" << std::setw(8) << "25%" << std::setw(8) << "50%";
  }
};

std::ostream& operator<<(std::ostream& out, distribution const& x) {
  return out << std::fixed << std::setprecision(2) << std::setw(8) << x.min
             << std::showpos << std::setw(8) << x.percentile - x.min
             << std::setw(8) << x.decile - x.min << std::setw(8)
             << x.quartile - x.min << std::setw(8) << x.median - x.min
             << std::noshowpos << std::defaultfloat;
}

distribution operator*(double a, distribution x) {
  return {a * x.min, a * x.percentile, a * x.decile, a * x.quartile, a * x.median};
}

distribution operator+(distribution x, double b) {
  return {x.min + b, x.percentile + b, x.decile + b, x.quartile + b, x.median + b};
}

__declspec(noinline) distribution benchmark(double (BENCHMARK_CALLING_CONVENTION* f)(double)) {
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
  return {
      durations[0],
      durations[k / 100],
      durations[k / 10],
      durations[k / 4],
      durations[k / 2],
  };
}

BENCHMARK_FUNCTION_WITH_NAME(
    "Cbrt no FMA faithful",
    principia::numerics::_cbrt::internal::method_3²ᴄZ5¹::Cbrt<
        principia::numerics::_cbrt::internal::Rounding::Faithful>);
BENCHMARK_FUNCTION_WITH_NAME(
    "Cbrt no FMA correct",
    principia::numerics::_cbrt::internal::method_3²ᴄZ5¹::Cbrt<
        principia::numerics::_cbrt::internal::Rounding::Correct>);
BENCHMARK_FUNCTION_WITH_NAME(
    "Cbrt FMA faithful",
    principia::numerics::_cbrt::internal::method_5²Z4¹FMA::Cbrt<
        principia::numerics::_cbrt::internal::Rounding::Faithful>);
BENCHMARK_FUNCTION_WITH_NAME(
    "Cbrt FMA correct",
    principia::numerics::_cbrt::internal::method_5²Z4¹FMA::Cbrt<
        principia::numerics::_cbrt::internal::Rounding::Correct>);
BENCHMARK_FUNCTION_WITH_NAME("Cbrt",
    principia::numerics::_cbrt::Cbrt);

int __cdecl main(int argc, const char** argv) {
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
  std::map<double(BENCHMARK_CALLING_CONVENTION*)(double), distribution>
      reference_measurements;
  std::cout << std::setw(name_width + 1) << "RAW TSC:" << distribution::heading
            << "\n";
  for (auto const& [function, _] : reference_functions) {
    auto const result = benchmark(function);
    reference_measurements.emplace(function, result);
    auto const& name =
        std::ranges::find(function_registry, function, [&](auto pair) {
          return pair.second;
        })->first;
    std::cout << std::setw(name_width + 1) << name << result << "\n";
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
  std::cout << std::setw(name_width + 1) << "Cycles:" << distribution::heading
            << "\n";
  auto bm_cycles = [&](double(BENCHMARK_CALLING_CONVENTION * f)(double)) {
    return a * benchmark(f) + b;
  };
  for (auto const& [name, f] : function_registry) {
    std::cout << (std::ranges::contains(reference_functions,
                                        f,
                                        [&](auto pair) { return pair.first; })
                      ? "R"
                      : " ")
              << std::setw(name_width + 1) << name << a * benchmark(f) + b
              << "\n";
  }
}