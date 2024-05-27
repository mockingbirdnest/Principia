// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=BM_Experiment  // NOLINT(whitespace/line_length)

#include <array>
#include <cmath>
#include <cstdint>
#include <random>

#include "benchmark/benchmark.h"
#include "numerics/double_precision.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"  // 🧙 For π.

namespace principia {
namespace functions {

// TODO(phl): The polynomials in this file should use class
// |PolynomialInMonomialBasis|.

using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_elementary_functions;

using Value = double;
using Argument = double;

constexpr Argument x_min = π / 6;  // The sinus is greater than 1/2.
constexpr Argument x_max = π / 4;  // Upper bound after argument reduction.
static constexpr std::int64_t number_of_iterations = 1000;

// A helper class to benchmark implementations with various table spacings.
template<Argument table_spacing>
class TableSpacingImplementation {
 public:
  void Initialize();

  Value Sin(Argument x);
  Value Cos(Argument x);

 private:
  // Despite the name these are not accurate values, but for the purposes of
  // benchmarking it doesn't matter.
  struct AccurateValues {
    Argument x;
    Value sin_x;
    Value cos_x;
  };

  Value SinPolynomial(Argument x);
  Value CosPolynomial(Argument x);

  std::array<AccurateValues,
             static_cast<std::int64_t>(x_max / table_spacing) + 1>
      accurate_values_;
};

template<Argument max_table_spacing, std::int64_t number_of_tables>
class MultiTableImplementation {
 public:
  void Initialize();

  Value Sin(Argument x);
  Value Cos(Argument x);

 private:
  // Despite the name these are not accurate values, but for the purposes of
  // benchmarking it doesn't matter.
  struct AccurateValues {
    Argument x;
    Value sin_x;
    Value cos_x;
  };

  //TODO(phl)templatize
  Value SinPolynomial(Argument x);
  Value CosPolynomial(Argument x);

  static constexpr std::int64_t table_size =
      static_cast<std::int64_t>(0.5 / max_table_spacing);

  std::array<std::int64_t, number_of_tables> one_over_table_spacings_;
  std::array<std::array<AccurateValues, table_size>, number_of_tables>
      accurate_values_;
};

template<Argument table_spacing>
void TableSpacingImplementation<table_spacing>::Initialize() {
  int i = 0;
  for (Argument x = table_spacing / 2;
       x <= x_max + table_spacing / 2;
       x += table_spacing, ++i) {
    accurate_values_[i] = {.x = x,
                           .sin_x = std::sin(x),
                           .cos_x = std::cos(x)};
  }
}

template<Argument table_spacing>
FORCE_INLINE(inline)
Value TableSpacingImplementation<table_spacing>::Sin(Argument const x) {
  auto const i = static_cast<std::int64_t>(x * (1.0 / table_spacing));
  auto const& accurate_values = accurate_values_[i];
  auto const& x₀ = accurate_values.x;
  auto const& sin_x₀ = accurate_values.sin_x;
  auto const& cos_x₀ = accurate_values.cos_x;
  auto const h = x - x₀;
  auto const h² = h * h;
  auto const h³ = h² * h;
  auto const sin_x₀_plus_h_cos_x₀ = TwoProductAdd(cos_x₀, h, sin_x₀);
  return sin_x₀_plus_h_cos_x₀.value +
         ((sin_x₀ * h² * CosPolynomial(h²) + cos_x₀ * h³ * SinPolynomial(h²)) +
          sin_x₀_plus_h_cos_x₀.error);
}

template<Argument table_spacing>
FORCE_INLINE(inline)
Value TableSpacingImplementation<table_spacing>::Cos(Argument const x) {
  auto const i = static_cast<std::int64_t>(x * (1.0 / table_spacing));
  auto const& accurate_values = accurate_values_[i];
  auto const& x₀ = accurate_values.x;
  auto const& sin_x₀ = accurate_values.sin_x;
  auto const& cos_x₀ = accurate_values.cos_x;
  auto const h = x - x₀;
  auto const h² = h * h;
  auto const h³ = h² * h;
  auto const cos_x₀_minus_h_sin_x₀ = TwoProductNegatedAdd(sin_x₀, h, cos_x₀);
  return cos_x₀_minus_h_sin_x₀.value +
         ((cos_x₀ * h² * CosPolynomial(h²) - sin_x₀ * h³ * SinPolynomial(h²)) +
          cos_x₀_minus_h_sin_x₀.error);
}

template<Argument table_spacing>
Value TableSpacingImplementation<table_spacing>::SinPolynomial(Argument const x) {
  if constexpr (table_spacing == 2.0 / 256.0) {
    // 71 bits.
    return -0.166666666666666666666421797625 +
           0.00833333057503280528178543245797 * x;
  } else if constexpr (table_spacing == 2.0 / 1024.0) {
    // 85 bits.
    return -0.166666666666666666666666651721 +
           0.00833333316093951937646271666739 * x;
  }
}

template<Argument table_spacing>
Value TableSpacingImplementation<table_spacing>::CosPolynomial(Argument const x) {
  if constexpr (table_spacing == 2.0 / 256.0) {
    // 77 bits.
    return -0.499999999999999999999999974543 +
           x * (0.0416666666666633318024480868405 -
                0.00138888829905860875255146938745 * x);
  } else if constexpr (table_spacing == 2.0 / 1024.0) {
    // 72 bits.
    return -0.499999999999999999999872434553 +
           0.0416666654823785864634569932662 * x;
  }
}

template<Argument max_table_spacing, std::int64_t number_of_tables>
void MultiTableImplementation<max_table_spacing, number_of_tables>::
Initialize() {
  Argument current_x_max = 2 * x_min;
  Argument current_x_min = x_min;
  Argument current_table_spacing = max_table_spacing;
  for (std::int64_t i = number_of_tables - 1; i >= 0; --i) {
    one_over_table_spacings_[i] = 1.0 / current_table_spacing;
    std::int64_t j = number_of_tables - 1;
    for (Argument x = current_x_max - current_table_spacing / 2;
         x > current_x_min;
         x -= current_table_spacing, --j) {
      accurate_values_[i][j] = {.x = x,
                                .sin_x = std::sin(x),
                                .cos_x = std::cos(x)};
    }
    current_x_max = current_x_min;
    current_x_min /= 2;
    current_table_spacing /= 2;
  }
}

template<Argument max_table_spacing, std::int64_t number_of_tables>
FORCE_INLINE(inline)
Value MultiTableImplementation<max_table_spacing, number_of_tables>::
Sin(Argument const x) {
  int x_exponent;
  auto const x_mantissa = std::frexp(x, &x_exponent);
  auto const i = number_of_tables + x_exponent;
  auto const j = static_cast<std::int64_t>((x_mantissa - 0.5) *
                                           one_over_table_spacings_[i]);
  auto const& accurate_values = accurate_values_[i][j];
  auto const& x₀ = accurate_values.x;
  auto const& sin_x₀ = accurate_values.sin_x;
  auto const& cos_x₀ = accurate_values.cos_x;
  auto const h = x - x₀;
  auto const h² = h * h;
  auto const h³ = h² * h;
  auto const sin_x₀_plus_h_cos_x₀ = TwoProductAdd(cos_x₀, h, sin_x₀);
  return sin_x₀_plus_h_cos_x₀.value +
         ((sin_x₀ * h² * CosPolynomial(h²) + cos_x₀ * h³ * SinPolynomial(h²)) +
          sin_x₀_plus_h_cos_x₀.error);
}

template<Argument max_table_spacing, std::int64_t number_of_tables>
FORCE_INLINE(inline)
Value MultiTableImplementation<max_table_spacing, number_of_tables>::
Cos(Argument const x) {
  int x_exponent;
  auto const x_mantissa = std::frexp(x, &x_exponent);
  auto const i = number_of_tables + x_exponent;
  auto const j = static_cast<std::int64_t>((x_mantissa - 0.5) *
                                           one_over_table_spacings_[i]);
  auto const& accurate_values = accurate_values_[i][j];
  auto const& x₀ = accurate_values.x;
  auto const& sin_x₀ = accurate_values.sin_x;
  auto const& cos_x₀ = accurate_values.cos_x;
  auto const h = x - x₀;
  auto const h² = h * h;
  auto const h³ = h² * h;
  auto const cos_x₀_minus_h_sin_x₀ = TwoProductNegatedAdd(sin_x₀, h, cos_x₀);
  return cos_x₀_minus_h_sin_x₀.value +
         ((cos_x₀ * h² * CosPolynomial(h²) - sin_x₀ * h³ * SinPolynomial(h²)) +
          cos_x₀_minus_h_sin_x₀.error);
}

template<Argument max_table_spacing, std::int64_t number_of_tables>
Value MultiTableImplementation<max_table_spacing, number_of_tables>::
SinPolynomial(Argument const x) {
  // 85 bits.
  return -0.166666666666666666666666651721 +
         0.00833333316093951937646271666739 * x;
}

template<Argument max_table_spacing, std::int64_t number_of_tables>
Value MultiTableImplementation<max_table_spacing, number_of_tables>::
CosPolynomial(Argument const x) {
  // 72 bits.
  return -0.499999999999999999999872434553 +
          0.0416666654823785864634569932662 * x;
}

template<Argument table_spacing>
void BM_ExperimentSinTableSpacing(benchmark::State& state) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(x_min, x_max);

  TableSpacingImplementation<table_spacing> implementation;
  implementation.Initialize();

  Argument a[number_of_iterations];
  for (std::int64_t i = 0; i < number_of_iterations; ++i) {
    a[i] = uniformly_at(random);
  }

  Value v[number_of_iterations];
  while (state.KeepRunningBatch(number_of_iterations)) {
    for (std::int64_t i = 0; i < number_of_iterations; ++i) {
      using namespace principia::quantities;
      v[i] = implementation.Sin(a[i]);
#if _DEBUG
      // The implementation is not accurate, but let's check that it's not
      // broken.
      auto const absolute_error = Abs(v[i] - std::sin(a[i]));
      CHECK_LT(absolute_error, 5.6e-17);
#endif
    }
    benchmark::DoNotOptimize(v);
  }
}

template<Argument table_spacing>
void BM_ExperimentCosTableSpacing(benchmark::State& state) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(x_min, x_max);

  TableSpacingImplementation<table_spacing> implementation;
  implementation.Initialize();

  Argument a[number_of_iterations];
  for (std::int64_t i = 0; i < number_of_iterations; ++i) {
    a[i] = uniformly_at(random);
  }

  Value v[number_of_iterations];
  while (state.KeepRunningBatch(number_of_iterations)) {
    for (std::int64_t i = 0; i < number_of_iterations; ++i) {
      using namespace principia::quantities;
      v[i] = implementation.Cos(a[i]);
#if _DEBUG
      // The implementation is not accurate, but let's check that it's not
      // broken.
      auto const absolute_error = Abs(v[i] - std::cos(a[i]));
      CHECK_LT(absolute_error, 1.2e-16);
#endif
    }
    benchmark::DoNotOptimize(v);
  }
}

template<Argument max_table_spacing, std::int64_t number_of_tables>
void BM_ExperimentSinMultiTable(benchmark::State& state) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(
      std::scalbln(x_min, 1 - number_of_tables), x_max);

  MultiTableImplementation<max_table_spacing, number_of_tables> implementation;
  implementation.Initialize();

  Argument a[number_of_iterations];
  for (std::int64_t i = 0; i < number_of_iterations; ++i) {
    a[i] = uniformly_at(random);
  }

  Value v[number_of_iterations];
  while (state.KeepRunningBatch(number_of_iterations)) {
    for (std::int64_t i = 0; i < number_of_iterations; ++i) {
      using namespace principia::quantities;
      v[i] = implementation.Sin(a[i]);
#if _DEBUG
      // The implementation is not accurate, but let's check that it's not
      // broken.
      auto const absolute_error = Abs(v[i] - std::sin(a[i]));
      CHECK_LT(absolute_error, 5.6e-17);
#endif
    }
    benchmark::DoNotOptimize(v);
  }
}

template<Argument max_table_spacing, std::int64_t number_of_tables>
void BM_ExperimentCosMultiTable(benchmark::State& state) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(
      std::scalbln(x_min, 1 - number_of_tables), x_max);

  MultiTableImplementation<max_table_spacing, number_of_tables> implementation;
  implementation.Initialize();

  Argument a[number_of_iterations];
  for (std::int64_t i = 0; i < number_of_iterations; ++i) {
    a[i] = uniformly_at(random);
  }

  Value v[number_of_iterations];
  while (state.KeepRunningBatch(number_of_iterations)) {
    for (std::int64_t i = 0; i < number_of_iterations; ++i) {
      using namespace principia::quantities;
      v[i] = implementation.Cos(a[i]);
#if _DEBUG
      // The implementation is not accurate, but let's check that it's not
      // broken.
      auto const absolute_error = Abs(v[i] - std::cos(a[i]));
      CHECK_LT(absolute_error, 1.2e-16);
#endif
    }
    benchmark::DoNotOptimize(v);
  }
}

BENCHMARK_TEMPLATE(BM_ExperimentSinTableSpacing, 2.0 / 256.0)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_ExperimentSinTableSpacing, 2.0 / 1024.0)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_ExperimentCosTableSpacing, 2.0 / 256.0)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_ExperimentCosTableSpacing, 2.0 / 1024.0)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_ExperimentSinMultiTable, 2.0 / 1024.0, 9)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_ExperimentCosMultiTable, 2.0 / 256.0, 9)
    ->Unit(benchmark::kNanosecond);

}  // namespace functions
}  // namespace principia
