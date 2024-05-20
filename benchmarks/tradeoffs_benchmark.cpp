#include <array>
#include <cmath>
#include <cstdint>
#include <random>

#include "benchmark/benchmark.h"
#include "numerics/double_precision.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"

namespace principia {
namespace functions {

using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_elementary_functions;

using Value = double;
using Argument = double;

constexpr Argument x_min = π / 12;
constexpr Argument x_max = π / 6;
static constexpr std::int64_t number_of_iterations = 1000;

template<Argument table_spacing>
class Implementation {
 public:
  void Initialize();

  Value Sin(Argument x);
  Value Cos(Argument x);

 private:
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

template<Argument table_spacing>
void Implementation<table_spacing>::Initialize() {
  using namespace principia::quantities;
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
Value Implementation<table_spacing>::Sin(Argument const x) {
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
Value Implementation<table_spacing>::Cos(Argument const x) {
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

//TODO(phl):polynomials

template<Argument table_spacing>
Value Implementation<table_spacing>::SinPolynomial(Argument const x) {
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
Value Implementation<table_spacing>::CosPolynomial(Argument const x) {
  if constexpr (table_spacing == 2.0 / 1024.0) {
    // 72 bits.
    return -0.499999999999999999999872434553 +
           0.0416666654823785864634569932662 * x;
  } else if constexpr (table_spacing == 2.0 / 256.0) {
    // 77 bits.
    return -0.499999999999999999999999974543 +
           x * (0.0416666666666633318024480868405 -
                0.00138888829905860875255146938745 * x);
  }
}

template<Argument table_spacing>
void BM_EvaluateSinTradeOffs(benchmark::State& state) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(x_min, x_max);

  Implementation<table_spacing> implementation;
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
      auto const absolute_error = Abs(v[i] - std::sin(a[i]));
      CHECK_LT(absolute_error, 5.6e-17);
#endif
    }
    benchmark::DoNotOptimize(v);
  }
}

template<Argument table_spacing>
void BM_EvaluateCosTradeOffs(benchmark::State& state) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(x_min, x_max);

  Implementation<table_spacing> implementation;
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
      auto const absolute_error = Abs(v[i] - std::cos(a[i]));
      CHECK_LT(absolute_error, 1.2e-16);
#endif
    }
    benchmark::DoNotOptimize(v);
  }
}

BENCHMARK_TEMPLATE(BM_EvaluateSinTradeOffs, 1.0 / 128.0)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateSinTradeOffs, 1.0 / 512.0)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateCosTradeOffs, 1.0 / 128.0)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluateCosTradeOffs, 1.0 / 512.0)
    ->Unit(benchmark::kNanosecond);

}  // namespace functions
}  // namespace principia
