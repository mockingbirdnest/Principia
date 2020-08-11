
#include "numerics/frequency_analysis.hpp"

#include <random>
#include <vector>

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/apodization.hpp"
#include "numerics/fast_fourier_transform.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {

using geometry::Instant;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::Pow;
using quantities::Product;
using quantities::Sin;
using quantities::Square;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::operator""_⑴;
using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;
namespace si = quantities::si;

class DotImplementation {
 public:
  DotImplementation(Instant const& t_min, Instant const& t_max);

  template<typename LFunction, typename RFunction, typename Weight>
  Product<std::invoke_result_t<LFunction, Instant>,
          std::invoke_result_t<RFunction, Instant>>
  operator()(LFunction const& left,
             RFunction const& right,
             Weight const& weight) const;

 private:
  Instant const t_min_;
  Instant const t_max_;
};

DotImplementation::DotImplementation(Instant const& t_min,
                                     Instant const& t_max)
    : t_min_(t_min),
      t_max_(t_max) {}

template<typename LFunction, typename RFunction, typename Weight>
Product<std::invoke_result_t<LFunction, Instant>,
        std::invoke_result_t<RFunction, Instant>>
DotImplementation::operator()(LFunction const& left,
                              RFunction const& right,
                              Weight const& weight) const {
  return Dot(left, right, weight, t_min_, t_max_);
}

class FrequencyAnalysisTest : public ::testing::Test {};

TEST_F(FrequencyAnalysisTest, PreciseMode) {
  using FFT = FastFourierTransform<Length, 1 << 16>;
  Instant const t0;
  AngularFrequency const ω = 666.543 * π / FFT::size * Radian / Second;
  Time const Δt = 1 * Second;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> amplitude_distribution(-0.1, 0.1);
  std::uniform_real_distribution<> frequency_distribution(-100.0, 100.0);

  using Series = PoissonSeries<Length, 0, HornerEvaluator>;
  Series::PolynomialsByAngularFrequency polynomials;

  // Main harmonic.
  polynomials.emplace(
      ω,
      Series::Polynomials{/*sin=*/Series::Polynomial({1 * Metre}, t0),
                          /*cos=*/Series::Polynomial({0 * Metre}, t0)});

  // Noise with lower amplitude and higher frequency.
  for (int i = 0; i < 10; ++i) {
    auto const sin_amplitude = amplitude_distribution(random) * Metre;
    auto const cos_amplitude = amplitude_distribution(random) * Metre;
    polynomials.emplace(
        ω * frequency_distribution(random),
        Series::Polynomials{/*sin=*/Series::Polynomial({sin_amplitude}, t0),
                            /*cos=*/Series::Polynomial({cos_amplitude}, t0)});
  }
  Series const sin(
      Series::Polynomial({amplitude_distribution(random) * Metre}, t0),
      polynomials);

  Instant const t_min = t0;
  Instant const t_max = t0 + (FFT::size - 1) * Δt;
  std::vector<Length> signal;
  for (int n = 0; n < FFT::size; ++n) {
    signal.push_back(sin(t0 + n * Δt));
  }

  // Won't fit on the stack.
  auto transform = std::make_unique<FFT>(signal, Δt);

  // The FFT gives us an accuracy which is of the order of the number of points.
  auto const mode = transform->Mode();
  EXPECT_THAT(mode.midpoint(), RelativeErrorFrom(ω, IsNear(8.1e-4_⑴)));

  DotImplementation dot(t_min, t_max);

  // The precise analysis is only limited by our ability to pinpoint the
  // maximum.
  auto const precise_mode = PreciseMode(
      mode, sin, apodization::Hann<HornerEvaluator>(t_min, t_max), dot);
  EXPECT_THAT(precise_mode, RelativeErrorFrom(ω, IsNear(6.4e-11_⑴)));
}

TEST_F(FrequencyAnalysisTest, PoissonSeriesProjection) {
  Instant const t0;
  AngularFrequency const ω = 666.543 * π * Radian / Second;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> amplitude_distribution(-10.0, 10.0);

  using Series = PoissonSeries<Length, 4, HornerEvaluator>;

  auto random_polynomial = [&amplitude_distribution, &random, &t0]() {
    auto const c0 = amplitude_distribution(random) * Metre;
    auto const c1 = amplitude_distribution(random) * Metre / Second;
    auto const c2 = amplitude_distribution(random) * Metre / Pow<2>(Second);
    auto const c3 = amplitude_distribution(random) * Metre / Pow<3>(Second);
    auto const c4 = amplitude_distribution(random) * Metre / Pow<4>(Second);

    return Series::Polynomial({c0, c1, c2, c3, c4}, t0);
  };

  auto const sin = random_polynomial();
  auto const cos = random_polynomial();
  Series const series(
      Series::Polynomial(Series::Polynomial::Coefficients{}, t0),
      {{ω, Series::Polynomials{sin, cos}}});

  Instant const t_min = t0;
  Instant const t_max = t0 + 100 * Radian / ω;
  DotImplementation const dot(t_min, t_max);

  // Projection on a 4-th degree basis accurately reconstructs the function.
  auto const projection4 = Projection<4>(
      ω, series, apodization::Hann<HornerEvaluator>(t_min, t_max), dot);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(projection4(t0 + i * Radian / ω),
                AlmostEquals(series(t0 + i * Radian / ω), 0, 2688));
  }

  // Projection on a 5-th degree basis is also accurate.
  auto const projection5 = Projection<5>(
      ω, series, apodization::Hann<HornerEvaluator>(t_min, t_max), dot);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(projection5(t0 + i * Radian / ω),
                AlmostEquals(series(t0 + i * Radian / ω), 0, 8000));
  }

  // Projection on a 3-rd degree basis introduces significant errors.
  auto const projection3 = Projection<3>(
      ω, series, apodization::Hann<HornerEvaluator>(t_min, t_max), dot);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(projection3(t0 + i * Radian / ω),
                RelativeErrorFrom(series(t0 + i * Radian / ω),
                                  AllOf(Gt(3.6e-13), Lt(9.0e-6))));
  }
}

TEST_F(FrequencyAnalysisTest, PiecewisePoissonSeriesProjection) {
  Instant const t0;
  AngularFrequency const ω = 666.543 * π * Radian / Second;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> amplitude_distribution(-10.0, 10.0);
  std::uniform_real_distribution<> perturbation_distribution(-1e-6, 1e-6);

  using Series = PoissonSeries<Length, 4, HornerEvaluator>;
  using PiecewiseSeries = PiecewisePoissonSeries<Length, 4, HornerEvaluator>;

  auto random_polynomial =
      [&random, &t0](std::uniform_real_distribution<> const distribution) {
        auto const c0 = distribution(random) * Metre;
        auto const c1 = distribution(random) * Metre / Second;
        auto const c2 = distribution(random) * Metre / Pow<2>(Second);
        auto const c3 = distribution(random) * Metre / Pow<3>(Second);
        auto const c4 = distribution(random) * Metre / Pow<4>(Second);

        return Series::Polynomial({c0, c1, c2, c3, c4}, t0);
      };

  auto const sin = random_polynomial(amplitude_distribution);
  auto const cos = random_polynomial(amplitude_distribution);
  Series const series(
      Series::Polynomial(Series::Polynomial::Coefficients{}, t0),
      {{ω, Series::Polynomials{sin, cos}}});

  // Build a series that is based on |series| with different perturbations over
  // different intervals.
  PiecewiseSeries piecewise_series({t0, t0 + 1 * Second}, series);
  for (int i = 1; i < 3; ++i) {
    auto const perturbation_sin = random_polynomial(perturbation_distribution);
    auto const perturbation_cos = random_polynomial(perturbation_distribution);
    Series const perturbation_series(
        Series::Polynomial(Series::Polynomial::Coefficients{}, t0),
        {{ω, Series::Polynomials{perturbation_sin, perturbation_cos}}});
    piecewise_series.Append({t0 + i * Second, t0 + (i + 1) * Second},
                            series + perturbation_series);
  }

  Instant const t_min = piecewise_series.t_min();
  Instant const t_max = piecewise_series.t_max();
  DotImplementation const dot(t_min, t_max);

  // Projection on a 4-th degree basis.  The errors are of the order of the
  // perturbation.
  auto const projection4 =
      Projection<4>(ω,
                    piecewise_series,
                    apodization::Hann<HornerEvaluator>(t_min, t_max),
                    dot);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(projection4(t0 + i * Radian / ω),
                RelativeErrorFrom(series(t0 + i * Radian / ω),
                                  AllOf(Gt(2.1e-7), Lt(8.8e-4))));
  }
}

}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
