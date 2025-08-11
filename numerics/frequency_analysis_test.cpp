#include "numerics/frequency_analysis.hpp"

#include <algorithm>
#include <functional>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/apodization.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/fast_fourier_transform.hpp"
#include "numerics/piecewise_poisson_series.hpp"
#include "numerics/poisson_series.hpp"
#include "numerics/polynomial_evaluators.hpp"
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

using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Gt;
using ::testing::Lt;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::numerics::_apodization;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_fast_fourier_transform;
using namespace principia::numerics::_frequency_analysis;
using namespace principia::numerics::_piecewise_poisson_series;
using namespace principia::numerics::_poisson_series;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics_matchers;

// Constructs a piecewise Poisson series that has the given number of pieces
// covering [t_min, t_max] and that matches `series` over that interval.
template<typename Piecewise>
Piecewise Slice(typename Piecewise::Series const& series,
                int const pieces,
                Instant const& t_min,
                Instant const& t_max) {
  Time const Δt = (t_max - t_min) / pieces;
  Piecewise piecewise({t_min, t_min + Δt}, series);
  for (int i = 1; i < pieces; ++i) {
    piecewise.Append({t_min + i * Δt, t_min + (i + 1) * Δt}, series);
  }
  return piecewise;
}

class FrequencyAnalysisTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  using Series0 = PoissonSeries<Length, 0, 0>;
  using Series4 = PoissonSeries<Length, 4, 4>;
  using Polynomial4 = Series4::AperiodicPolynomial;

  FrequencyAnalysisTest()
      : random_polynomial4_([](Instant const& origin,
                               std::mt19937_64& random,
                               std::uniform_real_distribution<>& distribution) {
          auto const c0 = distribution(random) * Metre;
          auto const c1 = distribution(random) * Metre / Second;
          auto const c2 = distribution(random) * Metre / Pow<2>(Second);
          auto const c3 = distribution(random) * Metre / Pow<3>(Second);
          auto const c4 = distribution(random) * Metre / Pow<4>(Second);

          return Polynomial4({c0, c1, c2, c3, c4}, origin);
        }) {}

  Instant const t0_;
  std::function<Polynomial4(
      Instant const& origin,
      std::mt19937_64& random,
      std::uniform_real_distribution<>& distribution)>
      random_polynomial4_;
};

TEST_F(FrequencyAnalysisTest, PreciseModeScalar) {
  using FFT = FastFourierTransform<Length, Instant, 1 << 16>;
  AngularFrequency const ω = 666.543 * π / FFT::size * Radian / Second;
  Time const Δt = 1 * Second;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> amplitude_distribution(-0.1, 0.1);
  std::uniform_real_distribution<> frequency_distribution(-100.0, 100.0);

  using PiecewiseSeries0 = PiecewisePoissonSeries<Length, 0, 0>;
  using Series0 = PiecewiseSeries0::Series;
  Series0::PolynomialsByAngularFrequency polynomials;

  // Main harmonic.
  polynomials.emplace_back(
      ω,
      Series0::Polynomials{
          .sin = Series0::PeriodicPolynomial({1 * Metre}, t0_),
          .cos = Series0::PeriodicPolynomial({0 * Metre}, t0_)});

  // Noise with lower amplitude and higher frequency.
  for (int i = 0; i < 10; ++i) {
    auto const sin_amplitude = amplitude_distribution(random) * Metre;
    auto const cos_amplitude = amplitude_distribution(random) * Metre;
    polynomials.emplace_back(
        ω * frequency_distribution(random),
        Series0::Polynomials{
            .sin = Series0::PeriodicPolynomial({sin_amplitude}, t0_),
            .cos = Series0::PeriodicPolynomial({cos_amplitude}, t0_)});
  }
  Series0 const sin(Series0::AperiodicPolynomial(
                        {amplitude_distribution(random) * Metre}, t0_),
                    polynomials);

  Instant const t_min = t0_;
  Instant const t_max = t0_ + (FFT::size - 1) * Δt;
  PiecewiseSeries0 const piecewise_sin =
      Slice<PiecewiseSeries0>(sin, /*pieces=*/1000, t_min, t_max);

  std::vector<Length> signal;
  for (int n = 0; n < FFT::size; ++n) {
    signal.push_back(piecewise_sin(t0_ + n * Δt));
  }

  // Won't fit on the stack.
  auto transform = std::make_unique<FFT>(signal, Δt);

  // The FFT gives us an accuracy which is of the order of the number of points.
  auto const mode = transform->Mode(AngularFrequency{},
                                    Infinity<AngularFrequency>);
  EXPECT_THAT(mode.midpoint(), RelativeErrorFrom(ω, IsNear(8.1e-4_(1))));

  // The precise analysis is only limited by our ability to pinpoint the
  // maximum.
  auto const precise_mode =
      PreciseMode(mode,
                  piecewise_sin,
                  _apodization::Hann(t_min, t_max));
  EXPECT_THAT(precise_mode, RelativeErrorFrom(ω, IsNear(2.6e-8_(1))));
}

TEST_F(FrequencyAnalysisTest, PreciseModeVector) {
  using FFT = FastFourierTransform<Displacement<World>, Instant, 1 << 16>;
  AngularFrequency const ω = 666.543 * π / FFT::size * Radian / Second;
  Time const Δt = 1 * Second;

  using PiecewiseSeries0 =
      PiecewisePoissonSeries<Displacement<World>, 0, 0>;
  using Series0 = PiecewiseSeries0::Series;
  Series0::PolynomialsByAngularFrequency polynomials;

  // Main harmonic.
  polynomials.emplace_back(
      ω,
      Series0::Polynomials{
          .sin = Series0::PeriodicPolynomial(
              {Displacement<World>({1 * Metre, 2 * Metre, 3 * Metre})}, t0_),
          .cos = Series0::PeriodicPolynomial(
              {Displacement<World>({-5 * Metre, 7 * Metre, 11 * Metre})},
              t0_)});
  Series0 const sin(Series0::AperiodicPolynomial(Displacement<World>(), t0_),
                    polynomials);

  Instant const t_min = t0_;
  Instant const t_max = t0_ + (FFT::size - 1) * Δt;
  PiecewiseSeries0 const piecewise_sin =
      Slice<PiecewiseSeries0>(sin, /*pieces=*/1000, t_min, t_max);

  std::vector<Displacement<World>> signal;
  for (int n = 0; n < FFT::size; ++n) {
    signal.push_back(piecewise_sin(t0_ + n * Δt));
  }

  // Won't fit on the stack.
  auto transform = std::make_unique<FFT>(signal, Δt);

  // The FFT gives us an accuracy which is of the order of the number of points.
  auto const mode = transform->Mode(AngularFrequency{},
                                    Infinity<AngularFrequency>);
  EXPECT_THAT(mode.midpoint(), RelativeErrorFrom(ω, IsNear(8.1e-4_(1))));

  // The precise analysis is only limited by our ability to pinpoint the
  // maximum.
  auto const precise_mode =
      PreciseMode(mode,
                  piecewise_sin,
                  _apodization::Hann(t_min, t_max));
  EXPECT_THAT(precise_mode, RelativeErrorFrom(ω, IsNear(4.2e-11_(1))));
}

TEST_F(FrequencyAnalysisTest, PoissonSeriesScalarProjection) {
  AngularFrequency const ω = 666.543 * π * Radian / Second;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> amplitude_distribution(-10.0, 10.0);

  Instant const t_min = t0_;
  Instant const t_mid = t0_ + 50 * Radian / ω;
  Instant const t_max = t0_ + 100 * Radian / ω;

  auto const sin = random_polynomial4_(t_mid, random, amplitude_distribution);
  auto const cos = random_polynomial4_(t_mid, random, amplitude_distribution);
  Series4 const series(
      Series4::AperiodicPolynomial({}, t_mid),
      {{ω, Series4::Polynomials{sin, cos}}});

  // Projection on a 4th degree basis accurately reconstructs the function.
  auto const projection4 =
      Projection<4, 4>(series,
                       ω,
                       _apodization::Hann(t_min, t_max),
                       t_min, t_max);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(projection4(t_min + i * Radian / ω),
                AlmostEquals(series(t_min + i * Radian / ω), 0, 1536));
  }

  // Projection on a 5th degree basis is also accurate.
  auto const projection5 =
      Projection<5, 5>(series,
                       ω,
                       _apodization::Hann(t_min, t_max),
                       t_min, t_max);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(projection5(t_min + i * Radian / ω),
                AlmostEquals(series(t_min + i * Radian / ω), 0, 1536));
  }

  // Projection on a 3rd degree basis introduces significant errors.
  auto const projection3 =
      Projection<3, 3>(series,
                       ω,
                       _apodization::Hann(t_min, t_max),
                       t_min, t_max);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(projection3(t_min + i * Radian / ω),
                RelativeErrorFrom(series(t_min + i * Radian / ω),
                                  AllOf(Gt(3.6e-13), Lt(9.0e-6))));
  }
}

TEST_F(FrequencyAnalysisTest, PoissonSeriesVectorProjection) {
  AngularFrequency const ω = 666.543 * π * Radian / Second;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> amplitude_distribution(-10.0, 10.0);
  using VectorSeries4 = PoissonSeries<Vector<Length, World>, 4, 4>;
  using Polynomial4 = VectorSeries4::AperiodicPolynomial;

  Instant const t_min = t0_;
  Instant const t_mid = t0_ + 50 * Radian / ω;
  Instant const t_max = t0_ + 100 * Radian / ω;

  auto random_polynomial4 = [](Instant const& origin,
                               std::mt19937_64& random,
                               std::uniform_real_distribution<>& distribution) {
    auto const c0x = distribution(random) * Metre;
    auto const c1x = distribution(random) * Metre / Second;
    auto const c2x = distribution(random) * Metre / Pow<2>(Second);
    auto const c3x = distribution(random) * Metre / Pow<3>(Second);
    auto const c4x = distribution(random) * Metre / Pow<4>(Second);
    auto const c0y = distribution(random) * Metre;
    auto const c1y = distribution(random) * Metre / Second;
    auto const c2y = distribution(random) * Metre / Pow<2>(Second);
    auto const c3y = distribution(random) * Metre / Pow<3>(Second);
    auto const c4y = distribution(random) * Metre / Pow<4>(Second);
    auto const c0z = distribution(random) * Metre;
    auto const c1z = distribution(random) * Metre / Second;
    auto const c2z = distribution(random) * Metre / Pow<2>(Second);
    auto const c3z = distribution(random) * Metre / Pow<3>(Second);
    auto const c4z = distribution(random) * Metre / Pow<4>(Second);
    Vector<Length, World> const v0({c0x, c0y, c0z});
    Vector<Speed, World> const v1({c1x, c1y, c1z});
    Vector<Acceleration, World> const v2({c2x, c2y, c2z});
    Vector<Jerk, World> const v3({c3x, c3y, c3z});
    Vector<Snap, World> const v4({c4x, c4y, c4z});

    return Polynomial4({v0, v1, v2, v3, v4}, origin);
  };

  auto const sin = random_polynomial4(t_mid, random, amplitude_distribution);
  auto const cos = random_polynomial4(t_mid, random, amplitude_distribution);
  VectorSeries4 const series(
      VectorSeries4::AperiodicPolynomial({}, t_mid),
      {{ω, VectorSeries4::Polynomials{sin, cos}}});

  // Projection on a 4th degree basis accurately reconstructs the function.
  auto const projection4 =
      Projection<4, 4>(series,
                       ω,
                       _apodization::Hann(t_min, t_max),
                       t_min, t_max);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(projection4(t_min + i * Radian / ω),
                AlmostEquals(series(t_min + i * Radian / ω), 0, 5120));
  }

  // Projection on a 5th degree basis is also accurate.
  auto const projection5 =
      Projection<5, 5>(series,
                       ω,
                       _apodization::Hann(t_min, t_max),
                       t_min, t_max);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(projection5(t_min + i * Radian / ω),
                AlmostEquals(series(t_min + i * Radian / ω), 0, 5120));
  }

  // Projection on a 3rd degree basis introduces significant errors.
  auto const projection3 =
      Projection<3, 3>(series,
                       ω,
                       _apodization::Hann(t_min, t_max),
                       t_min, t_max);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(projection3(t_min + i * Radian / ω),
                RelativeErrorFrom(series(t_min + i * Radian / ω),
                                  AllOf(Gt(7.3e-11), Lt(2.7e-7))));
  }
}

TEST_F(FrequencyAnalysisTest, PiecewisePoissonSeriesProjection) {
  AngularFrequency const ω = 0.0566543 * π * Radian / Second;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> amplitude_distribution(-10.0, 10.0);
  std::uniform_real_distribution<> perturbation_distribution(-1e-6, 1e-6);

  Instant const t_min = t0_;
  Instant const t_mid = t0_ + 5 * Second;
  Instant const t_max = t0_ + 10 * Second;

  using PiecewiseSeries4 = PiecewisePoissonSeries<Length, 4, 4>;

  auto const sin = random_polynomial4_(t_mid, random, amplitude_distribution);
  auto const cos = random_polynomial4_(t_mid, random, amplitude_distribution);
  Series4 const series(
      Series4::AperiodicPolynomial({}, t_mid),
      {{ω, Series4::Polynomials{sin, cos}}});

  // Build a series that is based on `series` with different perturbations over
  // different intervals.
  PiecewiseSeries4 piecewise_series({t_min, t_min + 1 * Second}, series);
  for (int i = 1; i < 10; ++i) {
    auto const perturbation_sin =
        random_polynomial4_(t_mid, random, perturbation_distribution);
    auto const perturbation_cos =
        random_polynomial4_(t_mid, random, perturbation_distribution);
    Series4 const perturbation_series(
        Series4::AperiodicPolynomial({}, t_mid),
        {{ω, Series4::Polynomials{perturbation_sin, perturbation_cos}}});
    piecewise_series.Append({t_min + i * Second, t_min + (i + 1) * Second},
                            series + perturbation_series);
  }

  // Projection on a 4th degree basis.  The approximation is reasonably
  // accurate.
  auto const projection4 =
      Projection<4, 4>(piecewise_series,
                       ω,
                       _apodization::Dirichlet(t_min, t_max),
                       t_min, t_max);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(
        projection4(t_min + i * (t_max - t_min) / 100),
        RelativeErrorFrom(series(t_min + i * (t_max - t_min) / 100),
                          AllOf(Gt(2.3e-10), Lt(9.9e-5))));
  }
}

#if !defined(_DEBUG)

TEST_F(FrequencyAnalysisTest, PoissonSeriesIncrementalProjectionNoSecular) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> frequency_distribution(2000.0, 3000.0);

  std::vector<AngularFrequency> ωs;
  for (int i = 0; i < 3; ++i) {
    ωs.push_back(frequency_distribution(random) * Radian / Second);
  }

  Instant const t_min = t0_;
  Instant const t_mid =
      t0_ + 100 * Radian / *std::max_element(ωs.cbegin(), ωs.cend());
  Instant const t_max =
      t0_ + 200 * Radian / *std::max_element(ωs.cbegin(), ωs.cend());

  std::optional<Series4> series;
  for (int i = 0; i < 3; ++i) {
    std::uniform_real_distribution<> amplitude_distribution(-(8 >> i),
                                                            (8 >> i));
    auto const sin = random_polynomial4_(t_mid, random, amplitude_distribution);
    auto const cos = random_polynomial4_(t_mid, random, amplitude_distribution);
    Series4 const s(
        Series4::AperiodicPolynomial({}, t_mid),
        {{ωs[i], Series4::Polynomials{sin, cos}}});
    if (series.has_value()) {
      series.value() += s;
    } else {
      series = s;
    }
  }

  // A perfect calculator for the frequencies of the series.
  int ω_index = 0;
  auto angular_frequency_calculator =
      [t_min, t_max, &ω_index, &ωs](
          auto const& residual) -> std::optional<AngularFrequency> {
    for (int i = 0; i <= 100; ++i) {
      EXPECT_THAT(
          Abs(residual(t_min + i * (t_max - t_min) / 100)),
          ω_index == 0
              ? AllOf(Ge(1.3e-1 * Metre), Lt(12 * Metre))
              : ω_index == 1
                    ? AllOf(Ge(4.8e-3 * Metre), Lt(7.8 * Metre))
                    : ω_index == 2
                          ? AllOf(Ge(5.6e-14 * Metre), Lt(6.2e-10 * Metre))
                          : AllOf(Ge(0 * Metre), Lt(1.1e-13 * Metre)))
          << ω_index;
    }
    if (ω_index == ωs.size()) {
      return std::nullopt;
    } else {
      return ωs[ω_index++];
    }
  };

  // Projection on a 4th degree basis reconstructs the function with a decent
  // accuracy.
  auto const projection4 =
      IncrementalProjection<4, 4>(
          series.value(),
          angular_frequency_calculator,
          _apodization::Hann(t_min, t_max),
          t_min, t_max);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(
        projection4(t_min + i * (t_max - t_min) / 100),
        RelativeErrorFrom(series.value()(t_min + i * (t_max - t_min) / 100),
                          AllOf(Ge(0), Lt(1.4e-13))));
  }
}

TEST_F(FrequencyAnalysisTest, PoissonSeriesIncrementalProjectionSecular) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> frequency_distribution(2000.0, 3000.0);
  std::uniform_real_distribution<> secular_distribution(-30.0, 30.0);

  std::vector<AngularFrequency> ωs = {AngularFrequency{}};
  for (int i = 1; i < 4; ++i) {
    ωs.push_back(frequency_distribution(random) * Radian / Second);
  }

  Instant const t_min = t0_;
  Instant const t_mid =
      t0_ + 100 * Radian / *std::max_element(ωs.cbegin(), ωs.cend());
  Instant const t_max =
      t0_ + 200 * Radian / *std::max_element(ωs.cbegin(), ωs.cend());

  Series4 series(
      random_polynomial4_(t_mid, random, secular_distribution), {});

  // Skip the zero frequency in the loop.
  for (int i = 1; i < 4; ++i) {
    std::uniform_real_distribution<> amplitude_distribution(-(16 >> i),
                                                            (16 >> i));
    auto const sin = random_polynomial4_(t_mid, random, amplitude_distribution);
    auto const cos = random_polynomial4_(t_mid, random, amplitude_distribution);
    series += Series4(
        Series4::AperiodicPolynomial({}, t_mid),
        {{ωs[i], Series4::Polynomials{sin, cos}}});
  }

  // A perfect calculator for the frequencies of the series.
  int ω_index = 0;
  auto angular_frequency_calculator =
      [t_min, t_max, &ω_index, &ωs](
          auto const& residual) -> std::optional<AngularFrequency> {
    for (int i = 0; i <= 100; ++i) {
      EXPECT_THAT(Abs(residual(t_min + i * (t_max - t_min) / 100)),
                  ω_index == 0
                      ? AllOf(Ge(12 * Metre), Lt(32 * Metre))
                      : ω_index == 1
                            ? AllOf(Ge(5.2e-3 * Metre), Lt(9.0 * Metre))
                            : ω_index == 2
                                  ? AllOf(Ge(7.0e-3 * Metre), Lt(3.1 * Metre))
                                  : ω_index == 3 ? AllOf(Ge(5.9e-15 * Metre),
                                                         Lt(1.7e-10 * Metre))
                                                 : AllOf(Ge(0 * Metre),
                                                         Lt(3.6e-14 * Metre)))
          << ω_index;
    }
    if (ω_index == ωs.size()) {
      return std::nullopt;
    } else {
      return ωs[ω_index++];
    }
  };

  // Projection on a 4th degree basis reconstructs the function with a decent
  // accuracy.
  auto const projection4 =
      IncrementalProjection<4, 4>(
          series,
          angular_frequency_calculator,
          _apodization::Dirichlet(t_min, t_max),
          t_min, t_max);
  for (int i = 0; i <= 100; ++i) {
    EXPECT_THAT(
        projection4(t_min + i * (t_max - t_min) / 100),
        RelativeErrorFrom(series(t_min + i * (t_max - t_min) / 100),
                          AllOf(Ge(0), Lt(1.7e-15))));
  }
}

#endif

}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
