
#include "numerics/piecewise_poisson_series.hpp"

#include <functional>
#include <limits>
#include <memory>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "numerics/apodization.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/quadrature.hpp"
#include "numerics/root_finders.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/numerics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics_matchers.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace numerics {

using geometry::Displacement;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::Instant;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Sin;
using quantities::Sqrt;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::EqualsProto;
using testing_utilities::IsNear;
using testing_utilities::VanishesBefore;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::operator""_⑴;

class PiecewisePoissonSeriesTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  using Degree0 = PiecewisePoissonSeries<double, 0, 0, HornerEvaluator>;

  PiecewisePoissonSeriesTest()
      : ω_(π / 2 * Radian / Second),
        p_(Degree0::Series::AperiodicPolynomial({1.5}, t0_),
           {{ω_,
             {.sin = Degree0::Series::PeriodicPolynomial({0.5}, t0_),
              .cos = Degree0::Series::PeriodicPolynomial({-1}, t0_)}}}),
        pp_({t0_, t0_ + 1 * Second},
            Degree0::Series(
                Degree0::Series::AperiodicPolynomial({1}, t0_),
                {{ω_,
                  {.sin = Degree0::Series::PeriodicPolynomial({-1}, t0_),
                   .cos = Degree0::Series::PeriodicPolynomial({0}, t0_)}}})) {
    pp_.Append(
        {t0_ + 1 * Second, t0_ + 2 * Second},
        Degree0::Series(
            Degree0::Series::AperiodicPolynomial({0}, t0_),
            {{ω_,
              {.sin = Degree0::Series::PeriodicPolynomial({0}, t0_),
               .cos = Degree0::Series::PeriodicPolynomial({1}, t0_)}}}));
  }

  Instant const t0_;
  AngularFrequency const ω_;
  // p[t_, t0_] := 3/2 - Cos[π(t - t0)/2] + 1/2 Sin[π(t - t0)/2]
  Degree0::Series p_;
  // pp[t_, t0_] := If[t < t0 + 1, 1 - Sin[π(t - t0)/2], Cos[π(t - t0)/2]]
  Degree0 pp_;
};

TEST_F(PiecewisePoissonSeriesTest, Evaluate) {
  double const ε = std::numeric_limits<double>::epsilon();
  EXPECT_THAT(pp_(t0_), AlmostEquals(1, 0));
  EXPECT_THAT(pp_(t0_ + 0.5 * Second), AlmostEquals(1 - Sqrt(0.5), 0, 2));
  EXPECT_THAT(pp_(t0_ + 1 * (1 - ε / 2) * Second), AlmostEquals(0, 0));
  EXPECT_THAT(pp_(t0_ + 1 * Second), VanishesBefore(1, 0));
  EXPECT_THAT(pp_(t0_ + 1 * (1 + ε) * Second), VanishesBefore(1, 3));
  EXPECT_THAT(pp_(t0_ + 1.5 * Second), AlmostEquals(-Sqrt(0.5), 1));
  EXPECT_THAT(pp_(t0_ + 2 * (1 - ε / 2) * Second), AlmostEquals(-1, 0));
  EXPECT_THAT(pp_(t0_ + 2 * Second), AlmostEquals(-1, 0));
}

TEST_F(PiecewisePoissonSeriesTest, VectorSpace) {
  {
    auto const pp = +pp_;
    EXPECT_THAT(pp(t0_ + 0.5 * Second), AlmostEquals(1 - Sqrt(0.5), 0, 2));
    EXPECT_THAT(pp(t0_ + 1.5 * Second), AlmostEquals(-Sqrt(0.5), 1));
  }
  {
    auto const pp = -pp_;
    EXPECT_THAT(pp(t0_ + 0.5 * Second), AlmostEquals(-1 + Sqrt(0.5), 0, 2));
    EXPECT_THAT(pp(t0_ + 1.5 * Second), AlmostEquals(Sqrt(0.5), 1));
  }
  {
    auto const pp = 2 * pp_;
    EXPECT_THAT(pp(t0_ + 0.5 * Second), AlmostEquals(2 - Sqrt(2), 0, 2));
    EXPECT_THAT(pp(t0_ + 1.5 * Second), AlmostEquals(-Sqrt(2), 1));
  }
  {
    auto const pp = pp_ * 3;
    EXPECT_THAT(pp(t0_ + 0.5 * Second), AlmostEquals(3 - 3 * Sqrt(0.5), 0, 4));
    EXPECT_THAT(pp(t0_ + 1.5 * Second), AlmostEquals(-3 * Sqrt(0.5), 1));
  }
  {
    auto const pp = pp_ / 4;
    EXPECT_THAT(pp(t0_ + 0.5 * Second), AlmostEquals((2 - Sqrt(2)) / 8, 0, 2));
    EXPECT_THAT(pp(t0_ + 1.5 * Second), AlmostEquals(-Sqrt(0.5) / 4, 1));
  }
}

TEST_F(PiecewisePoissonSeriesTest, Action) {
  {
    auto const s1 = p_ + pp_;
    auto const s2 = pp_ + p_;
    EXPECT_THAT(s1(t0_ + 0.5 * Second),
                AlmostEquals((10 - 3 * Sqrt(2)) / 4, 0));
    EXPECT_THAT(s1(t0_ + 1.5 * Second),
                AlmostEquals((6 + Sqrt(2)) / 4, 1));
    EXPECT_THAT(s2(t0_ + 0.5 * Second),
                AlmostEquals((10 - 3 * Sqrt(2)) / 4, 0));
    EXPECT_THAT(s2(t0_ + 1.5 * Second),
                AlmostEquals((6 + Sqrt(2)) / 4, 1));
  }
  {
    auto const d1 = p_ - pp_;
    auto const d2 = pp_ - p_;
    EXPECT_THAT(d1(t0_ + 0.5 * Second),
                AlmostEquals((2 + Sqrt(2)) / 4, 1, 2));
    EXPECT_THAT(d1(t0_ + 1.5 * Second),
                AlmostEquals((6 + 5 * Sqrt(2)) / 4, 1));
    EXPECT_THAT(d2(t0_ + 0.5 * Second),
                AlmostEquals((-2 - Sqrt(2)) / 4, 1, 2));
    EXPECT_THAT(d2(t0_ + 1.5 * Second),
                AlmostEquals((-6 - 5 * Sqrt(2)) / 4, 1));
  }
  {
    auto const p1 = p_ * pp_;
    auto const p2 = pp_ * p_;
    EXPECT_THAT(p1(t0_ + 0.5 * Second),
                AlmostEquals((7 - 4 * Sqrt(2)) / 4, 0, 4));
    EXPECT_THAT(p1(t0_ + 1.5 * Second),
                AlmostEquals((-3 - 3 * Sqrt(2)) / 4, 1));
    EXPECT_THAT(p2(t0_ + 0.5 * Second),
                AlmostEquals((7 - 4 * Sqrt(2)) / 4, 0, 4));
    EXPECT_THAT(p2(t0_ + 1.5 * Second),
                AlmostEquals((-3 - 3 * Sqrt(2)) / 4, 1));
  }
}

TEST_F(PiecewisePoissonSeriesTest, ActionMultiorigin) {
  auto const p = p_.AtOrigin(t0_ + 2 * Second);
  {
    auto const s1 = p + pp_;
    auto const s2 = pp_ + p;
    EXPECT_THAT(s1(t0_ + 0.5 * Second),
                AlmostEquals((10 - 3 * Sqrt(2)) / 4, 0, 1));
    EXPECT_THAT(s1(t0_ + 1.5 * Second),
                AlmostEquals((6 + Sqrt(2)) / 4, 1));
    EXPECT_THAT(s2(t0_ + 0.5 * Second),
                AlmostEquals((10 - 3 * Sqrt(2)) / 4, 0, 1));
    EXPECT_THAT(s2(t0_ + 1.5 * Second),
                AlmostEquals((6 + Sqrt(2)) / 4, 1));
  }
  {
    auto const d1 = p - pp_;
    auto const d2 = pp_ - p;
    EXPECT_THAT(d1(t0_ + 0.5 * Second),
                AlmostEquals((2 + Sqrt(2)) / 4, 0, 2));
    EXPECT_THAT(d1(t0_ + 1.5 * Second),
                AlmostEquals((6 + 5 * Sqrt(2)) / 4, 0, 1));
    EXPECT_THAT(d2(t0_ + 0.5 * Second),
                AlmostEquals((-2 - Sqrt(2)) / 4, 0, 2));
    EXPECT_THAT(d2(t0_ + 1.5 * Second),
                AlmostEquals((-6 - 5 * Sqrt(2)) / 4, 0, 1));
  }
  {
    auto const p1 = p * pp_;
    auto const p2 = pp_ * p;
    EXPECT_THAT(p1(t0_ + 0.5 * Second),
                AlmostEquals((7 - 4 * Sqrt(2)) / 4, 0, 4));
    EXPECT_THAT(p1(t0_ + 1.5 * Second),
                AlmostEquals((-3 - 3 * Sqrt(2)) / 4, 1));
    EXPECT_THAT(p2(t0_ + 0.5 * Second),
                AlmostEquals((7 - 4 * Sqrt(2)) / 4, 0, 4));
    EXPECT_THAT(p2(t0_ + 1.5 * Second),
                AlmostEquals((-3 - 3 * Sqrt(2)) / 4, 1));
  }
}

TEST_F(PiecewisePoissonSeriesTest, InnerProduct) {
  double const d1 = InnerProduct(
      pp_, p_,
      apodization::Dirichlet<HornerEvaluator>(t0_, t0_ + 2 * Second),
      /*max_points=*/1 << 20);
  double const d2 = InnerProduct(
      p_, pp_,
      apodization::Dirichlet<HornerEvaluator>(t0_, t0_ + 2 * Second),
      /*max_points=*/1 << 20);
  EXPECT_THAT(d1, RelativeErrorFrom((3 * π - 26) / (8 * π), IsNear(3e-11_⑴)));
  EXPECT_THAT(d2, RelativeErrorFrom((3 * π - 26) / (8 * π), IsNear(3e-11_⑴)));
}

TEST_F(PiecewisePoissonSeriesTest, InnerProductMultiorigin) {
  auto const p = p_.AtOrigin(t0_ + 2 * Second);
  double const d1 = InnerProduct(
      pp_, p,
      apodization::Dirichlet<HornerEvaluator>(t0_, t0_ + 2 * Second),
      /*max_points=*/1 << 20);
  double const d2 = InnerProduct(
      p, pp_,
      apodization::Dirichlet<HornerEvaluator>(t0_, t0_ + 2 * Second),
      /*max_points=*/1 << 20);
  EXPECT_THAT(d1, RelativeErrorFrom((3 * π - 26) / (8 * π), IsNear(3e-11_⑴)));
  EXPECT_THAT(d2, RelativeErrorFrom((3 * π - 26) / (8 * π), IsNear(3e-11_⑴)));
}

TEST_F(PiecewisePoissonSeriesTest, Fourier) {
  Degree0::Series::AperiodicPolynomial aperiodic_constant({1.0}, t0_);
  Degree0::Series::PeriodicPolynomial periodic_constant({1.0}, t0_);
  AngularFrequency const ω = 4 * Radian / Second;
  PoissonSeries<Displacement<World>, 0, 0, HornerEvaluator> signal(
      aperiodic_constant * Displacement<World>{},
      {{ω,
        {.sin = periodic_constant *
             Displacement<World>({2 * Metre, -3 * Metre, 5 * Metre}),
         .cos = periodic_constant *
             Displacement<World>({-7 * Metre, 11 * Metre, -13 * Metre})}}});
  // Slice our signal into segments short enough that one-point Gauss-Legendre
  // (also known as midpoint) does the job.
  constexpr int segments = 100;
  PiecewisePoissonSeries<Displacement<World>, 0, 0, HornerEvaluator> f(
      {t0_, t0_ + π * Second / segments}, signal);
  for (int k = 1; k < segments; ++k) {
    f.Append({t0_ + k * π * Second / segments,
              t0_ + (k + 1) * π * Second / segments},
             signal);
  }
  auto const f_fourier_transform = f.FourierTransform();
  auto const f_power_spectrum =
      [&f_fourier_transform](AngularFrequency const& ω) {
        return f_fourier_transform(ω).Norm²();
      };
  EXPECT_THAT(Brent(f_power_spectrum, 0.9 * ω, 1.1 * ω, std::greater<>{}),
              RelativeErrorFrom(ω, IsNear(0.03_⑴)));
  // A peak arising from the finite interval.
  EXPECT_THAT(Brent(f_power_spectrum,
                    0 * Radian / Second,
                    2 * Radian / Second,
                    std::greater<>{}),
              IsNear(1.209_⑴ * Radian / Second));

  auto const fw = f * apodization::Hann<HornerEvaluator>(f.t_min(), f.t_max());

  auto const fw_fourier_transform = fw.FourierTransform();
  auto const fw_power_spectrum =
      [&fw_fourier_transform](AngularFrequency const& ω) {
        return fw_fourier_transform(ω).Norm²();
      };
  EXPECT_THAT(Brent(fw_power_spectrum, 0.9 * ω, 1.1 * ω, std::greater<>{}),
              RelativeErrorFrom(ω, IsNear(0.005_⑴)));
  // Hann filters out the interval; this search for a second maximum converges
  // to its boundary.
  EXPECT_THAT(Brent(fw_power_spectrum,
                    0 * Radian / Second,
                    2 * Radian / Second,
                    std::greater<>{}),
              IsNear(1.9999999_⑴ * Radian / Second));
}

TEST_F(PiecewisePoissonSeriesTest, Serialization) {
  serialization::PiecewisePoissonSeries message;
  pp_ += p_;
  pp_.WriteToMessage(&message);
  EXPECT_EQ(3, message.bounds_size());
  EXPECT_EQ(2, message.series_size());
  EXPECT_TRUE(message.has_addend());

  auto const piecewise_poisson_series_read = Degree0::ReadFromMessage(message);
  EXPECT_THAT(
      pp_(t0_ + 0.5 * Second),
      AlmostEquals(piecewise_poisson_series_read(t0_ + 0.5 * Second), 0));
  EXPECT_THAT(
      pp_(t0_ + 1 * Second),
      AlmostEquals(piecewise_poisson_series_read(t0_ + 1 * Second), 0));
  EXPECT_THAT(
      pp_(t0_ + 1.5 * Second),
      AlmostEquals(piecewise_poisson_series_read(t0_ + 1.5 * Second), 0));

  serialization::PiecewisePoissonSeries message2;
  piecewise_poisson_series_read.WriteToMessage(&message2);
  EXPECT_THAT(message2, EqualsProto(message));
}

}  // namespace numerics
}  // namespace principia
