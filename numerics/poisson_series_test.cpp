
#include "numerics/poisson_series.hpp"

#include <functional>
#include <limits>
#include <memory>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
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
using geometry::Vector;
using geometry::Velocity;
using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::Pow;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::EqualsProto;
using testing_utilities::IsNear;
using testing_utilities::VanishesBefore;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::operator""_⑴;

class PoissonSeriesTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  using Degree0 = PoissonSeries<double, 0, HornerEvaluator>;
  using Degree1 = PoissonSeries<double, 1, HornerEvaluator>;

  PoissonSeriesTest()
      : ω0_(0 * Radian / Second),
        ω1_(1 * Radian / Second),
        ω2_(2 * Radian / Second),
        ω3_(-3 * Radian / Second) {
    Degree1::Polynomial pa0({0, 0 / Second}, t0_);
    Degree1::Polynomial psa0({100, 200 / Second}, t0_);
    Degree1::Polynomial pca0({1, 2 / Second}, t0_);
    Degree1::Polynomial pb0({3, 4 / Second}, t0_);

    Degree1::Polynomial psa1({5, 6 / Second}, t0_);
    Degree1::Polynomial pca1({7, 8 / Second}, t0_);
    Degree1::Polynomial psb1({9, 10 / Second}, t0_);
    Degree1::Polynomial pcb1({11, 12 / Second}, t0_);

    Degree1::Polynomial psa2({13, 14 / Second}, t0_);
    Degree1::Polynomial pca2({15, 16 / Second}, t0_);

    Degree1::Polynomial psb3({-17, -18 / Second}, t0_);
    Degree1::Polynomial pcb3({19, 20 / Second}, t0_);

    Degree1::Polynomials psca0{/*sin=*/psa0, /*cos=*/pca0};

    Degree1::Polynomials psca1{/*sin=*/psa1, /*cos=*/pca1};
    Degree1::Polynomials pscb1{/*sin=*/psb1, /*cos=*/pcb1};

    Degree1::Polynomials psca2{/*sin=*/psa2, /*cos=*/pca2};

    Degree1::Polynomials pscb3{/*sin=*/psb3, /*cos=*/pcb3};

    pa_ = std::make_unique<Degree1>(
        pa0,
        Degree1::PolynomialsByAngularFrequency{
            {ω0_, psca0}, {ω1_, psca1}, {ω2_, psca2}});
    pb_ = std::make_unique<Degree1>(
        pb0,
        Degree1::PolynomialsByAngularFrequency{{ω1_, pscb1}, {ω3_, pscb3}});
  }

  Instant const t0_;
  AngularFrequency const ω0_;
  AngularFrequency const ω1_;
  AngularFrequency const ω2_;
  AngularFrequency const ω3_;
  std::unique_ptr<Degree1> pa_;
  std::unique_ptr<Degree1> pb_;
};

TEST_F(PoissonSeriesTest, Evaluate) {
  EXPECT_THAT((*pa_)(t0_ + 1 * Second),
              AlmostEquals(3 + 11 * Sin(1 * Radian) + 15 * Cos(1 * Radian) +
                               27 * Sin(2 * Radian) + 31 * Cos(2 * Radian),
                           0, 1));
  EXPECT_THAT((*pb_)(t0_ + 1 * Second),
              AlmostEquals(7 + 19 * Sin(1 * Radian) + 23 * Cos(1 * Radian) +
                               35 * Sin(3 * Radian) + 39 * Cos(3 * Radian),
                           32));
}

TEST_F(PoissonSeriesTest, VectorSpace) {
  {
    auto const identity = +*pa_;
    EXPECT_THAT(identity(t0_ + 1 * Second),
                AlmostEquals((*pa_)(t0_ + 1 * Second), 0));
  }
  {
    auto const negated = -*pb_;
    EXPECT_THAT(negated(t0_ + 1 * Second),
                AlmostEquals(-(*pb_)(t0_ + 1 * Second), 0));
  }
  {
    auto const sum = *pa_ + *pb_;
    EXPECT_THAT(sum(t0_ + 1 * Second),
                AlmostEquals((*pa_)(t0_ + 1 * Second) +
                                 (*pb_)(t0_ + 1 * Second), 1));
  }
  {
    auto const difference = *pa_ - *pb_;
    EXPECT_THAT(difference(t0_ + 1 * Second),
                AlmostEquals((*pa_)(t0_ + 1 * Second) -
                                 (*pb_)(t0_ + 1 * Second), 0));
  }
  {
    auto const left_product = 3 * *pa_;
    EXPECT_THAT(left_product(t0_ + 1 * Second),
                AlmostEquals(3 * (*pa_)(t0_ + 1 * Second), 1));
  }
  {
    auto const right_product = *pb_ * 4;
    EXPECT_THAT(right_product(t0_ + 1 * Second),
                AlmostEquals((*pb_)(t0_ + 1 * Second) * 4, 0));
  }
  {
    auto const quotient = *pb_ / 1.5;
    EXPECT_THAT(quotient(t0_ + 1 * Second),
                AlmostEquals((*pb_)(t0_ + 1 * Second) / 1.5, 0, 32));
  }
}

TEST_F(PoissonSeriesTest, Algebra) {
  auto const product = *pa_ * *pb_;
  EXPECT_THAT(product(t0_ + 1 * Second),
              AlmostEquals((*pa_)(t0_ + 1 * Second) *
                               (*pb_)(t0_ + 1 * Second), 6, 38));
}

TEST_F(PoissonSeriesTest, AtOrigin) {
  auto const pa_at_origin = pa_->AtOrigin(t0_ + 2 * Second);
  for (int i = -5; i < 5; ++i) {
    Instant const t = t0_ + i * Second;
    EXPECT_THAT(pa_at_origin(t), AlmostEquals((*pa_)(t), 0, 45));
  }

  auto const pb_at_origin = pb_->AtOrigin(t0_ - 7 * Second);
  for (int i = -5; i < 5; ++i) {
    Instant const t = t0_ + i * Second;
    EXPECT_THAT(pb_at_origin(t), AlmostEquals((*pb_)(t), 0, 132));
  }
}

TEST_F(PoissonSeriesTest, PointwiseInnerProduct) {
  using Degree2 = PoissonSeries<Displacement<World>, 2, HornerEvaluator>;
  Degree2::Polynomial::Coefficients const coefficients_a({
      Displacement<World>({0 * Metre,
                            0 * Metre,
                            1 * Metre}),
      Velocity<World>({0 * Metre / Second,
                        1 * Metre / Second,
                        0 * Metre / Second}),
      Vector<Acceleration, World>({1 * Metre / Second / Second,
                                    0 * Metre / Second / Second,
                                    0 * Metre / Second / Second})});
  Degree2::Polynomial::Coefficients const coefficients_b({
      Displacement<World>({0 * Metre,
                           2 * Metre,
                           3 * Metre}),
      Velocity<World>({-1 * Metre / Second,
                       1 * Metre / Second,
                       0 * Metre / Second}),
      Vector<Acceleration, World>({1 * Metre / Second / Second,
                                   1 * Metre / Second / Second,
                                   -2 * Metre / Second / Second})});
  Degree2 const pa(Degree2::Polynomial({coefficients_a}, t0_), {{}});
  Degree2 const pb(Degree2::Polynomial({coefficients_b}, t0_), {{}});

  auto const product = PointwiseInnerProduct(pa, pb);
  EXPECT_THAT(product(t0_ + 1 * Second),
              AlmostEquals(5 * Metre * Metre, 0));
}

TEST_F(PoissonSeriesTest, Fourier) {
  Degree0::Polynomial constant({1.0}, t0_);
  AngularFrequency const ω = 4 * Radian / Second;
  PoissonSeries<Displacement<World>, 0, HornerEvaluator> signal(
      constant * Displacement<World>{},
      {{ω,
        {/*sin=*/constant *
             Displacement<World>({2 * Metre, -3 * Metre, 5 * Metre}),
         /*cos=*/constant *
             Displacement<World>({-7 * Metre, 11 * Metre, -13 * Metre})}}});
  // Slice our signal into segments short enough that one-point Gauss-Legendre
  // (also known as midpoint) does the job.
  constexpr int segments = 100;
  PiecewisePoissonSeries<Displacement<World>, 0, HornerEvaluator> f(
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

TEST_F(PoissonSeriesTest, Primitive) {
  auto const actual_primitive = pb_->Primitive();

  // The primitive was computed using Mathematica.
  auto const expected_primitive = [=](Time const& t) {
    auto const a0 = 3;
    auto const a1 = 4 / Second;
    auto const b0 = 9;
    auto const b1 = 10 / Second;
    auto const c0 = 11;
    auto const c1 = 12 / Second;
    auto const d0 = -17;
    auto const d1 = -18 / Second;
    auto const e0 = 19;
    auto const e1 = 20 / Second;
    return a0 * t + (a1 * t * t) / 2 +
           (c1 * Cos(ω1_ * t) * Radian * Radian) / (ω1_ * ω1_) -
           (b0 * Cos(ω1_ * t) * Radian) / ω1_ -
           (b1 * t * Cos(ω1_ * t) * Radian) / ω1_ +
           (e1 * Cos(ω3_ * t) * Radian * Radian) / (ω3_ * ω3_) -
           (d0 * Cos(ω3_ * t) * Radian) / ω3_ -
           (d1 * t * Cos(ω3_ * t) * Radian) / ω3_ +
           (b1 * Sin(ω1_ * t) * Radian * Radian) / (ω1_ * ω1_) +
           (c0 * Sin(ω1_ * t) * Radian) / ω1_ +
           (c1 * t * Sin(ω1_ * t) * Radian) / ω1_ +
           (d1 * Sin(ω3_ * t) * Radian * Radian) / (ω3_ * ω3_) +
           (e0 * Sin(ω3_ * t) * Radian) / ω3_ +
           (e1 * t * Sin(ω3_ * t) * Radian) / ω3_;
  };

  for (int i = -10; i < 10; ++i) {
    EXPECT_THAT(actual_primitive(t0_ + i * Second),
                AlmostEquals(expected_primitive(i * Second), 0, 6));
  }

  EXPECT_THAT(
      pb_->Integrate(t0_ + 5 * Second, t0_ + 13 * Second),
      AlmostEquals(expected_primitive(13 * Second) -
                   expected_primitive(5 * Second), 1, 2));
}

TEST_F(PoissonSeriesTest, InnerProduct) {
  Instant const t_min = t0_;
  Instant const t_max = t0_ + 3 * Second;
  // Computed using Mathematica.
  EXPECT_THAT(InnerProduct(*pa_,
                           *pb_,
                           apodization::Hann<HornerEvaluator>(t_min, t_max),
                           t_min,
                           t_max),
              AlmostEquals(-381.25522770148542400, 4, 7));
}

TEST_F(PoissonSeriesTest, PoorlyConditionedInnerProduct) {
  using Degree4 = PoissonSeries<Length, 4, HornerEvaluator>;
  using Degree5 = PoissonSeries<Length, 5, HornerEvaluator>;
  Instant const t_min = t0_;
  Instant const t_max = t0_ + 4.77553415434249021e-02 * Second;
  AngularFrequency const ω = 2.09400659210170170e+03 * Radian / Second;
  Degree4 const f(Degree4::Polynomial({}, t0_),
                  {{ω,
                    {/*sin=*/Degree4::Polynomial(
                         {+5.10311065909077932e+00 * Metre,
                          +2.78062787709394854e+00 * Metre / Second,
                          +5.04290401496053242e+00 * Metre / Pow<2>(Second),
                          -7.27454632735125806e+00 * Metre / Pow<3>(Second),
                          +8.06537932856756612e+00 * Metre / Pow<4>(Second)},
                         t0_),
                     /*cos=*/Degree4::Polynomial(
                         {-8.11863376474325804e+00 * Metre,
                          +1.49140608216528037e+00 * Metre / Second,
                          -2.54224601087630298e+00 * Metre / Pow<2>(Second),
                          -4.52251796525658367e+00 * Metre / Pow<3>(Second),
                          -2.19458237171412751e+00 * Metre / Pow<4>(Second)},
                         t0_)}}});
  Degree5 const q(Degree5::Polynomial({}, t0_),
                  {{ω,
                    {/*sin=*/Degree5::Polynomial(
                         {-4.41249783881549433e+01 * Metre,
                          +1.50208859053174347e+04 * Metre / Second,
                          -1.70674564621978020e+06 * Metre / Pow<2>(Second),
                          +8.52015772027946562e+07 * Metre / Pow<3>(Second),
                          -1.92799129073151779e+09 * Metre / Pow<4>(Second),
                          +1.61514557548221931e+10 * Metre / Pow<5>(Second)},
                         t0_),
                     /*cos=*/Degree5::Polynomial(
                         {-1.00752842659088765e-01 * Metre,
                          +2.25402995957193006e+01 * Metre / Second,
                          -1.66819064858902379e+03 * Metre / Pow<2>(Second),
                          +4.98682536071893774e+04 * Metre / Pow<3>(Second),
                          -5.18229522289936838e+05 * Metre / Pow<4>(Second),
                          0 * Metre / Pow<5>(Second)},
                         t0_)}}});

  // The integral is very small compared to the functions, so we end up in the
  // numerical noise, and adding more points would not help much.
  auto const product = InnerProduct(
      f, q, apodization::Hann<HornerEvaluator>(t_min, t_max), t_min, t_max);
  // Exact result obtained using Mathematica.
  EXPECT_THAT(product,
              RelativeErrorFrom(-4.848079980325297e-13 * Metre * Metre,
                                IsNear(0.05_⑴)));
}

TEST_F(PoissonSeriesTest, Output) {
  LOG(ERROR) << *pa_;
}

TEST_F(PoissonSeriesTest, Serialization) {
  serialization::PoissonSeries message;
  pa_->WriteToMessage(&message);
  EXPECT_TRUE(message.has_aperiodic());
  EXPECT_EQ(2, message.periodic_size());

  auto const poisson_series_read = Degree1::ReadFromMessage(message);
  EXPECT_THAT((*pa_)(t0_ + 1 * Second),
              AlmostEquals(poisson_series_read(t0_ + 1 * Second), 0));
  EXPECT_THAT((*pa_)(t0_ + 2 * Second),
              AlmostEquals(poisson_series_read(t0_ + 2 * Second), 0));
  EXPECT_THAT((*pa_)(t0_ + 3 * Second),
              AlmostEquals(poisson_series_read(t0_ + 3 * Second), 0));

  serialization::PoissonSeries message2;
  poisson_series_read.WriteToMessage(&message2);
  EXPECT_THAT(message2, EqualsProto(message));
}

class PiecewisePoissonSeriesTest : public ::testing::Test {
 protected:
  using Degree0 = PiecewisePoissonSeries<double, 0, HornerEvaluator>;

  PiecewisePoissonSeriesTest()
      : ω_(π / 2 * Radian / Second),
        p_(Degree0::Series::Polynomial({1.5}, t0_),
           {{ω_,
             {/*sin=*/Degree0::Series::Polynomial({0.5}, t0_),
              /*cos=*/Degree0::Series::Polynomial({-1}, t0_)}}}),
        pp_({t0_, t0_ + 1 * Second},
            Degree0::Series(
                Degree0::Series::Polynomial({1}, t0_),
                {{ω_,
                  {/*sin=*/Degree0::Series::Polynomial({-1}, t0_),
                   /*cos=*/Degree0::Series::Polynomial({0}, t0_)}}})) {
    pp_.Append(
        {t0_ + 1 * Second, t0_ + 2 * Second},
        Degree0::Series(Degree0::Series::Polynomial({0}, t0_),
                        {{ω_,
                          {/*sin=*/Degree0::Series::Polynomial({0}, t0_),
                           /*cos=*/Degree0::Series::Polynomial({1}, t0_)}}}));
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
                AlmostEquals((6 + Sqrt(2)) / 4, 0));
    EXPECT_THAT(s2(t0_ + 0.5 * Second),
                AlmostEquals((10 - 3 * Sqrt(2)) / 4, 0));
    EXPECT_THAT(s2(t0_ + 1.5 * Second),
                AlmostEquals((6 + Sqrt(2)) / 4, 0));
  }
  {
    auto const d1 = p_ - pp_;
    auto const d2 = pp_ - p_;
    EXPECT_THAT(d1(t0_ + 0.5 * Second),
                AlmostEquals((2 + Sqrt(2)) / 4, 1));
    EXPECT_THAT(d1(t0_ + 1.5 * Second),
                AlmostEquals((6 + 5 * Sqrt(2)) / 4, 0));
    EXPECT_THAT(d2(t0_ + 0.5 * Second),
                AlmostEquals((-2 - Sqrt(2)) / 4, 1));
    EXPECT_THAT(d2(t0_ + 1.5 * Second),
                AlmostEquals((-6 - 5 * Sqrt(2)) / 4, 0));
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
                AlmostEquals((10 - 3 * Sqrt(2)) / 4, 1));
    EXPECT_THAT(s1(t0_ + 1.5 * Second),
                AlmostEquals((6 + Sqrt(2)) / 4, 0));
    EXPECT_THAT(s2(t0_ + 0.5 * Second),
                AlmostEquals((10 - 3 * Sqrt(2)) / 4, 1));
    EXPECT_THAT(s2(t0_ + 1.5 * Second),
                AlmostEquals((6 + Sqrt(2)) / 4, 0));
  }
  {
    auto const d1 = p - pp_;
    auto const d2 = pp_ - p;
    EXPECT_THAT(d1(t0_ + 0.5 * Second),
                AlmostEquals((2 + Sqrt(2)) / 4, 0, 2));
    EXPECT_THAT(d1(t0_ + 1.5 * Second),
                AlmostEquals((6 + 5 * Sqrt(2)) / 4, 0));
    EXPECT_THAT(d2(t0_ + 0.5 * Second),
                AlmostEquals((-2 - Sqrt(2)) / 4, 0, 2));
    EXPECT_THAT(d2(t0_ + 1.5 * Second),
                AlmostEquals((-6 - 5 * Sqrt(2)) / 4, 0));
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
  double const d1 = InnerProduct<double, double, 0, 0, 0, HornerEvaluator, 8>(
      pp_, p_, apodization::Dirichlet<HornerEvaluator>(t0_, t0_ + 2 * Second));
  double const d2 = InnerProduct<double, double, 0, 0, 0, HornerEvaluator, 8>(
      p_, pp_, apodization::Dirichlet<HornerEvaluator>(t0_, t0_ + 2 * Second));
  EXPECT_THAT(d1, AlmostEquals((3 * π - 26) / (8 * π), 0));
  EXPECT_THAT(d2, AlmostEquals((3 * π - 26) / (8 * π), 0));
}

TEST_F(PiecewisePoissonSeriesTest, InnerProductMultiorigin) {
  auto const p = p_.AtOrigin(t0_ + 2 * Second);
  double const d1 = InnerProduct<double, double, 0, 0, 0, HornerEvaluator, 8>(
      pp_, p, apodization::Dirichlet<HornerEvaluator>(t0_, t0_ + 2 * Second));
  double const d2 = InnerProduct<double, double, 0, 0, 0, HornerEvaluator, 8>(
      p, pp_, apodization::Dirichlet<HornerEvaluator>(t0_, t0_ + 2 * Second));
  EXPECT_THAT(d1, AlmostEquals((3 * π - 26) / (8 * π), 0));
  EXPECT_THAT(d2, AlmostEquals((3 * π - 26) / (8 * π), 0));
}

TEST_F(PiecewisePoissonSeriesTest, Serialization) {
  serialization::PiecewisePoissonSeries message;
  pp_.WriteToMessage(&message);
  EXPECT_EQ(3, message.bounds_size());
  EXPECT_EQ(2, message.series_size());

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
