
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
using testing_utilities::operator""_;
using ::testing::AnyOf;

class PoissonSeriesTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  using Degree0 = PoissonSeries<double, 0, 0, HornerEvaluator>;
  using Degree1 = PoissonSeries<double, 1, 1, HornerEvaluator>;

  PoissonSeriesTest()
      : ω0_(0 * Radian / Second),
        ω1_(1 * Radian / Second),
        ω2_(2 * Radian / Second),
        ω3_(-3 * Radian / Second) {
    Degree1::AperiodicPolynomial pa0({0, 0 / Second}, t0_);
    Degree1::PeriodicPolynomial psa0({100, 200 / Second}, t0_);
    Degree1::PeriodicPolynomial pca0({1, 2 / Second}, t0_);
    Degree1::AperiodicPolynomial pb0({3, 4 / Second}, t0_);

    Degree1::PeriodicPolynomial psa1({5, 6 / Second}, t0_);
    Degree1::PeriodicPolynomial pca1({7, 8 / Second}, t0_);
    Degree1::PeriodicPolynomial psb1({9, 10 / Second}, t0_);
    Degree1::PeriodicPolynomial pcb1({11, 12 / Second}, t0_);

    Degree1::PeriodicPolynomial psa2({13, 14 / Second}, t0_);
    Degree1::PeriodicPolynomial pca2({15, 16 / Second}, t0_);

    Degree1::PeriodicPolynomial psb3({-17, -18 / Second}, t0_);
    Degree1::PeriodicPolynomial pcb3({19, 20 / Second}, t0_);

    Degree1::Polynomials psca0{.sin = psa0, .cos = pca0};

    Degree1::Polynomials psca1{.sin = psa1, .cos = pca1};
    Degree1::Polynomials pscb1{.sin = psb1, .cos = pcb1};

    Degree1::Polynomials psca2{.sin = psa2, .cos = pca2};

    Degree1::Polynomials pscb3{.sin = psb3, .cos = pcb3};

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

TEST_F(PoissonSeriesTest, Conversion) {
  using Degree3 = PoissonSeries<double, 3, 3, HornerEvaluator>;
  Degree3 const pa3 = Degree3(*pa_);
  EXPECT_THAT(pa3(t0_ + 1 * Second),
              AlmostEquals(3 + 11 * Sin(1 * Radian) + 15 * Cos(1 * Radian) +
                               27 * Sin(2 * Radian) + 31 * Cos(2 * Radian),
                           0, 1));
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
  using Degree2 = PoissonSeries<Displacement<World>, 2, 0, HornerEvaluator>;
  Degree2::AperiodicPolynomial::Coefficients const coefficients_a({
      Displacement<World>({0 * Metre,
                            0 * Metre,
                            1 * Metre}),
      Velocity<World>({0 * Metre / Second,
                        1 * Metre / Second,
                        0 * Metre / Second}),
      Vector<Acceleration, World>({1 * Metre / Second / Second,
                                    0 * Metre / Second / Second,
                                    0 * Metre / Second / Second})});
  Degree2::AperiodicPolynomial::Coefficients const coefficients_b({
      Displacement<World>({0 * Metre,
                           2 * Metre,
                           3 * Metre}),
      Velocity<World>({-1 * Metre / Second,
                       1 * Metre / Second,
                       0 * Metre / Second}),
      Vector<Acceleration, World>({1 * Metre / Second / Second,
                                   1 * Metre / Second / Second,
                                   -2 * Metre / Second / Second})});
  Degree2 const pa(Degree2::AperiodicPolynomial({coefficients_a}, t0_), {{}});
  Degree2 const pb(Degree2::AperiodicPolynomial({coefficients_b}, t0_), {{}});

  auto const product = PointwiseInnerProduct(pa, pb);
  EXPECT_THAT(product(t0_ + 1 * Second),
              AlmostEquals(5 * Metre * Metre, 0));
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
                   expected_primitive(5 * Second), 0));
}

TEST_F(PoissonSeriesTest, InnerProduct) {
  Instant const t_min = t0_;
  Instant const t_mid = t0_ + 1.5 * Second;
  Instant const t_max = t0_ + 3 * Second;
  // Computed using Mathematica.
  EXPECT_THAT(InnerProduct(pa_->AtOrigin(t_mid),
                           pb_->AtOrigin(t_mid),
                           apodization::Hann<HornerEvaluator>(t_min, t_max),
                           t_min,
                           t_max),
              AlmostEquals(-381.25522770148542400, 0, 7));
}

TEST_F(PoissonSeriesTest, PoorlyConditionedInnerProduct1) {
  using Degree4 = PoissonSeries<Length, 0, 4, HornerEvaluator>;
  using Degree5 = PoissonSeries<Length, 0, 5, HornerEvaluator>;
  Time const duration = 4.77553415434249021e-02 * Second;
  Instant const t_min = t0_;
  Instant const t_mid = t0_ + duration / 2;
  Instant const t_max = t0_ + duration;
  AngularFrequency const ω = 2.09400659210170170e+03 * Radian / Second;
  Degree4 const f(Degree4::AperiodicPolynomial({}, t0_),
                  {{ω,
                    {.sin = Degree4::PeriodicPolynomial(
                         {+5.10311065909077932e+00 * Metre,
                          +2.78062787709394854e+00 * Metre / Second,
                          +5.04290401496053242e+00 * Metre / Pow<2>(Second),
                          -7.27454632735125806e+00 * Metre / Pow<3>(Second),
                          +8.06537932856756612e+00 * Metre / Pow<4>(Second)},
                         t0_),
                     .cos = Degree4::PeriodicPolynomial(
                         {-8.11863376474325804e+00 * Metre,
                          +1.49140608216528037e+00 * Metre / Second,
                          -2.54224601087630298e+00 * Metre / Pow<2>(Second),
                          -4.52251796525658367e+00 * Metre / Pow<3>(Second),
                          -2.19458237171412751e+00 * Metre / Pow<4>(Second)},
                         t0_)}}});
  Degree5 const q(Degree5::AperiodicPolynomial({}, t0_),
                  {{ω,
                    {.sin = Degree5::PeriodicPolynomial(
                         {-4.41249783881549433e+01 * Metre,
                          +1.50208859053174347e+04 * Metre / Second,
                          -1.70674564621978020e+06 * Metre / Pow<2>(Second),
                          +8.52015772027946562e+07 * Metre / Pow<3>(Second),
                          -1.92799129073151779e+09 * Metre / Pow<4>(Second),
                          +1.61514557548221931e+10 * Metre / Pow<5>(Second)},
                         t0_),
                     .cos = Degree5::PeriodicPolynomial(
                         {-1.00752842659088765e-01 * Metre,
                          +2.25402995957193006e+01 * Metre / Second,
                          -1.66819064858902379e+03 * Metre / Pow<2>(Second),
                          +4.98682536071893774e+04 * Metre / Pow<3>(Second),
                          -5.18229522289936838e+05 * Metre / Pow<4>(Second),
                          0 * Metre / Pow<5>(Second)},
                         t0_)}}});

  // The integral is very small compared to the functions, so we end up in the
  // numerical noise, and adding more points would not help much.
  auto const product =
      InnerProduct(f.AtOrigin(t_mid),
                   q.AtOrigin(t_mid),
                   apodization::Hann<HornerEvaluator>(t_min, t_max),
                   t_min,
                   t_max);
  // Exact result obtained using Mathematica.
  EXPECT_THAT(product,
              RelativeErrorFrom(-4.848079980325297e-13 * Metre * Metre,
                                IsNear(0.19_(1))));
}

TEST_F(PoissonSeriesTest, PoorlyConditionedInnerProduct2) {
  using Degree0 = PoissonSeries<double, 0, 0, HornerEvaluator>;
  using Degree4 = PoissonSeries<double, 4, 4, HornerEvaluator>;
  Time const duration = +3.62955915932496390e-02 * Second;
  Instant const t_min = t0_;
  Instant const t_mid = t0_ + duration;
  Instant const t_max = t0_ + 2 * duration;
  AngularFrequency const ω1 = 2.63903139385469740e+03 * Radian / Second;
  AngularFrequency const ω2 = 2.75515553295453901e+03 * Radian / Second;
  AngularFrequency const ω3 = 2.75214520074802658e+03 * Radian / Second;
  Degree4 const f(
      Degree4::AperiodicPolynomial({-6.84494802606519098e-02,
                                    0 / Second,
                                    +7.62776153984628195e+02 / Pow<2>(Second),
                                    0 / Pow<3>(Second),
                                    -9.09000157727785874e+05 / Pow<4>(Second)},
                                   t_mid),
      {{ω1,
        {.sin = Degree4::PeriodicPolynomial(
             {0,
              -1.91877905970852938e+06 / Second,
              0 / Pow<2>(Second),
              +6.41616364182685137e+08 / Pow<3>(Second),
              0 / Pow<4>(Second)},
             t_mid),
         .cos = Degree4::PeriodicPolynomial(
             {+3.30073762193190414e+04,
              0 / Second,
              -4.82683293024840727e+07 / Pow<2>(Second),
              0 / Pow<3>(Second),
              +4.17242083788431835e+09 / Pow<4>(Second)},
             t_mid)}},
       {ω2,
        {.sin = Degree4::PeriodicPolynomial(
             {0,
              -1.91373931192037021e+06 / Second,
              0 / Pow<2>(Second),
              +6.34735730739231467e+08 / Pow<3>(Second),
              0 / Pow<4>(Second)},
             t_mid),
         .cos = Degree4::PeriodicPolynomial(
             {-3.30072258387235779e+04,
              0 / Second,
              +4.79716446803979352e+07 / Pow<2>(Second),
              0 / Pow<3>(Second),
              -4.08946453458770609e+09 / Pow<4>(Second)},
             t_mid)}}});

  Degree0 const g(Degree0::AperiodicPolynomial({}, t_mid),
                  {{ω3,
                    {.sin = Degree0::PeriodicPolynomial({}, t_mid),
                     .cos = Degree0::PeriodicPolynomial({1}, t_mid)}}});

  {
    auto const product = InnerProduct(f, g,
                     apodization::Dirichlet<HornerEvaluator>(t_min, t_max),
                     t_min, t_max);
    EXPECT_THAT(
        product,
        RelativeErrorFrom(
            +2.0267451184776034270e-11,
            AnyOf(IsNear(0.26_(1)), IsNear(0.33_(1)), IsNear(0.38_(1)))));
  }
  {
    auto const product = (PointwiseInnerProduct(f, g) *
                          apodization::Dirichlet<HornerEvaluator>(t_min, t_max))
                             .Integrate(t_min, t_max) /
                         (t_max - t_min);
    EXPECT_THAT(
        product,
        RelativeErrorFrom(+2.0267451184776034270e-11, IsNear(4.0e3_(1))));
  }
}

// This product occurs when orthogonalizing the basis for the Moon.
TEST_F(PoissonSeriesTest, PoorlyConditionedInnerProduct3) {
  using Series = PoissonSeries<double, 5, 3, EstrinEvaluator>;
  Time const duration = 3'945'600 * Second;
  Instant const t_min = t0_;
  Instant const t_mid = t0_ + duration;
  Instant const t_max = t0_ + 2 * duration;
  AngularFrequency const ω1 = 2.54633324931485220e-06 * Radian / Second;
  AngularFrequency const ω2 = 5.18914422934776064e-06 * Radian / Second;
  AngularFrequency const ω3 = 7.57219765654625247e-06 * Radian / Second;
  AngularFrequency const ω4 = 1.02018287773409324e-05 * Radian / Second;
  AngularFrequency const ω5 = 1.28742687976556365e-05 * Radian / Second;
  AngularFrequency const ω6 = 1.55555632099122024e-05 * Radian / Second;
  Series const f(
      Series::AperiodicPolynomial({0,
                                   +1.39108537630392846e+01 / Second,
                                   0 / Pow<2>(Second),
                                   -4.28766872062564245e-12 / Pow<3>(Second),
                                   0 / Pow<4>(Second),
                                   +2.19546428440080956e-25 / Pow<5>(Second)},
                                  t_mid),
      {{ω1,
        {.sin = Series::PeriodicPolynomial(
             {-8.35023222649048083e+06,
              0 / Second,
              +3.14056500770589297e-06 / Pow<2>(Second),
              0 / Pow<3>(Second)},
             t_mid),
         .cos = Series::PeriodicPolynomial(
             {0,
              +1.01168828383876246e+01 / Second,
              0 / Pow<2>(Second),
              -1.04312630151309924e-12 / Pow<3>(Second)},
             t_mid)}},
       {ω2,
        {.sin = Series::PeriodicPolynomial(
             {-4.91207898169697961e+05,
              0 / Second,
              +2.37420192215719412e-07 / Pow<2>(Second),
              0 / Pow<3>(Second)},
             t_mid),
         .cos = Series::PeriodicPolynomial(
             {0,
              +1.06674312857622944e+00 / Second,
              0 / Pow<2>(Second),
              -1.52459927383197132e-13 / Pow<3>(Second)},
             t_mid)}},
       {ω3,
        {.sin = Series::PeriodicPolynomial(
             {-1.85065842467988288e+05,
              0 / Second,
              +7.57044590566114320e-08 / Pow<2>(Second),
              0 / Pow<3>(Second)},
             t_mid),
         .cos = Series::PeriodicPolynomial(
             {0,
              +2.22354757856304569e-01 / Second,
              0 / Pow<2>(Second),
              -2.34730179758300617e-14 / Pow<3>(Second)},
             t_mid)}},
       {ω4,
        {.sin = Series::PeriodicPolynomial(
             {-1.09194717716096820e+04,
              0 / Second,
              +4.85897411923250085e-09 / Pow<2>(Second),
              0 / Pow<3>(Second)},
             t_mid),
         .cos = Series::PeriodicPolynomial(
             {0,
              +1.31475484344253241e-02 / Second,
              0 / Pow<2>(Second),
              -1.41217171514512962e-15 / Pow<3>(Second)},
             t_mid)}},
       {ω5,
        {.sin = Series::PeriodicPolynomial(
             {-4.86759403394113406e+02,
              0 / Second,
              +2.03337116963348534e-10 / Pow<2>(Second),
              0 / Pow<3>(Second)},
             t_mid),
         .cos = Series::PeriodicPolynomial(
             {0,
              +5.16201144683406409e-04 / Second,
              0 / Pow<2>(Second),
              -4.29187706226779285e-17 / Pow<3>(Second)},
             t_mid)}},
       {ω6,
        {.sin = Series::PeriodicPolynomial({-2.66335757891223146e+00,
                                             0 / Second,
                                             0 / Pow<2>(Second),
                                             0 / Pow<3>(Second)},
                                            t_mid),
         .cos = Series::PeriodicPolynomial({0,
                                             +3.60057656319038812e-06 / Second,
                                             0 / Pow<2>(Second),
                                             0 / Pow<3>(Second)},
                                            t_mid)}}});

  Series const g(
      Series::AperiodicPolynomial({0,
                                   +2.63050245340104294e-01 / Second,
                                   0 / Pow<2>(Second),
                                   -8.09823291316016917e-14 / Pow<3>(Second),
                                   0 / Pow<4>(Second),
                                   +4.14056618803521400e-27 / Pow<5>(Second)},
                                  t_mid),
      {{ω1,
        {.sin = Series::PeriodicPolynomial(
             {-1.57708056520423561e+05,
              0 / Second,
              +5.93414912295899838e-08 / Pow<2>(Second),
              0 / Pow<3>(Second)},
             t_mid),
         .cos = Series::PeriodicPolynomial(
             {0,
              +1.91902315754011399e-01 / Second,
              0 / Pow<2>(Second),
              -1.98598703231757471e-14 / Pow<3>(Second)},
             t_mid)}},
       {ω2,
        {.sin = Series::PeriodicPolynomial(
             {-9.31144597269088263e+03,
              0 / Second,
              +4.51121125348543404e-09 / Pow<2>(Second),
              0 / Pow<3>(Second)},
             t_mid),
         .cos = Series::PeriodicPolynomial(
             {0,
              +2.08585690999323684e-02 / Second,
              0 / Pow<2>(Second),
              -2.99826421429525608e-15 / Pow<3>(Second)},
             t_mid)}},
       {ω3,
        {.sin = Series::PeriodicPolynomial(
             {-3.71934064860903527e+03,
              0 / Second,
              +1.52991642654451839e-09 / Pow<2>(Second),
              0 / Pow<3>(Second)},
             t_mid),
         .cos = Series::PeriodicPolynomial(
             {0,
              +4.54325428407370037e-03 / Second,
              0 / Pow<2>(Second),
              -4.88811294904941287e-16 / Pow<3>(Second)},
             t_mid)}},
       {ω4,
        {.sin = Series::PeriodicPolynomial(
             {-2.38782465236558210e+02,
              0 / Second,
              +1.07675550289214021e-10 / Pow<2>(Second),
              0 / Pow<3>(Second)},
             t_mid),
         .cos = Series::PeriodicPolynomial(
             {0,
              +2.96344251787576066e-04 / Second,
              0 / Pow<2>(Second),
              -3.31546592884930502e-17 / Pow<3>(Second)},
             t_mid)}},
       {ω5,
        {.sin = Series::PeriodicPolynomial(
             {-1.29376251439425527e+01,
              0 / Second,
              +5.63467102431514781e-12 / Pow<2>(Second),
              0 / Pow<3>(Second)},
             t_mid),
         .cos = Series::PeriodicPolynomial(
             {0,
              +1.43747865681285978e-05 / Second,
              0 / Pow<2>(Second),
              -1.34490730453117390e-18 / Pow<3>(Second)},
             t_mid)}},
       {ω6,
        {.sin = Series::PeriodicPolynomial(
             {-2.29407203043806407e-01,
              0 / Second,
              +6.42353197319720082e-14 / Pow<2>(Second),
              0 / Pow<3>(Second)},
             t_mid),
         .cos = Series::PeriodicPolynomial(
             {0,
              +1.76202789977736386e-07 / Second,
              0 / Pow<2>(Second),
              0 / Pow<3>(Second)},
             t_mid)}}});

  double const expected_product = -9.93743685469196454745e-9;
  {
    auto const product =
        InnerProduct(f, g,
                     apodization::Dirichlet<EstrinEvaluator>(t_min, t_max),
                     t_min, t_max);
    EXPECT_THAT(product,
                RelativeErrorFrom(expected_product,
                                  AnyOf(IsNear(0.00069_(1)),
                                        IsNear(0.0013_(1)),
                                        IsNear(0.0015_(1)))));
  }
  // This test demonstrates how bad Integrate can be, for products that arise in
  // practice.  Exact integration of the result of PointwiseInnerProduct yields
  // the correct value to 7 digits, but Integrate suffers from enormous
  // cancellations: the intermediate terms may be as large as 6.2e12,
  // effectively losing 69 bits.  In other words, there is no hope of computing
  // this product using double.
  {
    auto const product = (PointwiseInnerProduct(f, g) *
                          apodization::Dirichlet<EstrinEvaluator>(t_min, t_max))
                             .Integrate(t_min, t_max) /
                         (t_max - t_min);
    EXPECT_THAT(
        product,
        RelativeErrorFrom(
            expected_product,
            AnyOf(IsNear(7.7e6_(1)), IsNear(8.8e6_(1)), IsNear(1.0e7_(1)))));
  }
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

}  // namespace numerics
}  // namespace principia
