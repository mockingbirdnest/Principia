
#include "numerics/fit_hermite_spline.hpp"

#include <list>
#include <vector>

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/double_precision.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {

using base::Range;
using geometry::Instant;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Sin;
using quantities::Length;
using quantities::Speed;
using quantities::si::Centi;
using quantities::si::Metre;
using quantities::si::Micro;
using quantities::si::Milli;
using quantities::si::Nano;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::IsNear;
using testing_utilities::operator""_⑴;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::ResultOf;

namespace numerics {

class FitHermiteSplineTest : public ::testing::Test {
 protected:
  struct Sample {
    Instant t;
    Length x;
    Speed v;
  };
  Instant const t0_;
};

using FitHermiteSplineDeathTest = FitHermiteSplineTest;

TEST_F(FitHermiteSplineTest, Sinusoid) {
  AngularFrequency const ω = 1 * Radian / Second;
  auto const f = [ω, this](Instant const& t) {
    return Cos(ω * (t - t0_)) * Metre;
  };
  auto const df = [ω, this](Instant const& t) {
    return -ω * Sin(ω *(t - t0_)) * Metre / Radian;
  };
  std::vector<Sample> samples;
  {
    auto t = DoublePrecision<Instant>(t0_);
    for (; t.value < t0_ + π / 2 * Second; t.Increment(100 * Milli(Second))) {
      samples.push_back({t.value, f(t.value), df(t.value)});
    }
    for (; t.value < t0_ + π * Second; t.Increment(20 * Milli(Second))) {
      samples.push_back({t.value, f(t.value), df(t.value)});
    }
  }
  std::list<std::vector<Sample>::const_iterator> const interpolation_points =
      FitHermiteSpline<Instant, Length>(
          samples,
          [](auto&& sample) -> auto&& { return sample.t; },
          [](auto&& sample) -> auto&& { return sample.x; },
          [](auto&& sample) -> auto&& { return sample.v; },
          1 * Centi(Metre));

  // Note that gmock doesn't do decltypes, so we can't pass a λ directly.
  // Also note that |Pointee| doesn't work with iterators, so
  // |Pointee(Field(&Sample::t, _))| is not an option.
  std::function<Instant(std::vector<Sample>::const_iterator)> const get_time =
      [](auto const it) { return it->t; };
  EXPECT_THAT(interpolation_points,
              ElementsAre(ResultOf(get_time, Eq(t0_ + 1.5 * Second)),
                          ResultOf(get_time, Eq(t0_ + 3.06 * Second))));

  auto lower_bound = samples.cbegin();
  auto upper_bound = interpolation_points.front();
  Hermite3<Instant, Length> first_polynomial({lower_bound->t, upper_bound->t},
                                             {lower_bound->x, upper_bound->x},
                                             {lower_bound->v, upper_bound->v});
  EXPECT_THAT(first_polynomial.LInfinityError(
                  Range(lower_bound, upper_bound),
                  [](auto&& sample) -> auto&& { return sample.t; },
                  [](auto&& sample) -> auto&& { return sample.x; }),
              IsNear(9.3_⑴ * Milli(Metre)));
  lower_bound = upper_bound;
  upper_bound = interpolation_points.back();
  Hermite3<Instant, Length> second_polynomial({lower_bound->t, upper_bound->t},
                                              {lower_bound->x, upper_bound->x},
                                              {lower_bound->v, upper_bound->v});
  EXPECT_THAT(second_polynomial.LInfinityError(
                  Range(lower_bound, upper_bound),
                  [](auto&& sample) -> auto&& { return sample.t; },
                  [](auto&& sample) -> auto&& { return sample.x; }),
              IsNear(1.0_⑴ * Centi(Metre)));
  lower_bound = upper_bound;
  upper_bound = samples.cend() - 1;
  Hermite3<Instant, Length> tail_polynomial({lower_bound->t, upper_bound->t},
                                            {lower_bound->x, upper_bound->x},
                                            {lower_bound->v, upper_bound->v});
  EXPECT_THAT(tail_polynomial.LInfinityError(
                  Range(lower_bound, upper_bound),
                  [](auto&& sample) -> auto&& { return sample.t; },
                  [](auto&& sample) -> auto&& { return sample.x; }),
              IsNear(107_⑴ * Nano(Metre)));
}

TEST_F(FitHermiteSplineDeathTest, NoDownsampling) {
  AngularFrequency const ω = 1 * Radian / Second;
  auto const f = [ω, this](Instant const& t) {
    return Cos(ω * (t - t0_)) * Metre;
  };
  auto const df = [ω, this](Instant const& t) {
    return -ω * Sin(ω *(t - t0_)) * Metre / Radian;
  };
  std::vector<Sample> samples;
  {
    auto t = DoublePrecision<Instant>(t0_);
    for (; t.value < t0_ + π / 2 * Second; t.Increment(100 * Milli(Second))) {
      samples.push_back({t.value, f(t.value), df(t.value)});
    }
    for (; t.value < t0_ + π * Second; t.Increment(20 * Milli(Second))) {
      samples.push_back({t.value, f(t.value), df(t.value)});
    }
  }
  auto fit_hermite_spline = [&samples]() {
    return FitHermiteSpline<Instant, Length>(
        samples,
        [](auto&& sample) -> auto&& { return sample.t; },
        [](auto&& sample) -> auto&& { return sample.x; },
        [](auto&& sample) -> auto&& { return sample.v; },
        0 * Metre);
  };

  EXPECT_DEATH({
    std::list<std::vector<Sample>::const_iterator> const
        interpolation_points = fit_hermite_spline();
  }, "tail.size.*samples.size");
}

}  // namespace numerics
}  // namespace principia
