
#include "numerics/newhall.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace numerics {

using geometry::Instant;
using quantities::Abs;
using quantities::Length;
using quantities::Speed;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::ApproximateQuantity;
using testing_utilities::IsNear;
using testing_utilities::RelativeError;
using testing_utilities::operator""_⑴;

// The adapters wrap the result of the Newhall approximation so that they can be
// used consistently in this test.

template<int degree>
class ЧебышёвAdapter {
 public:
  static ЧебышёвAdapter NewhallApproximation(std::vector<Length> const& q,
                                             std::vector<Speed> const& v,
                                             Instant const& t_min,
                                             Instant const& t_max,
                                             Length& error_estimate);

  Length Evaluate(Instant const& t) const;
  Speed EvaluateDerivative(Instant const& t) const;

 private:
  explicit ЧебышёвAdapter(ЧебышёвSeries<Length> series);

  ЧебышёвSeries<Length> series_;
};

template<int degree>
class MonomialAdapter {
 public:
  static MonomialAdapter NewhallApproximation(std::vector<Length> const& q,
                                              std::vector<Speed> const& v,
                                              Instant const& t_min,
                                              Instant const& t_max,
                                              Length& error_estimate);

  Length Evaluate(Instant const& t) const;
  Speed EvaluateDerivative(Instant const& t) const;

 private:
  using P = PolynomialInMonomialBasis<Length, Instant, degree, EstrinEvaluator>;
  explicit MonomialAdapter(P const& polynomial);

  P polynomial_;
};

class NewhallTest : public ::testing::Test {
 protected:
  NewhallTest()
      : t_min_(t0_ - 1 * Second),
        t_max_(t0_ + 3 * Second),

        length_function_1_([this](Instant const t) -> Length {
          return 0.5 * Metre + 2 * Metre *
                                   std::sin((t - t_min_) / (0.3 * Second)) *
                                   std::exp((t - t_min_) / (1 * Second));
        }),
        speed_function_1_([this](Instant const t) -> Speed {
          return ((2 * Metre) / (0.3 * Second) *
                      std::cos((t - t_min_) / (0.3 * Second)) +
                  (2 * Metre / Second) *
                      std::sin((t - t_min_) / (0.3 * Second))) *
                 std::exp((t - t_min_) / (1 * Second));
        }),

        length_function_2_([this](Instant const t) -> Length {
          return 5 * Metre *
                 (1 + (t - t_min_) / (0.3 * Second) +
                  std::pow((t - t_min_) / (4 * Second), 7));
        }),
        speed_function_2_([this](Instant const t) -> Speed {
          return 5 * Metre *
                 (1 / (0.3 * Second) +
                  (7 / (4 * Second)) *
                      std::pow((t - t_min_) / (4 * Second), 6));
        }) {}

  template<typename Adapter>
  void CheckNewhallApproximationErrors(
      std::function<Length(Instant const)> length_function,
      std::function<Speed(Instant const)> speed_function,
      ApproximateQuantity<Length> const& expected_length_error_estimate,
      ApproximateQuantity<Length> const& expected_length_absolute_error,
      ApproximateQuantity<Speed> const& expected_speed_absolute_error) {
    std::vector<Length> lengths;
    std::vector<Speed> speeds;
    for (Instant t = t_min_; t <= t_max_; t += 0.5 * Second) {
      lengths.push_back(length_function(t));
      speeds.push_back(speed_function(t));
    }

    Length length_error_estimate;
    Adapter const approximation =
        Adapter::NewhallApproximation(
            lengths, speeds, t_min_, t_max_, length_error_estimate);

    // Compute the absolute error of both functions throughout the interval.
    Length length_absolute_error;
    Speed speed_absolute_error;
    for (Instant t = t_min_; t <= t_max_; t += 0.05 * Second) {
      Length const expected_length = length_function(t);
      Length const actual_length = approximation.Evaluate(t);
      Speed const expected_speed = speed_function(t);
      Speed const actual_speed = approximation.EvaluateDerivative(t);
      length_absolute_error =
          std::max(length_absolute_error,
                    AbsoluteError(expected_length, actual_length));
      speed_absolute_error =
          std::max(speed_absolute_error,
                    AbsoluteError(expected_speed, actual_speed));
    }

    EXPECT_THAT(Abs(length_error_estimate),
                IsNear(expected_length_error_estimate));
    EXPECT_THAT(length_absolute_error, IsNear(expected_length_absolute_error));
    EXPECT_THAT(speed_absolute_error, IsNear(expected_speed_absolute_error));
  }

  Instant const t0_;
  Instant t_min_;
  Instant t_max_;
  std::function<Length(Instant const& t)> length_function_1_;
  std::function<Speed(Instant const& t)> speed_function_1_;
  std::function<Length(Instant const& t)> length_function_2_;
  std::function<Speed(Instant const& t)> speed_function_2_;
};

template<int degree>
ЧебышёвAdapter<degree> ЧебышёвAdapter<degree>::NewhallApproximation(
    std::vector<Length> const& q,
    std::vector<Speed> const& v,
    Instant const& t_min,
    Instant const& t_max,
    Length& error_estimate) {
  return ЧебышёвAdapter(NewhallApproximationInЧебышёвBasis(degree,
                                                           q, v,
                                                           t_min, t_max,
                                                           error_estimate));
}

template<int degree>
Length ЧебышёвAdapter<degree>::Evaluate(Instant const& t) const {
  return series_.Evaluate(t);
}

template<int degree>
Speed ЧебышёвAdapter<degree>::EvaluateDerivative(Instant const& t) const {
  return series_.EvaluateDerivative(t);
}

template<int degree>
ЧебышёвAdapter<degree>::ЧебышёвAdapter(ЧебышёвSeries<Length> series)
    : series_(std::move(series)) {}

template<int degree>
MonomialAdapter<degree> MonomialAdapter<degree>::NewhallApproximation(
    std::vector<Length> const& q,
    std::vector<Speed> const& v,
    Instant const& t_min,
    Instant const& t_max,
    Length& error_estimate) {
  return MonomialAdapter(
      NewhallApproximationInMonomialBasis<Length, degree, EstrinEvaluator>(
          q, v,
          t_min, t_max,
          error_estimate));
}

template<int degree>
Length MonomialAdapter<degree>::Evaluate(Instant const& t) const {
  return polynomial_.Evaluate(t);
}

template<int degree>
Speed MonomialAdapter<degree>::EvaluateDerivative(Instant const& t) const {
  return polynomial_.EvaluateDerivative(t);
}

template<int degree>
MonomialAdapter<degree>::MonomialAdapter(P const& polynomial)
    : polynomial_(polynomial) {}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_3) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<3>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/3.9e1_⑴ * Metre,
      /*expected_length_absolute_error=*/1.7e2_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.3e2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_4) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<4>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/1.8e1_⑴ * Metre,
      /*expected_length_absolute_error=*/4.7e1_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.2e2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_5) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<5>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/8.3e0_⑴ * Metre,
      /*expected_length_absolute_error=*/4.3e1_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.2e2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_6) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<6>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/1.4e0_⑴ * Metre,
      /*expected_length_absolute_error=*/3.8e1_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.0e2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_7) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<7>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/7.5e0_⑴ * Metre,
      /*expected_length_absolute_error=*/1.4e1_⑴ * Metre,
      /*expected_speed_absolute_error=*/4.5e1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_8) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<8>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/4.2e0_⑴ * Metre,
      /*expected_length_absolute_error=*/6.3e0_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.8e1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_9) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<9>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/3.1e-1_⑴ * Metre,
      /*expected_length_absolute_error=*/4.9e0_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.2e1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_10) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<10>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/7.5e-1_⑴ * Metre,
      /*expected_length_absolute_error=*/6.5e-1_⑴ * Metre,
      /*expected_speed_absolute_error=*/3.6e0_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_11) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<11>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/1.8e-1_⑴ * Metre,
      /*expected_length_absolute_error=*/2.0e-1_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.6e0_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_12) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<12>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/4.9e-2_⑴ * Metre,
      /*expected_length_absolute_error=*/7.9e-2_⑴ * Metre,
      /*expected_speed_absolute_error=*/7.3e-1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_13) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<13>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/2.6e-2_⑴ * Metre,
      /*expected_length_absolute_error=*/1.3e-2_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.2e-1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_14) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<14>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/8.9e-4_⑴ * Metre,
      /*expected_length_absolute_error=*/1.5e-2_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.5e-1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_15) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<15>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/1.6e-3_⑴ * Metre,
      /*expected_length_absolute_error=*/4.3e-3_⑴ * Metre,
      /*expected_speed_absolute_error=*/4.4e-2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_16) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<16>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/2.9e-4_⑴ * Metre,
      /*expected_length_absolute_error=*/1.7e-3_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.8e-2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_17) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<17>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/3.7e-5_⑴ * Metre,
      /*expected_length_absolute_error=*/7.6e-4_⑴ * Metre,
      /*expected_speed_absolute_error=*/8.2e-3_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_3) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<3>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/7.9e-1_⑴ * Metre,
      /*expected_length_absolute_error=*/2.0e0_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.8e0_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_4) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<4>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/2.5e-1_⑴ * Metre,
      /*expected_length_absolute_error=*/2.9e-1_⑴ * Metre,
      /*expected_speed_absolute_error=*/4.6e-1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_5) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<5>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/5.7e-2_⑴ * Metre,
      /*expected_length_absolute_error=*/3.6e-2_⑴ * Metre,
      /*expected_speed_absolute_error=*/7.4e-2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_6) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<6>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/8.6e-3_⑴ * Metre,
      /*expected_length_absolute_error=*/2.3e-3_⑴ * Metre,
      /*expected_speed_absolute_error=*/6.0e-3_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_7) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<7>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/6.2e-4_⑴ * Metre,
      /*expected_length_absolute_error=*/2.9e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.5e-14_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_8) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<8>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/4.5e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/2.9e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.2e-14_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_9) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<9>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/4.8e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/2.9e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.2e-14_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_10) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<10>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/4.8e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/4.3e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.9e-14_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_11) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<11>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/1.1e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/3.6e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.2e-14_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_12) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<12>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/7.3e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/2.9e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.9e-14_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_13) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<13>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/1.3e-15_⑴ * Metre,
      /*expected_length_absolute_error=*/1.4e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/6.1e-14_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_14) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<14>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/2.1e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/1.4e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/6.8e-14_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_15) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<15>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/3.2e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/2.9e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/3.5e-13_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_16) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<16>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/1.1e-15_⑴ * Metre,
      /*expected_length_absolute_error=*/2.9e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/8.5e-13_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_17) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter<17>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/4.8e-15_⑴ * Metre,
      /*expected_length_absolute_error=*/7.2e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.2e-12_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_3) {
  CheckNewhallApproximationErrors<MonomialAdapter<3>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/3.9e1_⑴ * Metre,
      /*expected_length_absolute_error=*/1.7e2_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.3e2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_4) {
  CheckNewhallApproximationErrors<MonomialAdapter<4>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/1.8e1_⑴ * Metre,
      /*expected_length_absolute_error=*/4.7e1_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.2e2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_5) {
  CheckNewhallApproximationErrors<MonomialAdapter<5>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/8.3e0_⑴ * Metre,
      /*expected_length_absolute_error=*/4.3e1_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.2e2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_6) {
  CheckNewhallApproximationErrors<MonomialAdapter<6>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/1.4e0_⑴ * Metre,
      /*expected_length_absolute_error=*/3.8e1_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.0e2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_7) {
  CheckNewhallApproximationErrors<MonomialAdapter<7>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/7.5e0_⑴ * Metre,
      /*expected_length_absolute_error=*/1.4e1_⑴ * Metre,
      /*expected_speed_absolute_error=*/4.5e1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_8) {
  CheckNewhallApproximationErrors<MonomialAdapter<8>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/4.2e0_⑴ * Metre,
      /*expected_length_absolute_error=*/6.3e0_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.8e1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_9) {
  CheckNewhallApproximationErrors<MonomialAdapter<9>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/3.1e-1_⑴ * Metre,
      /*expected_length_absolute_error=*/4.9e0_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.2e1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_10) {
  CheckNewhallApproximationErrors<MonomialAdapter<10>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/7.5e-1_⑴ * Metre,
      /*expected_length_absolute_error=*/6.5e-1_⑴ * Metre,
      /*expected_speed_absolute_error=*/3.6e0_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_11) {
  CheckNewhallApproximationErrors<MonomialAdapter<11>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/1.8e-1_⑴ * Metre,
      /*expected_length_absolute_error=*/2.0e-1_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.6e0_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_12) {
  CheckNewhallApproximationErrors<MonomialAdapter<12>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/4.9e-2_⑴ * Metre,
      /*expected_length_absolute_error=*/7.9e-2_⑴ * Metre,
      /*expected_speed_absolute_error=*/7.3e-1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_13) {
  CheckNewhallApproximationErrors<MonomialAdapter<13>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/2.6e-2_⑴ * Metre,
      /*expected_length_absolute_error=*/1.3e-2_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.2e-1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_14) {
  CheckNewhallApproximationErrors<MonomialAdapter<14>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/8.9e-4_⑴ * Metre,
      /*expected_length_absolute_error=*/1.5e-2_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.5e-1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_15) {
  CheckNewhallApproximationErrors<MonomialAdapter<15>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/1.6e-3_⑴ * Metre,
      /*expected_length_absolute_error=*/4.3e-3_⑴ * Metre,
      /*expected_speed_absolute_error=*/4.4e-2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_16) {
  CheckNewhallApproximationErrors<MonomialAdapter<16>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/2.9e-4_⑴ * Metre,
      /*expected_length_absolute_error=*/1.7e-3_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.8e-2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_17) {
  CheckNewhallApproximationErrors<MonomialAdapter<17>>(
      length_function_1_,
      speed_function_1_,
      /*expected_length_error_estimate=*/3.7e-5_⑴ * Metre,
      /*expected_length_absolute_error=*/7.6e-4_⑴ * Metre,
      /*expected_speed_absolute_error=*/8.2e-3_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_3) {
  CheckNewhallApproximationErrors<MonomialAdapter<3>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/7.9e-1_⑴ * Metre,
      /*expected_length_absolute_error=*/2.0e0_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.8e0_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_4) {
  CheckNewhallApproximationErrors<MonomialAdapter<4>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/2.5e-1_⑴ * Metre,
      /*expected_length_absolute_error=*/2.9e-1_⑴ * Metre,
      /*expected_speed_absolute_error=*/4.6e-1_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_5) {
  CheckNewhallApproximationErrors<MonomialAdapter<5>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/5.7e-2_⑴ * Metre,
      /*expected_length_absolute_error=*/3.6e-2_⑴ * Metre,
      /*expected_speed_absolute_error=*/7.4e-2_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_6) {
  CheckNewhallApproximationErrors<MonomialAdapter<6>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/8.6e-3_⑴ * Metre,
      /*expected_length_absolute_error=*/2.3e-3_⑴ * Metre,
      /*expected_speed_absolute_error=*/6.0e-3_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_7) {
  CheckNewhallApproximationErrors<MonomialAdapter<7>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/6.2e-4_⑴ * Metre,
      /*expected_length_absolute_error=*/2.9e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.5e-14_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_8) {
  CheckNewhallApproximationErrors<MonomialAdapter<8>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/4.5e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/3.9e-14_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.2e-13_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_9) {
  CheckNewhallApproximationErrors<MonomialAdapter<9>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/4.8e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/1.0e-13_⑴ * Metre,
      /*expected_speed_absolute_error=*/4.8e-13_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_10) {
  CheckNewhallApproximationErrors<MonomialAdapter<10>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/4.8e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/3.5e-13_⑴ * Metre,
      /*expected_speed_absolute_error=*/7.1e-13_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_11) {
  CheckNewhallApproximationErrors<MonomialAdapter<11>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/1.1e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/5.6e-13_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.7e-12_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_12) {
  CheckNewhallApproximationErrors<MonomialAdapter<12>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/7.3e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/1.3e-11_⑴ * Metre,
      /*expected_speed_absolute_error=*/5.4e-11_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_13) {
  CheckNewhallApproximationErrors<MonomialAdapter<13>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/1.3e-15_⑴ * Metre,
      /*expected_length_absolute_error=*/1.9e-11_⑴ * Metre,
      /*expected_speed_absolute_error=*/9.5e-11_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_14) {
  CheckNewhallApproximationErrors<MonomialAdapter<14>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/2.1e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/8.1e-12_⑴ * Metre,
      /*expected_speed_absolute_error=*/5.1e-11_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_15) {
  CheckNewhallApproximationErrors<MonomialAdapter<15>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/3.2e-16_⑴ * Metre,
      /*expected_length_absolute_error=*/3.0e-10_⑴ * Metre,
      /*expected_speed_absolute_error=*/1.6e-9_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_16) {
  CheckNewhallApproximationErrors<MonomialAdapter<16>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/1.1e-15_⑴ * Metre,
      /*expected_length_absolute_error=*/5.8e-10_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.8e-9_⑴ * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_17) {
  CheckNewhallApproximationErrors<MonomialAdapter<17>>(
      length_function_2_,
      speed_function_2_,
      /*expected_length_error_estimate=*/4.8e-15_⑴ * Metre,
      /*expected_length_absolute_error=*/3.7e-9_⑴ * Metre,
      /*expected_speed_absolute_error=*/2.2e-8_⑴ * Metre / Second);
}

TEST_F(NewhallTest, NonConstantDegree) {
    std::vector<Length> lengths;
    std::vector<Speed> speeds;
    for (Instant t = t_min_; t <= t_max_; t += 0.5 * Second) {
      lengths.push_back(length_function_1_(t));
      speeds.push_back(speed_function_1_(t));
    }

    Length length_error_estimate;
    auto const approximation =
        NewhallApproximationInMonomialBasis<Length, EstrinEvaluator>(
            /*degree=*/10,
            lengths, speeds, t_min_, t_max_, length_error_estimate);

    EXPECT_THAT(RelativeError(approximation->Evaluate(t_min_),
                              length_function_1_(t_min_)), IsNear(9e-13_⑴));
}

}  // namespace numerics
}  // namespace principia
