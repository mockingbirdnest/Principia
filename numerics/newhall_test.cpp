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

using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::ApproximateQuantity;
using testing_utilities::IsNear;
using testing_utilities::RelativeError;
using namespace principia::geometry::_named_quantities;
using namespace principia::numerics::_newhall;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using testing_utilities::operator""_;

// The adapters wrap the result of the Newhall approximation so that they can be
// used consistently in this test.

template<typename Value, int degree>
class ЧебышёвAdapter {
 public:
  static ЧебышёвAdapter NewhallApproximation(
      std::vector<Value> const& q,
      std::vector<Variation<Value>> const& v,
      Instant const& t_min,
      Instant const& t_max,
      Difference<Value>& error_estimate);

  Value Evaluate(Instant const& t) const;
  Variation<Value> EvaluateDerivative(Instant const& t) const;

 private:
  explicit ЧебышёвAdapter(ЧебышёвSeries<Value> series);

  ЧебышёвSeries<Value> series_;
};

template<typename Value, int degree>
class MonomialAdapter {
 public:
  static MonomialAdapter NewhallApproximation(
      std::vector<Value> const& q,
      std::vector<Variation<Value>> const& v,
      Instant const& t_min,
      Instant const& t_max,
      Difference<Value>& error_estimate);

  Value Evaluate(Instant const& t) const;
  Variation<Value> EvaluateDerivative(Instant const& t) const;

 private:
  using P =
      PolynomialInMonomialBasis<Value, Instant, degree, EstrinEvaluator>;
  explicit MonomialAdapter(P const& polynomial);

  P polynomial_;
};

template<typename Value, int degree>
ЧебышёвAdapter<Value, degree>
ЧебышёвAdapter<Value, degree>::NewhallApproximation(
    std::vector<Value> const& q,
    std::vector<Variation<Value>> const& v,
    Instant const& t_min,
    Instant const& t_max,
    Difference<Value>& error_estimate) {
  return ЧебышёвAdapter(NewhallApproximationInЧебышёвBasis(degree,
                                                           q, v,
                                                           t_min, t_max,
                                                           error_estimate));
}

template<typename Value, int degree>
Value ЧебышёвAdapter<Value, degree>::Evaluate(Instant const& t) const {
  return series_.Evaluate(t);
}

template<typename Value, int degree>
Variation<Value>
ЧебышёвAdapter<Value, degree>::EvaluateDerivative(Instant const& t) const {
  return series_.EvaluateDerivative(t);
}

template<typename Value, int degree>
ЧебышёвAdapter<Value, degree>::ЧебышёвAdapter(ЧебышёвSeries<Value> series)
    : series_(std::move(series)) {}

template<typename Value, int degree>
MonomialAdapter<Value, degree>
MonomialAdapter<Value, degree>::NewhallApproximation(
    std::vector<Value> const& q,
    std::vector<Variation<Value>> const& v,
    Instant const& t_min,
    Instant const& t_max,
    Difference<Value>& error_estimate) {
  return MonomialAdapter(
      NewhallApproximationInMonomialBasis<Value, degree, EstrinEvaluator>(
          q, v,
          t_min, t_max,
          error_estimate));
}

template<typename Value, int degree>
Value MonomialAdapter<Value, degree>::Evaluate(Instant const& t) const {
  return polynomial_(t);
}

template<typename Value, int degree>
Variation<Value>
MonomialAdapter<Value, degree>::EvaluateDerivative(Instant const& t) const {
  return polynomial_.EvaluateDerivative(t);
}

template<typename Value, int degree>
MonomialAdapter<Value, degree>::MonomialAdapter(P const& polynomial)
    : polynomial_(polynomial) {}

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

  template<template<typename, int> typename Adapter,
           typename Value, int degree>
  void CheckNewhallApproximationErrors(
      std::function<Value(Instant const)> length_function,
      std::function<Variation<Value>(Instant const)> speed_function,
      ApproximateQuantity<Difference<Value>> const&
          expected_value_error_estimate,
      ApproximateQuantity<Difference<Value>> const&
          expected_value_absolute_error,
      ApproximateQuantity<Variation<Value>> const&
          expected_variation_absolute_error,
      std::optional<ApproximateQuantity<Difference<Value>>> const&
          expected_value_absolute_error_with_fma = std::nullopt) {
    std::vector<Value> arguments;
    std::vector<Variation<Value>> variations;
    for (Instant t = t_min_; t <= t_max_; t += 0.5 * Second) {
      arguments.push_back(length_function(t));
      variations.push_back(speed_function(t));
    }

    Difference<Value> argument_error_estimate;
    auto const approximation =
        Adapter<Value, degree>::NewhallApproximation(
            arguments, variations, t_min_, t_max_, argument_error_estimate);

    // Compute the absolute error of both functions throughout the interval.
    Difference<Value> argument_absolute_error;
    Variation<Value> variation_absolute_error;
    for (Instant t = t_min_; t <= t_max_; t += 0.05 * Second) {
      Value const expected_value = length_function(t);
      Value const actual_value = approximation.Evaluate(t);
      Variation<Value> const expected_variation = speed_function(t);
      Variation<Value> const actual_variation =
          approximation.EvaluateDerivative(t);
      argument_absolute_error =
          std::max(argument_absolute_error,
                   AbsoluteError(expected_value, actual_value));
      variation_absolute_error =
          std::max(variation_absolute_error,
                   AbsoluteError(expected_variation, actual_variation));
    }

    EXPECT_THAT(Abs(argument_error_estimate),
                IsNear(expected_value_error_estimate));
    if (expected_value_absolute_error_with_fma.has_value() && UseHardwareFMA) {
      EXPECT_THAT(argument_absolute_error,
                  IsNear(*expected_value_absolute_error_with_fma));
    } else {
      EXPECT_THAT(argument_absolute_error,
                  IsNear(expected_value_absolute_error));
    }
    EXPECT_THAT(variation_absolute_error,
                IsNear(expected_variation_absolute_error));
  }

  Instant const t0_;
  Instant t_min_;
  Instant t_max_;
  std::function<Length(Instant const& t)> length_function_1_;
  std::function<Speed(Instant const& t)> speed_function_1_;
  std::function<Length(Instant const& t)> length_function_2_;
  std::function<Speed(Instant const& t)> speed_function_2_;
};

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_3) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 3>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/3.9e1_(1) * Metre,
      /*expected_value_absolute_error=*/1.7e2_(1) * Metre,
      /*expected_variation_absolute_error=*/2.3e2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_4) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 4>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/1.8e1_(1) * Metre,
      /*expected_value_absolute_error=*/4.7e1_(1) * Metre,
      /*expected_variation_absolute_error=*/1.2e2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_5) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 5>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/8.3e0_(1) * Metre,
      /*expected_value_absolute_error=*/4.3e1_(1) * Metre,
      /*expected_variation_absolute_error=*/1.2e2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_6) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 6>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/1.4e0_(1) * Metre,
      /*expected_value_absolute_error=*/3.8e1_(1) * Metre,
      /*expected_variation_absolute_error=*/1.0e2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_7) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 7>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/7.5e0_(1) * Metre,
      /*expected_value_absolute_error=*/1.4e1_(1) * Metre,
      /*expected_variation_absolute_error=*/4.5e1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_8) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 8>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/4.2e0_(1) * Metre,
      /*expected_value_absolute_error=*/6.3e0_(1) * Metre,
      /*expected_variation_absolute_error=*/2.8e1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_9) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 9>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/3.1e-1_(1) * Metre,
      /*expected_value_absolute_error=*/4.9e0_(1) * Metre,
      /*expected_variation_absolute_error=*/2.2e1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_10) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 10>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/7.5e-1_(1) * Metre,
      /*expected_value_absolute_error=*/6.5e-1_(1) * Metre,
      /*expected_variation_absolute_error=*/3.6e0_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_11) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 11>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/1.8e-1_(1) * Metre,
      /*expected_value_absolute_error=*/2.0e-1_(1) * Metre,
      /*expected_variation_absolute_error=*/1.6e0_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_12) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 12>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/4.9e-2_(1) * Metre,
      /*expected_value_absolute_error=*/7.9e-2_(1) * Metre,
      /*expected_variation_absolute_error=*/7.3e-1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_13) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 13>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/2.6e-2_(1) * Metre,
      /*expected_value_absolute_error=*/1.3e-2_(1) * Metre,
      /*expected_variation_absolute_error=*/1.2e-1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_14) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 14>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/8.9e-4_(1) * Metre,
      /*expected_value_absolute_error=*/1.5e-2_(1) * Metre,
      /*expected_variation_absolute_error=*/1.5e-1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_15) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 15>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/1.6e-3_(1) * Metre,
      /*expected_value_absolute_error=*/4.3e-3_(1) * Metre,
      /*expected_variation_absolute_error=*/4.4e-2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_16) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 16>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/2.9e-4_(1) * Metre,
      /*expected_value_absolute_error=*/1.7e-3_(1) * Metre,
      /*expected_variation_absolute_error=*/1.8e-2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_1_17) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 17>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/3.7e-5_(1) * Metre,
      /*expected_value_absolute_error=*/7.6e-4_(1) * Metre,
      /*expected_variation_absolute_error=*/8.2e-3_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_3) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 3>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/7.9e-1_(1) * Metre,
      /*expected_value_absolute_error=*/2.0e0_(1) * Metre,
      /*expected_variation_absolute_error=*/1.8e0_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_4) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 4>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/2.5e-1_(1) * Metre,
      /*expected_value_absolute_error=*/2.9e-1_(1) * Metre,
      /*expected_variation_absolute_error=*/4.6e-1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_5) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 5>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/5.7e-2_(1) * Metre,
      /*expected_value_absolute_error=*/3.6e-2_(1) * Metre,
      /*expected_variation_absolute_error=*/7.4e-2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_6) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 6>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/8.6e-3_(1) * Metre,
      /*expected_value_absolute_error=*/2.3e-3_(1) * Metre,
      /*expected_variation_absolute_error=*/6.0e-3_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_7) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 7>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/6.2e-4_(1) * Metre,
      /*expected_value_absolute_error=*/2.9e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/2.5e-14_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_8) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 8>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/4.5e-16_(1) * Metre,
      /*expected_value_absolute_error=*/2.9e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/2.2e-14_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_9) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 9>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/4.8e-16_(1) * Metre,
      /*expected_value_absolute_error=*/2.9e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/2.2e-14_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_10) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 10>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/4.8e-16_(1) * Metre,
      /*expected_value_absolute_error=*/4.3e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/2.9e-14_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_11) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 11>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/1.1e-16_(1) * Metre,
      /*expected_value_absolute_error=*/3.6e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/2.2e-14_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_12) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 12>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/7.3e-16_(1) * Metre,
      /*expected_value_absolute_error=*/2.9e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/2.9e-14_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_13) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 13>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/1.3e-15_(1) * Metre,
      /*expected_value_absolute_error=*/1.4e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/6.1e-14_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_14) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 14>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/2.1e-16_(1) * Metre,
      /*expected_value_absolute_error=*/1.4e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/6.8e-14_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_15) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 15>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/3.2e-16_(1) * Metre,
      /*expected_value_absolute_error=*/2.9e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/3.5e-13_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_16) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 16>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/1.1e-15_(1) * Metre,
      /*expected_value_absolute_error=*/2.9e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/8.5e-13_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInЧебышёвBasis_2_17) {
  CheckNewhallApproximationErrors<ЧебышёвAdapter, Length, 17>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/4.8e-15_(1) * Metre,
      /*expected_value_absolute_error=*/7.2e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/1.2e-12_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_3) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 3>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/3.9e1_(1) * Metre,
      /*expected_value_absolute_error=*/1.7e2_(1) * Metre,
      /*expected_variation_absolute_error=*/2.3e2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_4) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 4>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/1.8e1_(1) * Metre,
      /*expected_value_absolute_error=*/4.7e1_(1) * Metre,
      /*expected_variation_absolute_error=*/1.2e2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_5) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 5>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/8.3e0_(1) * Metre,
      /*expected_value_absolute_error=*/4.3e1_(1) * Metre,
      /*expected_variation_absolute_error=*/1.2e2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_6) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 6>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/1.4e0_(1) * Metre,
      /*expected_value_absolute_error=*/3.8e1_(1) * Metre,
      /*expected_variation_absolute_error=*/1.0e2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_7) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 7>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/7.5e0_(1) * Metre,
      /*expected_value_absolute_error=*/1.4e1_(1) * Metre,
      /*expected_variation_absolute_error=*/4.5e1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_8) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 8>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/4.2e0_(1) * Metre,
      /*expected_value_absolute_error=*/6.3e0_(1) * Metre,
      /*expected_variation_absolute_error=*/2.8e1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_9) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 9>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/3.1e-1_(1) * Metre,
      /*expected_value_absolute_error=*/4.9e0_(1) * Metre,
      /*expected_variation_absolute_error=*/2.2e1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_10) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 10>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/7.5e-1_(1) * Metre,
      /*expected_value_absolute_error=*/6.5e-1_(1) * Metre,
      /*expected_variation_absolute_error=*/3.6e0_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_11) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 11>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/1.8e-1_(1) * Metre,
      /*expected_value_absolute_error=*/2.0e-1_(1) * Metre,
      /*expected_variation_absolute_error=*/1.6e0_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_12) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 12>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/4.9e-2_(1) * Metre,
      /*expected_value_absolute_error=*/7.9e-2_(1) * Metre,
      /*expected_variation_absolute_error=*/7.3e-1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_13) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 13>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/2.6e-2_(1) * Metre,
      /*expected_value_absolute_error=*/1.3e-2_(1) * Metre,
      /*expected_variation_absolute_error=*/1.2e-1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_14) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 14>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/8.9e-4_(1) * Metre,
      /*expected_value_absolute_error=*/1.5e-2_(1) * Metre,
      /*expected_variation_absolute_error=*/1.5e-1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_15) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 15>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/1.6e-3_(1) * Metre,
      /*expected_value_absolute_error=*/4.3e-3_(1) * Metre,
      /*expected_variation_absolute_error=*/4.4e-2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_16) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 16>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/2.9e-4_(1) * Metre,
      /*expected_value_absolute_error=*/1.7e-3_(1) * Metre,
      /*expected_variation_absolute_error=*/1.8e-2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_1_17) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 17>(
      length_function_1_,
      speed_function_1_,
      /*expected_value_error_estimate=*/3.7e-5_(1) * Metre,
      /*expected_value_absolute_error=*/7.6e-4_(1) * Metre,
      /*expected_variation_absolute_error=*/8.2e-3_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_3) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 3>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/7.9e-1_(1) * Metre,
      /*expected_value_absolute_error=*/2.0e0_(1) * Metre,
      /*expected_variation_absolute_error=*/1.8e0_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_4) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 4>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/2.5e-1_(1) * Metre,
      /*expected_value_absolute_error=*/2.9e-1_(1) * Metre,
      /*expected_variation_absolute_error=*/4.6e-1_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_5) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 5>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/5.7e-2_(1) * Metre,
      /*expected_value_absolute_error=*/3.6e-2_(1) * Metre,
      /*expected_variation_absolute_error=*/7.4e-2_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_6) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 6>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/8.6e-3_(1) * Metre,
      /*expected_value_absolute_error=*/2.3e-3_(1) * Metre,
      /*expected_variation_absolute_error=*/6.0e-3_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_7) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 7>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/6.2e-4_(1) * Metre,
      /*expected_value_absolute_error=*/2.9e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/2.5e-14_(1) * Metre / Second,
      /*expected_value_absolute_error_with_fma=*/1.7e-14_(1) * Metre);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_8) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 8>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/4.5e-16_(1) * Metre,
      /*expected_value_absolute_error=*/3.9e-14_(1) * Metre,
      /*expected_variation_absolute_error=*/1.2e-13_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_9) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 9>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/4.8e-16_(1) * Metre,
      /*expected_value_absolute_error=*/1.0e-13_(1) * Metre,
      /*expected_variation_absolute_error=*/4.8e-13_(1) * Metre / Second,
      /*expected_value_absolute_error_with_fma=*/1.1e-13_(1) * Metre);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_10) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 10>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/4.8e-16_(1) * Metre,
      /*expected_value_absolute_error=*/3.5e-13_(1) * Metre,
      /*expected_variation_absolute_error=*/7.1e-13_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_11) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 11>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/1.1e-16_(1) * Metre,
      /*expected_value_absolute_error=*/5.6e-13_(1) * Metre,
      /*expected_variation_absolute_error=*/1.7e-12_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_12) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 12>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/7.3e-16_(1) * Metre,
      /*expected_value_absolute_error=*/1.3e-11_(1) * Metre,
      /*expected_variation_absolute_error=*/5.4e-11_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_13) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 13>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/1.3e-15_(1) * Metre,
      /*expected_value_absolute_error=*/1.9e-11_(1) * Metre,
      /*expected_variation_absolute_error=*/9.5e-11_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_14) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 14>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/2.1e-16_(1) * Metre,
      /*expected_value_absolute_error=*/8.1e-12_(1) * Metre,
      /*expected_variation_absolute_error=*/5.1e-11_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_15) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 15>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/3.2e-16_(1) * Metre,
      /*expected_value_absolute_error=*/3.0e-10_(1) * Metre,
      /*expected_variation_absolute_error=*/1.6e-9_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_16) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 16>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/1.1e-15_(1) * Metre,
      /*expected_value_absolute_error=*/5.8e-10_(1) * Metre,
      /*expected_variation_absolute_error=*/2.8e-9_(1) * Metre / Second);
}

TEST_F(NewhallTest, ApproximationInMonomialBasis_2_17) {
  CheckNewhallApproximationErrors<MonomialAdapter, Length, 17>(
      length_function_2_,
      speed_function_2_,
      /*expected_value_error_estimate=*/4.8e-15_(1) * Metre,
      /*expected_value_absolute_error=*/3.7e-9_(1) * Metre,
      /*expected_variation_absolute_error=*/2.2e-8_(1) * Metre / Second);
}

TEST_F(NewhallTest, Affine) {
  auto instant_function = [](Instant const t) -> Instant {
    return t;
  };
  auto double_function = [](Instant const t) -> double {
    return 1;
  };

  CheckNewhallApproximationErrors<MonomialAdapter, Instant, 3>(
      instant_function,
      double_function,
      /*expected_value_error_estimate=*/0.0_(1) * Second,
      /*expected_value_absolute_error=*/1.1e-16_(1) * Second,
      /*expected_variation_absolute_error=*/0.0_(1));
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

    EXPECT_THAT(RelativeError((*approximation)(t_min_),
                              length_function_1_(t_min_)), IsNear(9e-13_(1)));
}

}  // namespace numerics
}  // namespace principia
