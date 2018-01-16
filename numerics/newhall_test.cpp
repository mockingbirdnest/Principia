
#include "numerics/newhall.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace numerics {
namespace newhall {

using geometry::Instant;
using quantities::Length;
using quantities::Speed;
using quantities::si::Metre;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using ::testing::AllOf;
using ::testing::ElementsAre;
using ::testing::Gt;
using ::testing::Lt;

class NewhallTest : public ::testing::Test {
 protected:
  NewhallTest()
      : t_min_(t0_ - 1 * Second),
        t_max_(t0_ + 3 * Second) {}

  void NewhallApproximationErrors(
      std::function<Length(Instant const)> length_function,
      std::function<Speed(Instant const)> speed_function,
      std::vector<Length>& length_absolute_errors,
      std::vector<Speed>& speed_absolute_errors) {
    std::vector<Length> lengths;
    std::vector<Speed> speeds;
    for (Instant t = t_min_; t <= t_max_; t += 0.5 * Second) {
      lengths.push_back(length_function(t));
      speeds.push_back(speed_function(t));
    }

    length_absolute_errors.clear();
    speed_absolute_errors.clear();
    for (int degree = 3; degree <= 17; ++degree) {
      ЧебышёвSeries<Length> const approximation =
          ApproximationInЧебышёвBasis<Length>(
              degree, lengths, speeds, t_min_, t_max_);

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
      length_absolute_errors.push_back(length_absolute_error);
      speed_absolute_errors.push_back(speed_absolute_error);

      // Check the conditions at the bounds.
      EXPECT_THAT(approximation.Evaluate(t_min_),
                  AlmostEquals(length_function(t_min_), 0, 248));
      EXPECT_THAT(approximation.Evaluate(t_max_),
                  AlmostEquals(length_function(t_max_), 0, 3));
      EXPECT_THAT(approximation.EvaluateDerivative(t_min_),
                  AlmostEquals(speed_function(t_min_), 1, 1185));
      EXPECT_THAT(approximation.EvaluateDerivative(t_max_),
                  AlmostEquals(speed_function(t_max_), 0, 339));
    }
  }

  // A helper that splits an array in two chunks and applies distinct matchers
  // to the chunks.  Necessary because ElementsAre only supports 10 elements and
  // ElementsAreArray does not support matchers as arguments.
  template<typename Vector, typename Matcher1, typename Matcher2>
  void ExpectMultipart(std::vector<Vector> const& v,
                       Matcher1 const& matcher1,
                       Matcher2 const& matcher2) {
    std::vector<Vector> v_0_9(10);
    std::copy(v.begin(), v.begin() + 10, v_0_9.begin());
    EXPECT_THAT(v_0_9, matcher1);
    std::vector<Vector> v_10_end(v.size() - 10);
    std::copy(v.begin() + 10, v.end(), v_10_end.begin());
    EXPECT_THAT(v_10_end, matcher2);
  }

  Instant const t0_;
  Instant t_min_;
  Instant t_max_;
};

TEST_F(NewhallTest, ApproximationInЧебышёвBasis) {
  std::vector<Length> length_absolute_errors;
  std::vector<Speed> speed_absolute_errors;

  auto near_length = [](Length const& length) {
    return AllOf(Gt(0.9 * length), Lt(length));
  };
  auto near_speed = [](Speed const& speed) {
    return AllOf(Gt(0.9 * speed), Lt(speed));
  };

  {
    auto length_function = [this](Instant const t) -> Length {
      return 0.5 * Metre + 2 * Metre *
             std::sin((t - t_min_) / (0.3 * Second)) *
             std::exp((t - t_min_) / (1 * Second));
    };
    auto speed_function = [this](Instant const t) -> Speed {
      return ((2 * Metre) / (0.3 * Second) *
                   std::cos((t - t_min_) / (0.3 * Second)) +
              (2 * Metre / Second) *
                   std::sin((t - t_min_) / (0.3 * Second))) *
                       std::exp((t - t_min_) / (1 * Second));
    };

    NewhallApproximationErrors(length_function,
                               speed_function,
                               length_absolute_errors,
                               speed_absolute_errors);
  }

  ExpectMultipart(length_absolute_errors,
                  ElementsAre(near_length(1.7e2 * Metre),
                              near_length(4.7e1 * Metre),
                              near_length(4.3e1 * Metre),
                              near_length(3.8e1 * Metre),
                              near_length(1.5e1 * Metre),
                              near_length(6.3 * Metre),
                              near_length(4.9 * Metre),
                              near_length(6.5e-1 * Metre),
                              near_length(2.0e-1 * Metre),
                              near_length(7.9e-2 * Metre)),
                  ElementsAre(near_length(1.3e-2 * Metre),
                              near_length(1.6e-2 * Metre),
                              near_length(4.3e-3 * Metre),
                              near_length(1.7e-3 * Metre),
                              near_length(7.6e-4 * Metre)));
  ExpectMultipart(speed_absolute_errors,
                  ElementsAre(near_speed(2.3e2 * Metre / Second),
                              near_speed(1.3e2 * Metre / Second),
                              near_speed(1.2e2 * Metre / Second),
                              near_speed(1.1e2 * Metre / Second),
                              near_speed(4.5e1 * Metre / Second),
                              near_speed(2.8e1 * Metre / Second),
                              near_speed(2.2e1 * Metre / Second),
                              near_speed(3.6 * Metre / Second),
                              near_speed(1.6 * Metre / Second),
                              near_speed(7.3e-1 * Metre / Second)),
                  ElementsAre(near_speed(1.3e-1 * Metre / Second),
                              near_speed(1.5e-1 * Metre / Second),
                              near_speed(4.4e-2 * Metre / Second),
                              near_speed(1.8e-2 * Metre / Second),
                              near_speed(8.2e-3 * Metre / Second)));

  {
    auto length_function = [this](Instant const t) -> Length {
      return 5 * Metre * (1 + (t - t_min_) / (0.3 * Second) +
                          std::pow((t - t_min_) / (4 * Second), 7));
    };
    auto speed_function = [this](Instant const t) -> Speed {
      return 5 * Metre * (1 / (0.3 * Second) +
                          (7 / (4 * Second)) *
                              std::pow((t - t_min_) / (4 * Second), 6));
    };

    NewhallApproximationErrors(length_function,
                               speed_function,
                               length_absolute_errors,
                               speed_absolute_errors);
  }

  ExpectMultipart(length_absolute_errors,
                  ElementsAre(near_length(2.0 * Metre),
                              near_length(2.9e-1 * Metre),
                              near_length(3.6e-2 * Metre),
                              near_length(2.3e-3 * Metre),
                              near_length(2.9e-14 * Metre),
                              near_length(2.9e-14 * Metre),
                              near_length(2.9e-14 * Metre),
                              near_length(4.3e-14 * Metre),
                              near_length(3.6e-14 * Metre),
                              near_length(2.9e-14 * Metre)),
                  ElementsAre(near_length(1.5e-14 * Metre),
                              near_length(1.5e-14 * Metre),
                              near_length(2.9e-14 * Metre),
                              near_length(2.9e-14 * Metre),
                              near_length(7.2e-14 * Metre)));
  ExpectMultipart(speed_absolute_errors,
                  ElementsAre(near_speed(1.8 * Metre / Second),
                              near_speed(4.6e-1 * Metre / Second),
                              near_speed(7.4e-2 * Metre / Second),
                              near_speed(6.0e-3 * Metre / Second),
                              near_speed(2.5e-14 * Metre / Second),
                              near_speed(2.2e-14 * Metre / Second),
                              near_speed(2.2e-14 * Metre / Second),
                              near_speed(2.9e-14 * Metre / Second),
                              near_speed(2.2e-14 * Metre / Second),
                              near_speed(2.9e-14 * Metre / Second)),
                  ElementsAre(AllOf(Gt(5.3e-14 * Metre / Second),
                                    Lt(6.1e-14 * Metre / Second)),
                              near_speed(6.8e-14 * Metre / Second),
                              near_speed(3.5e-13 * Metre / Second),
                              near_speed(8.5e-13 * Metre / Second),
                              near_speed(1.3e-12 * Metre / Second)));
}

}  // namespace newhall
}  // namespace numerics
}  // namespace principia
