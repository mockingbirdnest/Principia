#pragma once

#include "gmock/gmock.h"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_numerics_matchers {

using quantities::Difference;
using ::testing::Matcher;

template<typename Value>
Matcher<Value const&> DifferenceFrom(
    Value const& expected,
    Matcher<Difference<Value>> const& error_matcher);

template<typename Value, typename ErrorMatcher>
Matcher<Value const&> AbsoluteErrorFrom(Value const& expected,
                                        ErrorMatcher const& error_matcher);

template<typename Value>
Matcher<Value const&> RelativeErrorFrom(Value const& expected,
                                        Matcher<double> const& error_matcher);

}  // namespace internal_numerics_matchers

using internal_numerics_matchers::AbsoluteErrorFrom;
using internal_numerics_matchers::DifferenceFrom;
using internal_numerics_matchers::RelativeErrorFrom;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/numerics_body.hpp"
#include "testing_utilities/numerics_matchers_body.hpp"
