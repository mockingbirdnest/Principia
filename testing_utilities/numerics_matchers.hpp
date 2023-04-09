#pragma once

#include "gmock/gmock.h"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace _numerics_matchers {
namespace internal {

using ::testing::Matcher;
using namespace principia::quantities::_named_quantities;

template<typename Value>
Matcher<Value> DifferenceFrom(
    Value const& expected,
    Matcher<Difference<Value>> const& error_matcher);

template<typename Value, typename ErrorMatcher>
Matcher<Value> AbsoluteErrorFrom(Value const& expected,
                                 ErrorMatcher const& error_matcher);

template<typename Value>
Matcher<Value> RelativeErrorFrom(Value const& expected,
                                 Matcher<double> const& error_matcher);

}  // namespace internal

using internal::AbsoluteErrorFrom;
using internal::DifferenceFrom;
using internal::RelativeErrorFrom;

}  // namespace _numerics_matchers
}  // namespace testing_utilities
}  // namespace principia

namespace principia::testing_utilities {
using namespace principia::testing_utilities::_numerics_matchers;
}  // namespace principia::testing_utilities

#include "testing_utilities/numerics_matchers_body.hpp"
