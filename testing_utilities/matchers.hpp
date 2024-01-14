#pragma once

#include <pmmintrin.h>

#include <string>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "gmock/gmock.h"
#include "google/protobuf/message.h"
#include "google/protobuf/util/message_differencer.h"
#include "quantities/quantities.hpp"  // 🧙 For the clients.
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace testing_utilities {
namespace _matchers {
namespace internal {

// This is not defined in base/status_utilities.hpp to avoid pulling gmock in
// non-test code.
#define EXPECT_OK(value) \
  EXPECT_THAT((value), ::principia::testing_utilities::_matchers::IsOk());

inline absl::Status StatusOf(absl::Status const& s) {
  return s;
}

template<typename T>
absl::Status StatusOf(absl::StatusOr<T> const& s) {
  return s.status();
}

MATCHER_P(EqualsProto,
          expected,
          std::string(negation ? "is not" : "is") + " equal to:\n" +
              expected.ShortDebugString()) {
  std::string result;
  ::google::protobuf::util::MessageDifferencer differencer;
  differencer.ReportDifferencesToString(&result);
  if (differencer.Compare(expected, arg)) {
    return true;
  }
  *result_listener << "\nDifference found: " << result;
  return false;
}

MATCHER(IsOk,
        std::string(negation ? "is not" : "is") + " ok") {
  if (!arg.ok()) {
    *result_listener << "Status is " << StatusOf(arg);
    return false;
  }
  return true;
}

MATCHER_P(IsOkAndHolds,
          value_matcher,
          std::string(negation ? "is not" : "is") + " ok and holds a value") {
  return arg.ok() && ::testing::ExplainMatchResult(
                         value_matcher, arg.value(), result_listener);
}

MATCHER_P(StatusIs,
          code,
          std::string(negation ? "does not have" : "has") + " code: " +
              (std::stringstream{} << code).str()) {
  return arg.code() == code;
}

MATCHER_P2(
    IntervalMatches,
    min_matcher, max_matcher,
    "has a min that " +
        ::testing::DescribeMatcher<decltype(std::declval<arg_type>().min)>(
            min_matcher,
            negation) +
        " and a max that " +
        ::testing::DescribeMatcher<decltype(std::declval<arg_type>().max)>(
            max_matcher,
            negation)) {
  bool const min_matches = ::testing::Matches(min_matcher)(arg.min);
  bool const max_matches = ::testing::Matches(max_matcher)(arg.max);
  if (!min_matches) {
    *result_listener << "whose min doesn't match because ";
    ::testing::ExplainMatchResult(min_matcher, arg.min, result_listener);
  }
  if (!max_matches) {
    if (!min_matches) {
      *result_listener << " and ";
    }
    *result_listener << "whose max doesn't match because ";
    ::testing::ExplainMatchResult(max_matcher, arg.max, result_listener);
  }
  return min_matches && max_matches;
}

MATCHER_P(SSEHighHalfIs,
          value,
          std::string(negation ? "does not have" : "has") + " high half: " +
              ::principia::quantities::_quantities::DebugString(value)) {
  return _mm_cvtsd_f64(_mm_unpackhi_pd(arg, arg)) == value;
}

MATCHER_P(SSELowHalfIs,
          value,
          std::string(negation ? "does not have" : "has") + " low half: " +
              ::principia::quantities::_quantities::DebugString(value)) {
  return _mm_cvtsd_f64(arg) == value;
}

}  // namespace internal

using internal::EqualsProto;
using internal::IntervalMatches;
using internal::IsOk;
using internal::IsOkAndHolds;
using internal::SSEHighHalfIs;
using internal::SSELowHalfIs;
using internal::StatusIs;

}  // namespace _matchers
}  // namespace testing_utilities

namespace serialization {

inline void PrintTo(Part const& message, std::ostream* const os) {
  *os << message.ShortDebugString();
}

}  // namespace serialization

}  // namespace principia
