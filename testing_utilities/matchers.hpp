#pragma once

#include <pmmintrin.h>

#include <string>

#include "base/status.hpp"
#include "gmock/gmock.h"
#include "google/protobuf/message.h"
#include "google/protobuf/util/message_differencer.h"
#include "quantities/quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace testing_utilities {

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
  return arg.ok();
}

MATCHER_P(StatusIs,
          error,
          std::string(negation ? "does not have" : "has") + " error: " +
              ::principia::base::ErrorToString(error)) {
  return arg.error() == error;
}

MATCHER_P(SSEHighHalfIs,
          value,
          std::string(negation ? "does not have" : "has") + " high half: " +
              ::principia::quantities::DebugString(value)) {
  return _mm_cvtsd_f64(_mm_unpackhi_pd(arg, arg)) == value;
}

MATCHER_P(SSELowHalfIs,
          value,
          std::string(negation ? "does not have" : "has") + " low half: " +
              ::principia::quantities::DebugString(value)) {
  return _mm_cvtsd_f64(arg) == value;
}

}  // namespace testing_utilities

namespace serialization {

inline void PrintTo(Part const& message, std::ostream* const os) {
  *os << message.ShortDebugString();
}

}  // namespace serialization

}  // namespace principia
