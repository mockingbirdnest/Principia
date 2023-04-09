#pragma once

#include <cstdint>
#include <string>

#include "absl/strings/str_cat.h"
#include "absl/synchronization/mutex.h"
#include "glog/logging.h"

namespace principia {
namespace testing_utilities {
namespace _string_log_sink {
namespace internal {

class StringLogSink : google::LogSink {
 public:
  inline explicit StringLogSink(google::LogSeverity const minimal_severity);

  inline ~StringLogSink();

  inline void send(google::LogSeverity const severity,
                   char const* const full_filename,
                   char const* const base_filename,
                   int const line,
                   tm const* const tm_time,
                   const char* const message,
                   std::size_t const message_len) override;

  inline std::string& string();

 private:
  google::LogSeverity const minimal_severity_;
  absl::Mutex mutex_;
  std::string string_ GUARDED_BY(mutex_);
};

}  // namespace internal

using internal::StringLogSink;

}  // namespace _string_log_sink
}  // namespace testing_utilities
}  // namespace principia

namespace principia::testing_utilities {
using namespace principia::testing_utilities::_string_log_sink;
}  // namespace principia::testing_utilities

#include "testing_utilities/string_log_sink_body.hpp"
