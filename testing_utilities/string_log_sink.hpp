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
  explicit StringLogSink(google::LogSeverity minimal_severity);

  ~StringLogSink() override;

  void send(google::LogSeverity severity,
            char const* full_filename,
            char const* base_filename,
            int line,
            tm const* tm_time,
            const char* message,
            std::size_t message_len) override;

  std::string const& string() const;

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

#include "testing_utilities/string_log_sink_body.hpp"
