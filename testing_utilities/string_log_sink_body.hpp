#pragma once

#include "testing_utilities/string_log_sink.hpp"

#include <string>

namespace principia {
namespace testing_utilities {

StringLogSink::StringLogSink(google::LogSeverity const minimal_severity)
    : minimal_severity_(minimal_severity) {
  google::AddLogSink(this);
}

StringLogSink::~StringLogSink() {
  google::RemoveLogSink(this);
}

void StringLogSink::send(google::LogSeverity const severity,
                         char const* const full_filename,
                         char const* const base_filename,
                         int const line,
                         tm const* const tm_time,
                         const char* const message,
                         std::size_t const message_len) {
  if (severity < minimal_severity_) {
    return;
  }
  absl::MutexLock lock(&mutex_);
  absl::StrAppend(
      &string_,
      ToString(severity, base_filename, line, tm_time, message, message_len));
}

std::string& StringLogSink::string() {
  return string_;
}

}  // namespace testing_utilities
}  // namespace principia
