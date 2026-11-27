// Body of link:base/string_log_sink.hpp

module principia.base.string_log_sink;

import std;

import "absl/log/check.h";
import "absl/time/time.h";

namespace principia {
namespace base {
namespace _string_log_sink {
namespace internal {

void StringLogSink::Send(absl::LogEntry const& entry) {
  absl::MutexLock l(lock_);
  // FATAL messages are sent twice, first without and then with the
  // stacktrace.  Log the message only once.
  if (entry.stacktrace().empty()) {
    logs_.emplace_back(entry.text_message_with_prefix_and_newline());
  } else {
    logs_.emplace_back(entry.stacktrace());
  }
}

std::vector<std::string> StringLogSink::logs() const {
  absl::ReaderMutexLock l(lock_);
  return logs_;
}

}  // namespace internal
}  // namespace _string_log_sink
}  // namespace base
}  // namespace principia
