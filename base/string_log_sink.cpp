#include "base/string_log_sink.hpp"

#include <exception>
#include <string>
#include <vector>
#include <utility>

#include "absl/log/check.h"
#include "absl/time/time.h"

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
