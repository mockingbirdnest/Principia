export module principia.base.string_log_sink;

import std;

import "absl/log/log_sink.h";
import "absl/synchronization/mutex.h";

namespace principia {
namespace base {
namespace _string_log_sink {
namespace internal {
export {

class StringLogSink : public absl::LogSink {
 public:
  StringLogSink() = default;

  void Send(absl::LogEntry const& entry) override;

  std::vector<std::string> logs() const;

 private:
  mutable absl::Mutex lock_;
  std::vector<std::string> logs_ ABSL_GUARDED_BY(lock_);
};

}  // export
}  // namespace internal

export using internal::StringLogSink;

}  // namespace _string_log_sink
}  // namespace base
}  // namespace principia
