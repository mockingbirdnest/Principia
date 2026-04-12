#pragma once

#include <array>
#include <filesystem>
#include <fstream>
#include <string>

#include "absl/log/log_sink.h"

namespace principia {
namespace base {
namespace _file_log_sink {

class FileLogSink : public absl::LogSink {
 public:
  FileLogSink(std::filesystem::path path, std::string extension);

  void Send(absl::LogEntry const& entry) override;
  void Flush() override;

  // Files for this severity and lower are not immediately flushed.
  // By default, this is set to `absl::LogSeverityAtMost::kINFO`, so that
  // INFO logs are buffered, and WARNING and higher logs are immediately
  // flushed.
  // set_buffered_level(absl::LogSeverityAtMost::kNegativeInfinity) causes
  // INFO logs to be flushed.
  // Note that higher severity entries do not cause lower severity files to be
  // flushed: if buffered_level is `absl::LogSeverityAtMost::kINFO`, a log
  // WARNING will cause the WARNING file to be flushed, but not the INFO file.
  absl::LogSeverityAtMost buffered_level() const;
  void set_buffered_level(absl::LogSeverityAtMost severity);

 private:
  std::array<std::ofstream, absl::LogSeverities().size()> files_;
  std::filesystem::path const path_;
  std::string const extension_;
  absl::LogSeverityAtMost buffered_level_ = absl::LogSeverityAtMost::kInfo;
};

}  // namespace _file_log_sink
}  // namespace base
}  // namespace principia
