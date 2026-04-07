#pragma once

#include <array>
#include <filesystem>
#include <fstream>

#include "absl/log/log_sink.h"

namespace principia {
namespace base {
namespace file_log_sink {

class FileLogSink : public absl::LogSink {
 public:
  FileLogSink(std::filesystem::path path, std::string extension);

  void Send(const absl::LogEntry& entry) override;
  void Flush() override;

  // Files for this severity and lower are not immediately flushed.
  // This can be set to values for which there is no enumerator:
  // set_buffered_level(-1) causes INFO (severity=0) logs to be flushed.
  // Note that higher severity entries do not cause lower severity files to be
  // flushed: if buffered_level is INFO, a log WARNING will cause the WARNING
  // file to be flushed, but not the INFO file.
  absl::LogSeverity buffered_level() const;
  void set_buffered_level(absl::LogSeverity severity);

 private:
  std::array<std::ofstream, absl::LogSeverities().size()> files_;
  std::filesystem::path const path_;
  std::string const extension_;
  absl::LogSeverity buffered_level_ = absl::LogSeverity::kInfo;
};

}  // namespace file_log_sink
}  // namespace base
}  // namespace principia
