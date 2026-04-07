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

 private:
  std::array<std::ofstream, absl::LogSeverities().size()> file_;
  std::filesystem::path const path_;
  std::string const extension_;
};

}  // namespace file_log_sink
}  // namespace base
}  // namespace principia
