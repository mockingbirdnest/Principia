
#include "base/file_log_sink.hpp"

#include "absl/log/initialize.h"
#include "absl/log/log.h"
#include "absl/log/log_sink_registry.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using namespace principia::base::file_log_sink;

class FileLogSinkTest : public ::testing::Test {
};

TEST_F(FileLogSinkTest, RelativePath) {
  for (auto const severity : absl::LogSeverities()) {
    std::filesystem::path log = std::string("FileLogSinkTest") +
                                absl::LogSeverityName(severity) +
                                ".20000101T115855Z.log";
    if (std::filesystem::exists(log)) {
      std::filesystem::remove(log);
    }
  }
  absl::InitializeLog();
  static FileLogSink* sink =
      new FileLogSink("FileLogSinkTest.", ".log");
  absl::AddLogSink(sink);
  absl::Time const j2000 =
      absl::FromDateTime(2000, 1, 1, 11, 58, 55, absl::UTCTimeZone()) +
      absl::Milliseconds(816);
  LOG(INFO).WithTimestamp(j2000).WithThreadID(1729) << "info";
  LOG(WARNING).WithTimestamp(j2000).WithThreadID(1729) << "warning";
  LOG(ERROR).WithTimestamp(j2000).WithThreadID(1729) << "error";
}

}}
