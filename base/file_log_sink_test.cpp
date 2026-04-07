
#include "base/file_log_sink.hpp"

#include "absl/log/initialize.h"
#include "absl/log/log.h"
#include "absl/log/log_sink_registry.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using ::testing::HasSubstr;
using ::testing::Not;
using namespace principia::base::file_log_sink;

class FileLogSinkTest : public ::testing::Test {
protected:
  FileLogSinkTest()
      : j2000_(absl::FromDateTime(2000, 1, 1, 11, 58, 55, absl::UTCTimeZone()) +
               absl::Milliseconds(816)) {
    for (auto const severity : absl::LogSeverities()) {
      std::filesystem::path log = std::string("FileLogSinkTest") +
                                  absl::LogSeverityName(severity) +
                                  ".20000101T115855Z.log";
      if (std::filesystem::exists(log)) {
        std::filesystem::remove(log);
      }
    }
  }

  absl::Time const j2000_;
};

using FileLogSinkDeathTest = FileLogSinkTest;

TEST_F(FileLogSinkDeathTest, Flush) {
  EXPECT_DEATH(
      {
        absl::InitializeLog();
        static FileLogSink* sink = new FileLogSink("FileLogSinkTest.", ".log");
        absl::AddLogSink(sink);
        LOG(INFO).WithTimestamp(j2000_).WithThreadID(1729) << "info";
        LOG(WARNING).WithTimestamp(j2000_).WithThreadID(1729) << "warning";
        LOG(ERROR).WithTimestamp(j2000_).WithThreadID(1729) << "error";
        std::terminate();
      },
      "error");
  {
    std::ifstream info_file("FileLogSinkTest.INFO.20000101T115855Z.log");
    std::string info_log;
    info_log.assign(std::istreambuf_iterator(info_file),
                    std::istreambuf_iterator<char>());
    EXPECT_THAT(info_log, Not(HasSubstr("error")));
  }
  {
    std::ifstream warning_file("FileLogSinkTest.WARNING.20000101T115855Z.log");
    std::string warning_log;
    warning_log.assign(std::istreambuf_iterator(warning_file),
                       std::istreambuf_iterator<char>());
    EXPECT_THAT(warning_log, HasSubstr("error"));
  }
}

TEST_F(FileLogSinkTest, RelativePath) {
  absl::InitializeLog();
  static FileLogSink* sink =
      new FileLogSink("FileLogSinkTest.", ".log");
  absl::AddLogSink(sink);
  LOG(INFO).WithTimestamp(j2000_).WithThreadID(1729) << "info";
  LOG(WARNING).WithTimestamp(j2000_).WithThreadID(1729) << "warning";
  LOG(ERROR).WithTimestamp(j2000_).WithThreadID(1729) << "error";
}

}  // namespace base
}  // namespace principia
