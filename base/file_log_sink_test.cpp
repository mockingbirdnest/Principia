
#include "base/file_log_sink.hpp"

#include "absl/log/initialize.h"
#include "absl/log/log.h"
#include "absl/log/log_sink_registry.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using ::testing::HasSubstr;
using ::testing::MatchesRegex;
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
    std::string const info_log{std::istreambuf_iterator(info_file), {}};
    EXPECT_THAT(info_log, Not(HasSubstr("error")));
  }
  {
    std::ifstream warning_file("FileLogSinkTest.WARNING.20000101T115855Z.log");
    std::string const warning_log{std::istreambuf_iterator(warning_file), {}};
    EXPECT_THAT(warning_log, HasSubstr("error"));
  }
  {
    std::ifstream error_file("FileLogSinkTest.ERROR.20000101T115855Z.log");
    std::string const error_log{std::istreambuf_iterator(error_file), {}};
    EXPECT_THAT(error_log, HasSubstr("error"));
  }
}

TEST_F(FileLogSinkDeathTest, SetBufferedLevel) {
  EXPECT_DEATH(
      {
        absl::InitializeLog();
        static FileLogSink* sink = new FileLogSink("FileLogSinkTest.", ".log");
        absl::AddLogSink(sink);
        sink->set_buffered_level(absl::LogSeverity::kWarning);
        LOG(INFO).WithTimestamp(j2000_).WithThreadID(1729) << "info";
        LOG(WARNING).WithTimestamp(j2000_).WithThreadID(1729) << "warning";
        LOG(ERROR).WithTimestamp(j2000_).WithThreadID(1729) << "error";
        std::terminate();
      },
      "error");
  {
    std::ifstream info_file("FileLogSinkTest.INFO.20000101T115855Z.log");
    std::string const info_log{std::istreambuf_iterator(info_file), {}};
    EXPECT_THAT(info_log, Not(HasSubstr("error")));
  }
  {
    std::ifstream warning_file("FileLogSinkTest.WARNING.20000101T115855Z.log");
    std::string const warning_log{std::istreambuf_iterator(warning_file), {}};
    EXPECT_THAT(warning_log, Not(HasSubstr("error")));
  }
  {
    std::ifstream error_file("FileLogSinkTest.ERROR.20000101T115855Z.log");
    std::string const error_log{std::istreambuf_iterator(error_file), {}};
    EXPECT_THAT(error_log, HasSubstr("error"));
  }
}

TEST_F(FileLogSinkTest, Nonfatal) {
  absl::InitializeLog();
  static FileLogSink* sink =
      new FileLogSink("FileLogSinkTest.", ".log");
  absl::AddLogSink(sink);
  LOG(INFO).WithTimestamp(j2000_).WithThreadID(1729) << "info";
  LOG(WARNING).WithTimestamp(j2000_).WithThreadID(1729) << "warning";
  LOG(ERROR).WithTimestamp(j2000_).WithThreadID(1729) << "error";
  absl::FlushLogSinks();
  {
    std::ifstream info_file("FileLogSinkTest.INFO.20000101T115855Z.log");
    std::string const info_log{std::istreambuf_iterator(info_file), {}};
    EXPECT_THAT(info_log, MatchesRegex(R"(^Log file created at: 20000101T115855Z
Running on machine: .*
Log line format: \[IWEF]mmdd hh:mm:ss.μμμμμμ threadid file:line] msg
I0101 11:58:55.816000    1729 file_log_sink_test.cpp:102] info
W0101 11:58:55.816000    1729 file_log_sink_test.cpp:103] warning
E0101 11:58:55.816000    1729 file_log_sink_test.cpp:104] error
$)"));
  }
  {
    std::ifstream warning_file("FileLogSinkTest.WARNING.20000101T115855Z.log");
    std::string const warning_log{std::istreambuf_iterator(warning_file), {}};
    EXPECT_THAT(warning_log, MatchesRegex(R"(^Log file created at: 20000101T115855Z
Running on machine: .*
Log line format: \[IWEF]mmdd hh:mm:ss.μμμμμμ threadid file:line] msg
W0101 11:58:55.816000    1729 file_log_sink_test.cpp:103] warning
E0101 11:58:55.816000    1729 file_log_sink_test.cpp:104] error
$)"));
  }
  {
    std::ifstream error_file("FileLogSinkTest.ERROR.20000101T115855Z.log");
    std::string const error_log{std::istreambuf_iterator(error_file), {}};
    EXPECT_THAT(error_log, MatchesRegex(R"(^Log file created at: 20000101T115855Z
Running on machine: .*
Log line format: \[IWEF]mmdd hh:mm:ss.μμμμμμ threadid file:line] msg
E0101 11:58:55.816000    1729 file_log_sink_test.cpp:104] error
$)"));
  }
}

}  // namespace base
}  // namespace principia
