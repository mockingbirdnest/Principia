#include "base/file_log_sink.hpp"

#include "absl/log/initialize.h"
#include "absl/log/log.h"
#include "re2/regexp.h"
#include "absl/log/log_sink_registry.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using ::testing::ContainsRegex;
using ::testing::HasSubstr;
using ::testing::MatchesRegex;
using ::testing::Not;
using ::testing::SizeIs;
using namespace principia::base::_file_log_sink;

class FileLogSinkTest : public ::testing::Test {
 protected:
  FileLogSinkTest()
      : j2000_(absl::FromDateTime(2000, 1, 1, 11, 58, 55, absl::UTCTimeZone()) +
               absl::Milliseconds(816)) {
    RemoveAllLogFiles();
  }

  void RemoveAllLogFiles() {
    for (auto const& entry : std::filesystem::directory_iterator(".")) {
      if (entry.is_regular_file() &&
          entry.path().filename().string().starts_with("FileLogSinkTest.") &&
          entry.path().extension() == ".log") {
        std::filesystem::remove(entry);
      }
    }
  }

  std::filesystem::path UniqueLogFile(absl::LogSeverity const severity) {
    const std::string prefix =
        std::string("FileLogSinkTest.") + absl::LogSeverityName(severity);
    std::vector<std::filesystem::path> candidates;
    for (auto const& entry : std::filesystem::directory_iterator(".")) {
      if (entry.is_regular_file() &&
          entry.path().filename().string().starts_with(prefix) &&
          entry.path().extension() == ".log") {
        candidates.push_back(entry.path());
      }
    }
    EXPECT_THAT(candidates, SizeIs(1));
    EXPECT_THAT(
        candidates.front().filename().string(),
        MatchesRegex(
            R"(FileLogSinkTest\.(INFO|WARNING|ERROR|FATAL)\.20000101-\d\d\d855\.\d+\.log)"));
    return candidates.front();
  }

  absl::Time const j2000_;
};

using FileLogSinkDeathTest = FileLogSinkTest;

// Check that files above the default buffered level of INFO are flushed.
TEST_F(FileLogSinkDeathTest, Flush) {
  EXPECT_DEATH(
      {
        static FileLogSink* sink = new FileLogSink("FileLogSinkTest.", ".log");
        absl::AddLogSink(sink);
        LOG(INFO).WithTimestamp(j2000_).WithThreadID(1729) << "info";
        LOG(WARNING).WithTimestamp(j2000_).WithThreadID(1729) << "warning";
        LOG(ERROR).WithTimestamp(j2000_).WithThreadID(1729) << "error";
        std::terminate();
      },
      "error");
  {
    std::ifstream info_file(UniqueLogFile(absl::LogSeverity::kInfo));
    std::string const info_log{std::istreambuf_iterator(info_file), {}};
    EXPECT_THAT(info_log, Not(HasSubstr("error")));
  }
  {
    std::ifstream warning_file(UniqueLogFile(absl::LogSeverity::kWarning));
    std::string const warning_log{std::istreambuf_iterator(warning_file), {}};
    EXPECT_THAT(warning_log, HasSubstr("error"));
  }
  {
    std::ifstream error_file(UniqueLogFile(absl::LogSeverity::kError));
    std::string const error_log{std::istreambuf_iterator(error_file), {}};
    EXPECT_THAT(error_log, HasSubstr("error"));
  }
}

TEST_F(FileLogSinkDeathTest, SetBufferedLevel) {
  EXPECT_DEATH(
      {
        static FileLogSink* sink = new FileLogSink("FileLogSinkTest.", ".log");
        absl::AddLogSink(sink);
        sink->set_buffered_level(absl::LogSeverityAtMost::kWarning);
        LOG(INFO).WithTimestamp(j2000_).WithThreadID(1729) << "info";
        LOG(WARNING).WithTimestamp(j2000_).WithThreadID(1729) << "warning";
        LOG(ERROR).WithTimestamp(j2000_).WithThreadID(1729) << "error";
        std::terminate();
      },
      "error");
  {
    std::ifstream info_file(UniqueLogFile(absl::LogSeverity::kInfo));
    std::string const info_log{std::istreambuf_iterator(info_file), {}};
    EXPECT_THAT(info_log, Not(HasSubstr("error")));
  }
  {
    std::ifstream warning_file(UniqueLogFile(absl::LogSeverity::kWarning));
    std::string const warning_log{std::istreambuf_iterator(warning_file), {}};
    EXPECT_THAT(warning_log, Not(HasSubstr("error")));
  }
  {
    std::ifstream error_file(UniqueLogFile(absl::LogSeverity::kError));
    std::string const error_log{std::istreambuf_iterator(error_file), {}};
    EXPECT_THAT(error_log, HasSubstr("error"));
  }
}

TEST_F(FileLogSinkDeathTest, Fatal) {
  EXPECT_DEATH(
      {
        static FileLogSink* sink = new FileLogSink("FileLogSinkTest.", ".log");
        absl::AddLogSink(sink);
        LOG(INFO).WithTimestamp(j2000_).WithThreadID(1729) << "info";
        LOG(WARNING).WithTimestamp(j2000_).WithThreadID(1729) << "warning";
        LOG(ERROR).WithTimestamp(j2000_).WithThreadID(1729) << "error";
        LOG(FATAL).WithTimestamp(j2000_).WithThreadID(1729) << "fatal";
      },
      "fatal");
  {
    std::ifstream info_file(UniqueLogFile(absl::LogSeverity::kInfo));
    std::string const info_log{std::istreambuf_iterator(info_file), {}};
    EXPECT_THAT(info_log, ContainsRegex(R"(^Log file created at: 2000-01-01T\d\d:\d8:55.816[+-]\d\d:\d0
Running on machine: .*
Log line format: \[IWEF]mmdd hh:mm:ss.μμμμμμ threadid file:line] msg
I0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] info
W0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] warning
E0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] error
F0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] fatal
\*\*\* Check failure stack trace: \*\*\*
    @   [0-9A-F]{16}  .* \(.*\.(cc|cpp):\d+\)
)"));
  }
  {
    std::ifstream warning_file(UniqueLogFile(absl::LogSeverity::kWarning));
    std::string const warning_log{std::istreambuf_iterator(warning_file), {}};
    EXPECT_THAT(warning_log, ContainsRegex(R"(^Log file created at: 2000-01-01T\d\d:\d8:55.816[+-]\d\d:\d0
Running on machine: .*
Log line format: \[IWEF]mmdd hh:mm:ss.μμμμμμ threadid file:line] msg
W0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] warning
E0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] error
F0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] fatal
\*\*\* Check failure stack trace: \*\*\*
    @   [0-9A-F]{16}  .* \(.*\.(cc|cpp):\d+\)
)"));
  }
  {
    std::ifstream error_file(UniqueLogFile(absl::LogSeverity::kError));
    std::string const error_log{std::istreambuf_iterator(error_file), {}};
    EXPECT_THAT(error_log, ContainsRegex(R"(^Log file created at: 2000-01-01T\d\d:\d8:55.816[+-]\d\d:\d0
Running on machine: .*
Log line format: \[IWEF]mmdd hh:mm:ss.μμμμμμ threadid file:line] msg
E0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] error
F0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] fatal
\*\*\* Check failure stack trace: \*\*\*
    @   [0-9A-F]{16}  .* \(.*\.(cc|cpp):\d+\)
)"));
  }
  {
    std::ifstream fatal_file(UniqueLogFile(absl::LogSeverity::kFatal));
    std::string const fatal_log{std::istreambuf_iterator(fatal_file), {}};
    EXPECT_THAT(fatal_log, ContainsRegex(R"(^Log file created at: 2000-01-01T\d\d:\d8:55.816[+-]\d\d:\d0
Running on machine: .*
Log line format: \[IWEF]mmdd hh:mm:ss.μμμμμμ threadid file:line] msg
F0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] fatal
\*\*\* Check failure stack trace: \*\*\*
    @   [0-9A-F]{16}  .* \(.*\.(cc|cpp):\d+\)
)"));
  }
}

TEST_F(FileLogSinkTest, Nonfatal) {
  static FileLogSink* sink =
      new FileLogSink("FileLogSinkTest.", ".log");
  absl::AddLogSink(sink);
  LOG(INFO).WithTimestamp(j2000_).WithThreadID(1729) << "info";
  LOG(WARNING).WithTimestamp(j2000_).WithThreadID(1729) << "warning";
  LOG(ERROR).WithTimestamp(j2000_).WithThreadID(1729) << "error";
  absl::FlushLogSinks();
  {
    std::ifstream info_file(UniqueLogFile(absl::LogSeverity::kInfo));
    std::string const info_log{std::istreambuf_iterator(info_file), {}};
    EXPECT_THAT(info_log, MatchesRegex(R"(Log file created at: 2000-01-01T\d\d:\d8:55.816[+-]\d\d:\d0
Running on machine: .*
Log line format: \[IWEF]mmdd hh:mm:ss.μμμμμμ threadid file:line] msg
I0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] info
W0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] warning
E0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] error
)"));
  }
  {
    std::ifstream warning_file(UniqueLogFile(absl::LogSeverity::kWarning));
    std::string const warning_log{std::istreambuf_iterator(warning_file), {}};
    EXPECT_THAT(warning_log, MatchesRegex(R"(Log file created at: 2000-01-01T\d\d:\d8:55.816[+-]\d\d:\d0
Running on machine: .*
Log line format: \[IWEF]mmdd hh:mm:ss.μμμμμμ threadid file:line] msg
W0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] warning
E0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] error
)"));
  }
  {
    std::ifstream error_file(UniqueLogFile(absl::LogSeverity::kError));
    std::string const error_log{std::istreambuf_iterator(error_file), {}};
    EXPECT_THAT(error_log, MatchesRegex(R"(Log file created at: 2000-01-01T\d\d:\d8:55.816[+-]\d\d:\d0
Running on machine: .*
Log line format: \[IWEF]mmdd hh:mm:ss.μμμμμμ threadid file:line] msg
E0101 \d\d:\d8:55.816000    1729 file_log_sink_test.cpp:\d+] error
)"));
  }
}

}  // namespace base
}  // namespace principia
