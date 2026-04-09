#include "base/file_log_sink.hpp"

#include <iostream>

#include "absl/log/check.h"
#include "absl/time/time.h"
#include "base/macros.hpp"  // 🧙 For OS_WIN.

#if OS_WIN
#include <process.h>
#include <windows.h>
#else
#include <sys/utsname.h>
#include <unistd.h>
#endif

namespace principia {
namespace base {
namespace _file_log_sink {

FileLogSink::FileLogSink(std::filesystem::path path, std::string extension)
    : path_(std::move(path)), extension_(std::move(extension)) {
  CHECK_EQ(extension_.front(), '.');
}

void FileLogSink::Send(absl::LogEntry const& entry) {
  for (auto const severity : absl::LogSeverities()) {
    if (severity > entry.log_severity()) {
      break;
    }
    auto& file = files_[static_cast<int>(severity)];
    if (!file.is_open()) {
#if OS_WIN
      std::int32_t const pid = _getpid();
#else
      std::int32_t const pid = getpid();
#endif
      file.open(
          std::filesystem::path(path_)
              .concat(absl::LogSeverityName(severity))
              .concat(".")
              .concat(absl::FormatTime(
                  "%Y%m%d-%H%M%S", entry.timestamp(), absl::LocalTimeZone()))
              .concat(".")
              .concat(std::to_string(pid))
              .concat(extension_));
      if (!file.good()) {
        std::cerr << "Failed To create log file for:\n";
        PrintTo(entry, &std::cerr);
        std::terminate();
      }
      file << "Log file created at: " << absl::FormatTime(entry.timestamp())
           << "\n";
#if OS_WIN
      char buf[MAX_COMPUTERNAME_LENGTH + 1];
      DWORD len = MAX_COMPUTERNAME_LENGTH + 1;
      GetComputerNameA(buf, &len);
      std::string_view const computer_name = {buf, len};
#else
      struct utsname buf;
      if (uname(&buf) != 0) {
        *buf.nodename = '\0';
      }
      std::string_view const computer_name = buf.nodename;
#endif
      file << "Running on machine: " << computer_name << "\n";
      file << "Log line format: [IWEF]mmdd hh:mm:ss.μμμμμμ threadid file:line] msg\n";
    }
    // FATAL messages are sent twice, first without and then with the
    // stacktrace.  Log the message only once.
    if (entry.stacktrace().empty()) {
      file << entry.text_message_with_prefix_and_newline();
    } else {
      file << entry.stacktrace();
    }
    if (severity > buffered_level_) {
      file.flush();
    }
  }
}

void FileLogSink::Flush() {
  for (auto& file : files_) {
    file.flush();
  }
}

absl::LogSeverityAtMost FileLogSink::buffered_level() const {
  return buffered_level_;
}

void FileLogSink::set_buffered_level(absl::LogSeverityAtMost const severity) {
  buffered_level_ = severity;
}

}  // namespace _file_log_sink
}  // namespace base
}  // namespace principia